# -*- coding: utf-8 -*-

"""Semantic VASP run parser built on the structured run document."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import re
from typing import Any
import xml.etree.ElementTree as ET

import numpy as np
from numpy.typing import NDArray

from quantas.interfaces.vasp.document import (
    VaspRunDocument,
    VaspRunSource,
    _as_float,
    _text,
)
from quantas.models.computation import (
    EnergyKind,
    EnergyRecord,
    RunTermination,
    RunTerminationStatus,
)
from quantas.models.structures import CrystalStructure


FloatArray = NDArray[np.float64]

_OUTCAR_NORMAL_TERMINATION = "General timing and accounting informations for this job:"
_OUTCAR_IONIC_CONVERGENCE = (
    "reached required accuracy - stopping structural energy minimisation"
)
_OUTCAR_ELECTRONIC_CONVERGENCE = "aborting loop because EDIFF is reached"

_OUTCAR_FREE_ENERGY_RE = re.compile(
    r"free\s+energy\s+TOTEN\s*=\s*([-+0-9.EeDd]+)\s+eV",
    flags=re.IGNORECASE,
)
_OUTCAR_ZERO_T_RE = re.compile(
    r"energy\s+without\s+entropy\s*=\s*([-+0-9.EeDd]+)\s+"
    r"energy\(sigma->0\)\s*=\s*([-+0-9.EeDd]+)",
    flags=re.IGNORECASE,
)

_ENERGY_NAMES = ("e_fr_energy", "e_wo_entrp", "e_0_energy")
_ENERGY_ATOL = 5.0e-7
_VOLUME_RTOL = 1.0e-8
_VOLUME_ATOL = 5.0e-7


@dataclass(slots=True)
class VaspEnergyComponents:
    """Resolved energy quantities for one VASP ionic state.

    VASP distinguishes the electronic free energy ``F`` (``e_fr_energy``),
    energy without entropy ``E`` (``e_wo_entrp``), and the energy extrapolated
    to zero smearing ``E0`` (``e_0_energy``).  The three values remain separate
    here so a later scientific adapter can choose the appropriate quantity
    explicitly.

    Parameters
    ----------
    free_energy : EnergyRecord
        VASP ``e_fr_energy`` in eV.
    energy_without_entropy : EnergyRecord
        VASP ``e_wo_entrp`` in eV.
    sigma_zero_energy : EnergyRecord
        VASP ``e_0_energy`` in eV.
    entropy_term : float or None, optional
        VASP ``eentropy`` value, when present, in eV.
    metadata : dict, optional
        Resolution and source provenance.
    """

    free_energy: EnergyRecord
    energy_without_entropy: EnergyRecord
    sigma_zero_energy: EnergyRecord
    entropy_term: float | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Normalize the optional entropy contribution."""
        if self.entropy_term is not None:
            self.entropy_term = float(self.entropy_term)


@dataclass(slots=True)
class VaspIonicStep:
    """One normalized ionic-state record from ``vasprun.xml``.

    Parameters
    ----------
    index : int
        One-based ionic-step index in source order.
    structure : CrystalStructure
        Cell and fractional coordinates reported for the step.
    energies : VaspEnergyComponents
        Resolved VASP energy quantities for the step.
    forces : array_like or None
        Hellmann--Feynman forces in eV/angstrom, shape ``(natoms, 3)``.
    stress : array_like or None
        VASP stress tensor in kbar, shape ``(3, 3)``.
    electronic_steps : int, optional
        Number of electronic ``scstep`` records associated with the ionic step.
    metadata : dict, optional
        VASP-specific source provenance and cross-check state.

    Raises
    ------
    ValueError
        If indices or optional array shapes are invalid.
    """

    index: int
    structure: CrystalStructure
    energies: VaspEnergyComponents
    forces: FloatArray | None = None
    stress: FloatArray | None = None
    electronic_steps: int = 0
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Normalize arrays and validate dimensions."""
        self.index = int(self.index)
        self.electronic_steps = int(self.electronic_steps)
        if self.index < 1:
            raise ValueError("VASP ionic-step index must be positive")
        if self.electronic_steps < 0:
            raise ValueError("VASP electronic-step count cannot be negative")
        if self.forces is not None:
            self.forces = np.asarray(self.forces, dtype=np.float64)
            if self.forces.shape != (self.structure.natoms, 3):
                raise ValueError("VASP forces must have shape (natoms, 3)")
        if self.stress is not None:
            self.stress = np.asarray(self.stress, dtype=np.float64)
            if self.stress.shape != (3, 3):
                raise ValueError("VASP stress must have shape (3, 3)")



class VaspOutputParser:
    """Interpret one VASP run without imposing a module-specific policy.

    Parameters
    ----------
    source : str or pathlib.Path
        VASP calculation directory, ``vasprun.xml``, or sibling ``OUTCAR``.

    Notes
    -----
    VASP 5.4.4 is known to mislabel ``e_wo_entrp`` and ``e_0_energy`` in the
    outer ``calculation/energy`` block.  For calculation-style XML this parser
    therefore resolves energy differences from the final electronic ``scstep``
    and transfers only the outer free-energy shift.  This also retains an
    additive correction present only in the outer ionic-state energy.
    """

    def __init__(self, source: str | Path) -> None:
        self.document = VaspRunDocument(source)

    @property
    def source(self) -> VaspRunSource:
        """Return the resolved VASP run source."""
        return self.document.source

    @property
    def version(self) -> str | None:
        """Return the reported VASP version, if available."""
        value = self.document.generator().get("version")
        return value or None

    def termination(self) -> RunTermination:
        """Return the observed VASP process-termination state.

        Returns
        -------
        RunTermination
            ``NORMAL`` when OUTCAR has the final timing section or, without an
            OUTCAR, a complete XML document contains ``finalpos``.  Otherwise
            the state is ``UNKNOWN`` because XML well-formedness alone does not
            establish why VASP stopped.
        """
        text = self.document.outcar_text
        if text is not None and _OUTCAR_NORMAL_TERMINATION in text:
            return RunTermination(
                status=RunTerminationStatus.NORMAL,
                metadata={"source": "OUTCAR", "marker": _OUTCAR_NORMAL_TERMINATION},
            )
        if text is None and self.document.root.find("structure[@name='finalpos']") is not None:
            return RunTermination(
                status=RunTerminationStatus.NORMAL,
                metadata={"source": "vasprun.xml", "marker": "finalpos"},
            )
        return RunTermination(
            status=RunTerminationStatus.UNKNOWN,
            metadata={"source": "OUTCAR" if text is not None else "vasprun.xml"},
        )

    def initial_structure(self) -> CrystalStructure:
        """Return the initial VASP structure in Quantas convention.

        Returns
        -------
        CrystalStructure
            Initial lattice, fractional coordinates, and atom ordering.
        """
        return _structure(
            self.document,
            self.document.structure_node("initialpos"),
            label="VASP initial structure",
            source_name="initialpos",
        )

    def final_structure(self) -> CrystalStructure:
        """Return the final VASP structure in Quantas convention.

        Returns
        -------
        CrystalStructure
            Final lattice, fractional coordinates, and atom ordering.
        """
        return _structure(
            self.document,
            self.document.structure_node("finalpos"),
            label="VASP final structure",
            source_name="finalpos",
        )

    def ionic_steps(self) -> tuple[VaspIonicStep, ...]:
        """Return all ionic states with structures, energies, forces, and stress.

        Returns
        -------
        tuple of VaspIonicStep
            Ionic states in VASP source order.

        Raises
        ------
        ValueError
            If a required structure or energy is absent, if XML and OUTCAR
            energy records disagree when they can be paired unambiguously, or
            if numerical array shapes are inconsistent.
        """
        containers = self.document.ionic_step_nodes()
        outcar_energies = _outcar_energy_components(self.document)
        validate_outcar = len(outcar_energies) == len(containers) and bool(containers)
        steps: list[VaspIonicStep] = []
        for index, container in enumerate(containers, start=1):
            structure_node = container.find("structure")
            if structure_node is None:
                raise ValueError(f"VASP ionic step {index} has no structure")
            structure = _structure(
                self.document,
                structure_node,
                label=f"VASP ionic step {index}",
                source_name="ionic_step",
            )
            energies = _resolve_energies(
                self.document, self.version, container, step_index=index
            )
            forces = _varray(container, "forces")
            stress = _varray(container, "stress")
            metadata: dict[str, Any] = {
                "source_layout": (
                    "calculation" if container.tag == "calculation" else "flat"
                ),
                "electronic_converged": _electronic_converged(
                    self.document, index, containers
                ),
                "ionic_converged": self.ionic_converged(),
            }
            if validate_outcar:
                outcar = outcar_energies[index - 1]
                _validate_energy_components(energies, outcar, index=index)
                metadata["outcar_energy_crosscheck"] = "matched"
            elif self.document.outcar_text is not None:
                metadata["outcar_energy_crosscheck"] = "not_uniquely_pairable"
            steps.append(
                VaspIonicStep(
                    index=index,
                    structure=structure,
                    energies=energies,
                    forces=forces,
                    stress=stress,
                    electronic_steps=len(container.findall("scstep")),
                    metadata=metadata,
                )
            )
        return tuple(steps)

    def ionic_converged(self) -> bool | None:
        """Return whether VASP explicitly reported ionic convergence.

        Returns
        -------
        bool or None
            ``True`` when OUTCAR contains VASP's required-accuracy marker;
            ``None`` when no OUTCAR is available or the marker is absent.  The
            parser does not infer a failed optimization merely from its absence.
        """
        text = self.document.outcar_text
        if text is None:
            return None
        if _OUTCAR_IONIC_CONVERGENCE in text:
            return True
        return None


def _structure(
    document: VaspRunDocument,
    node: ET.Element,
    *,
    label: str,
    source_name: str,
) -> CrystalStructure:
    """Convert one VASP structure node to ``CrystalStructure``."""
    basis = _varray(node.find("crystal"), "basis")
    positions = _varray(node, "positions")
    if basis is None or basis.shape != (3, 3):
        raise ValueError(f"VASP {source_name} basis must have shape (3, 3)")
    if positions is None or positions.ndim != 2 or positions.shape[1] != 3:
        raise ValueError(
            f"VASP {source_name} fractional positions must have shape (natoms, 3)"
        )
    numbers = document.atomic_numbers()
    if positions.shape[0] != numbers.size:
        raise ValueError(
            f"VASP {source_name} positions contain {positions.shape[0]} atoms "
            f"but atominfo contains {numbers.size}"
        )

    reported_volume = _named_float(node.find("crystal"), "volume")
    determinant_volume = abs(float(np.linalg.det(basis)))
    if reported_volume is not None and not np.isclose(
        determinant_volume,
        reported_volume,
        rtol=_VOLUME_RTOL,
        atol=_VOLUME_ATOL,
    ):
        raise ValueError(
            f"VASP {source_name} lattice volume {determinant_volume:.10g} Å^3 "
            f"disagrees with reported volume {reported_volume:.10g} Å^3"
        )
    metadata: dict[str, Any] = {
        "interface": "vasp",
        "vasp_structure_name": source_name,
        "coordinate_convention": "fractional",
    }
    if reported_volume is not None:
        metadata["reported_volume_angstrom3"] = float(reported_volume)
    return CrystalStructure(
        lattice=basis,
        fractional_positions=positions,
        atomic_numbers=numbers,
        label=label,
        metadata=metadata,
    )

def _resolve_energies(
    document: VaspRunDocument,
    version: str | None,
    container: ET.Element,
    *,
    step_index: int,
) -> VaspEnergyComponents:
    """Resolve one ionic state's VASP F, E, and E0 values."""
    scsteps = container.findall("scstep")
    inner_node = scsteps[-1].find("energy") if scsteps else None
    outer_node = container.find("energy")
    inner = _energy_values(inner_node)
    outer = _energy_values(outer_node)
    if not inner and not outer:
        raise ValueError(f"VASP ionic step {step_index} has no energy record")

    metadata: dict[str, Any] = {
        "ionic_step": int(step_index),
        "source_layout": (
            "calculation" if container.tag == "calculation" else "flat"
        ),
    }

    entropy_term: float | None = None
    if inner:
        entropy_term = inner.get("eentropy")
        metadata["raw_last_scstep_energy_eV"] = {
            key: float(value)
            for key, value in inner.items()
            if key in (*_ENERGY_NAMES, "eentropy")
        }

    if inner and all(name in inner for name in _ENERGY_NAMES):
        resolved = {name: float(inner[name]) for name in _ENERGY_NAMES}
        resolution = "last_scstep"
        shift = 0.0
        if outer and "e_fr_energy" in outer:
            shift = float(outer["e_fr_energy"]) - resolved["e_fr_energy"]
            for name in _ENERGY_NAMES:
                resolved[name] += shift
            resolution = "last_scstep_plus_outer_free_energy_shift"
            metadata["outer_free_energy_shift_eV"] = float(shift)
            metadata["outer_energy_tags_consistent"] = _outer_energy_consistent(
                outer,
                resolved,
            )
            metadata["known_vasp5_outer_energy_tag_bug"] = (
                _matches_vasp5_outer_energy_bug(outer, inner, resolved)
            )
            metadata["raw_outer_energy_eV"] = {
                key: float(value)
                for key, value in outer.items()
                if key in (*_ENERGY_NAMES, "eentropy")
            }
    elif outer and all(name in outer for name in _ENERGY_NAMES):
        resolved = {name: float(outer[name]) for name in _ENERGY_NAMES}
        resolution = "outer_ionic_energy"
        entropy_term = outer.get("eentropy")
    else:
        missing = [
            name
            for name in _ENERGY_NAMES
            if name not in (inner if inner else outer)
        ]
        raise ValueError(
            f"VASP ionic step {step_index} is missing required energy values: "
            + ", ".join(missing)
        )

    metadata["resolution"] = resolution
    common = {
        "ionic_step": int(step_index),
        "resolution": resolution,
        "source_file": str(document.source.vasprun_xml),
        "vasp_version": version,
    }
    return VaspEnergyComponents(
        free_energy=EnergyRecord(
            value=resolved["e_fr_energy"],
            unit="eV",
            kind=EnergyKind.DFT,
            metadata={
                **common,
                "vasp_tag": "e_fr_energy",
                "vasp_semantics": "electronic_free_energy",
            },
        ),
        energy_without_entropy=EnergyRecord(
            value=resolved["e_wo_entrp"],
            unit="eV",
            kind=EnergyKind.DFT,
            metadata={
                **common,
                "vasp_tag": "e_wo_entrp",
                "vasp_semantics": "energy_without_entropy",
            },
        ),
        sigma_zero_energy=EnergyRecord(
            value=resolved["e_0_energy"],
            unit="eV",
            kind=EnergyKind.DFT,
            metadata={
                **common,
                "vasp_tag": "e_0_energy",
                "vasp_semantics": "sigma_to_zero_energy",
            },
        ),
        entropy_term=entropy_term,
        metadata=metadata,
    )

def _outcar_energy_components(
    document: VaspRunDocument,
) -> tuple[VaspEnergyComponents, ...]:
    """Return paired OUTCAR energy triplets in source order."""
    text = document.outcar_text
    if text is None:
        return ()
    marker_re = re.compile(re.escape(_OUTCAR_ELECTRONIC_CONVERGENCE))
    markers = list(marker_re.finditer(text))
    pairs: list[tuple[float, tuple[float, float]]] = []
    for position, marker in enumerate(markers):
        stop = markers[position + 1].start() if position + 1 < len(markers) else len(text)
        free_match = _OUTCAR_FREE_ENERGY_RE.search(text, marker.end(), stop)
        zero_match = _OUTCAR_ZERO_T_RE.search(text, marker.end(), stop)
        if free_match is None or zero_match is None:
            return ()
        pairs.append(
            (
                _as_float(free_match.group(1)),
                (_as_float(zero_match.group(1)), _as_float(zero_match.group(2))),
            )
        )
    if not pairs:
        return ()
    records: list[VaspEnergyComponents] = []
    for index, (free, (without_entropy, sigma_zero)) in enumerate(
        pairs,
        start=1,
    ):
        common = {"source_file": str(document.source.outcar), "ionic_step": index}
        records.append(
            VaspEnergyComponents(
                free_energy=EnergyRecord(
                    free,
                    "eV",
                    EnergyKind.DFT,
                    {**common, "vasp_tag": "TOTEN"},
                ),
                energy_without_entropy=EnergyRecord(
                    without_entropy,
                    "eV",
                    EnergyKind.DFT,
                    {**common, "vasp_tag": "energy_without_entropy"},
                ),
                sigma_zero_energy=EnergyRecord(
                    sigma_zero,
                    "eV",
                    EnergyKind.DFT,
                    {**common, "vasp_tag": "energy(sigma->0)"},
                ),
                metadata={"resolution": "OUTCAR"},
            )
        )
    return tuple(records)

def _electronic_converged(
    document: VaspRunDocument,
    step_index: int,
    containers: tuple[ET.Element, ...],
) -> bool | None:
    """Return an explicit electronic-convergence observation when pairable."""
    text = document.outcar_text
    if text is None:
        return None
    count = text.count(_OUTCAR_ELECTRONIC_CONVERGENCE)
    if count == len(containers):
        return True
    if step_index <= count and count > 0:
        return None
    return None


def _varray(parent: ET.Element | None, name: str) -> FloatArray | None:
    """Return one VASP ``varray`` as a float64 matrix."""
    if parent is None:
        return None
    node = parent.find(f"varray[@name='{name}']")
    if node is None:
        return None
    rows: list[list[float]] = []
    for vector in node.findall("v"):
        values = [_as_float(item) for item in _text(vector).split()]
        rows.append(values)
    if not rows:
        return np.empty((0, 0), dtype=np.float64)
    width = len(rows[0])
    if width == 0 or any(len(row) != width for row in rows):
        raise ValueError(f"VASP varray {name!r} has inconsistent row widths")
    return np.asarray(rows, dtype=np.float64)


def _named_float(parent: ET.Element | None, name: str) -> float | None:
    """Return one named scalar float below ``parent`` when present."""
    if parent is None:
        return None
    node = parent.find(f"i[@name='{name}']")
    if node is None:
        return None
    return _as_float(_text(node))


def _energy_values(node: ET.Element | None) -> dict[str, float]:
    """Return named scalar energies from one VASP ``energy`` element."""
    if node is None:
        return {}
    values: dict[str, float] = {}
    for item in node.findall("i"):
        name = item.get("name")
        if not name:
            continue
        try:
            values[name] = _as_float(_text(item))
        except ValueError:
            continue
    return values


def _outer_energy_consistent(
    outer: dict[str, float],
    resolved: dict[str, float],
) -> bool | None:
    """Return whether complete outer energy tags match reconstructed values."""
    if not all(name in outer for name in _ENERGY_NAMES):
        return None
    return all(
        np.isclose(float(outer[name]), resolved[name], rtol=0.0, atol=_ENERGY_ATOL)
        for name in _ENERGY_NAMES
    )


def _matches_vasp5_outer_energy_bug(
    outer: dict[str, float],
    inner: dict[str, float],
    resolved: dict[str, float],
) -> bool:
    """Identify the documented VASP 5.4.4 outer-energy tag permutation."""
    if not all(name in outer for name in _ENERGY_NAMES):
        return False
    entropy = inner.get("eentropy")
    if entropy is None:
        return False
    return bool(
        np.isclose(
            outer["e_wo_entrp"],
            resolved["e_0_energy"],
            rtol=0.0,
            atol=_ENERGY_ATOL,
        )
        and np.isclose(
            outer["e_0_energy"],
            entropy,
            rtol=0.0,
            atol=_ENERGY_ATOL,
        )
    )


def _validate_energy_components(
    xml: VaspEnergyComponents,
    outcar: VaspEnergyComponents,
    *,
    index: int,
) -> None:
    """Require XML-resolved and OUTCAR VASP energies to agree."""
    pairs = (
        ("free energy", xml.free_energy.value, outcar.free_energy.value),
        (
            "energy without entropy",
            xml.energy_without_entropy.value,
            outcar.energy_without_entropy.value,
        ),
        ("sigma->0 energy", xml.sigma_zero_energy.value, outcar.sigma_zero_energy.value),
    )
    for label, xml_value, outcar_value in pairs:
        if not np.isclose(xml_value, outcar_value, rtol=0.0, atol=_ENERGY_ATOL):
            raise ValueError(
                f"VASP ionic step {index} {label} disagrees between "
                f"vasprun.xml ({xml_value:.10g} eV) and OUTCAR "
                f"({outcar_value:.10g} eV)"
            )
