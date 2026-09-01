# -*- coding: utf-8 -*-

"""Parse CRYSTAL phonon eigenvectors into backend-neutral mode data."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np

from quantas.core.chemistry import atomic_mass
from quantas.interfaces.crystal import markers, patterns
from quantas.interfaces.crystal.output import CrystalOutputDocument
from quantas.models.phonons import PhononModeData


@dataclass(slots=True)
class _ModeComponent:
    """One real-valued CRYSTAL eigenvector component before normalization."""

    frequencies: np.ndarray
    vectors: np.ndarray
    atom_symbols: tuple[str, ...]


def _as_float(value: str) -> float:
    """Convert a CRYSTAL/Fortran numeric literal to Python ``float``."""
    return float(value.replace("D", "E").replace("d", "e"))


def _numeric_values(text: str) -> np.ndarray:
    """Return all numeric fields from one structured CRYSTAL record."""
    return np.asarray(
        [_as_float(value) for value in patterns.FLOAT_RE.findall(text)],
        dtype=np.float64,
    )


class CrystalPhononModeParser:
    """Parse phonon frequencies and displacement eigenvectors from CRYSTAL.

    CRYSTAL prints eigenvectors in blocks containing at most six modes.  The
    final block may contain fewer columns.  At Gamma the vectors are real and
    are headed by ``NORMAL MODES NORMALIZED TO CLASSICAL AMPLITUDES``.  At a
    general q-point CRYSTAL prints separate in-phase and anti-phase components,
    which are combined as the real and imaginary parts of a complex vector.

    The source vectors are normalized to classical displacement amplitudes in
    bohr.  Quantas removes that amplitude and returns unit-norm mass-weighted
    eigenvectors.  This is the representation required for mode-overlap
    comparisons and is backend-neutral.

    Parameters
    ----------
    source : CrystalOutputDocument or str or Path or iterable of str
        CRYSTAL output document or source used to construct one.
    """

    def __init__(
        self,
        source: CrystalOutputDocument | str | Path | Iterable[str],
    ) -> None:
        """Initialize the parser from one CRYSTAL output source."""
        if isinstance(source, CrystalOutputDocument):
            self.document = source
        else:
            self.document = CrystalOutputDocument(source)

    @property
    def lines(self) -> list[str]:
        """Return source lines in original order."""
        return self.document.lines

    def parse(self, nmodes: int) -> PhononModeData | None:
        """Parse all printed phonon eigenvectors.

        Parameters
        ----------
        nmodes : int
            Expected number of vibrational modes at each q-point.  Quantas
            currently expects the complete ``3 * natoms`` mode set.

        Returns
        -------
        PhononModeData or None
            Parsed frequencies and unit-norm mass-weighted eigenvectors, or
            ``None`` when the output does not contain printed eigenvectors.

        Raises
        ------
        ValueError
            If CRYSTAL prints an incomplete or internally inconsistent mode
            block, q-point sequence, atom sequence, or phase decomposition.
        """
        nmodes = int(nmodes)
        if nmodes <= 0 or nmodes % 3 != 0:
            raise ValueError("nmodes must be a positive multiple of three")
        natoms = nmodes // 3

        sections = self._mode_sections()
        if not sections:
            return None

        frequencies: list[np.ndarray] = []
        classical_vectors: list[np.ndarray] = []
        atom_symbols: tuple[str, ...] | None = None

        for q_index, section in enumerate(sections, start=1):
            component = self._parse_qpoint_section(section, nmodes, natoms)
            if component is None:
                raise ValueError(
                    f"CRYSTAL q-point {q_index} has no readable eigenvectors"
                )
            if atom_symbols is None:
                atom_symbols = component.atom_symbols
            elif component.atom_symbols != atom_symbols:
                raise ValueError(
                    "CRYSTAL phonon eigenvector atom ordering changes between "
                    "q-points"
                )
            frequencies.append(component.frequencies)
            classical_vectors.append(component.vectors)

        if atom_symbols is None:
            return None

        source_vectors = np.asarray(classical_vectors, dtype=np.complex128)
        normalized = self._remove_classical_amplitude(source_vectors, atom_symbols)
        return PhononModeData(
            frequencies=np.asarray(frequencies, dtype=np.float64),
            eigenvectors=normalized,
            atom_symbols=atom_symbols,
            frequency_unit="cm^-1",
            eigenvector_normalization="mass-weighted-unit",
            metadata={
                "interface": "crystal",
                "source_normalization": "classical-amplitude",
                "source_length_unit": "bohr",
                "mass_source": "standard-atomic-mass",
            },
        )

    def _mode_sections(self) -> list[list[str]]:
        """Return one source-line section per printed phonon q-point."""
        qpoint_starts = [
            index
            for index, line in enumerate(self.lines)
            if markers.DISPERSION_QPOINT in line
        ]
        if not qpoint_starts:
            mode_starts = [
                index
                for index, line in enumerate(self.lines)
                if self._has_mode_header(line)
            ]
            if not mode_starts:
                return []
            return [
                self.lines[start : mode_starts[position + 1]]
                if position + 1 < len(mode_starts)
                else self.lines[start:]
                for position, start in enumerate(mode_starts)
            ]

        indices: list[int] = []
        sections: list[list[str]] = []
        for position, start in enumerate(qpoint_starts):
            match = patterns.DISPERSION_QPOINT_INDEX_RE.search(self.lines[start])
            if match is None:
                raise ValueError("unable to parse CRYSTAL dispersion q-point index")
            indices.append(int(match.group("index")))
            end = (
                qpoint_starts[position + 1]
                if position + 1 < len(qpoint_starts)
                else len(self.lines)
            )
            sections.append(self.lines[start:end])

        expected = list(range(1, len(indices) + 1))
        if indices != expected:
            raise ValueError("CRYSTAL phonon q-point indices are not consecutive")
        return sections

    @staticmethod
    def _has_mode_header(line: str) -> bool:
        """Return whether one line begins a CRYSTAL eigenvector section."""
        return (
            markers.NORMAL_MODES_CLASSICAL in line
            or markers.MODES_IN_PHASE in line
        )

    def _parse_qpoint_section(
        self,
        lines: Sequence[str],
        nmodes: int,
        natoms: int,
    ) -> _ModeComponent | None:
        """Parse one Gamma or complex-q eigenvector section."""
        normal_index = self._find_marker(lines, markers.NORMAL_MODES_CLASSICAL)
        if normal_index is not None:
            real = self._parse_component(lines[normal_index + 1 :], nmodes, natoms)
            return _ModeComponent(
                frequencies=real.frequencies,
                vectors=real.vectors.astype(np.complex128),
                atom_symbols=real.atom_symbols,
            )

        phase_index = self._find_marker(lines, markers.MODES_IN_PHASE)
        if phase_index is None:
            return None
        anti_index = self._find_marker(lines, markers.MODES_IN_ANTI_PHASE)
        if anti_index is None or anti_index <= phase_index:
            raise ValueError(
                "CRYSTAL complex phonon modes are missing the anti-phase block"
            )

        real = self._parse_component(
            lines[phase_index + 1 : anti_index],
            nmodes,
            natoms,
        )
        imaginary = self._parse_component(
            lines[anti_index + 1 :],
            nmodes,
            natoms,
        )
        if real.atom_symbols != imaginary.atom_symbols:
            raise ValueError(
                "CRYSTAL in-phase and anti-phase blocks use different atoms"
            )
        if not np.allclose(
            real.frequencies,
            imaginary.frequencies,
            rtol=0.0,
            atol=1.0e-8,
        ):
            raise ValueError(
                "CRYSTAL in-phase and anti-phase blocks use different frequencies"
            )
        return _ModeComponent(
            frequencies=real.frequencies,
            vectors=real.vectors + 1j * imaginary.vectors,
            atom_symbols=real.atom_symbols,
        )

    @staticmethod
    def _find_marker(lines: Sequence[str], marker: str) -> int | None:
        """Return the first index containing a literal CRYSTAL marker."""
        return next(
            (index for index, line in enumerate(lines) if marker in line),
            None,
        )

    def _parse_component(
        self,
        lines: Sequence[str],
        nmodes: int,
        natoms: int,
    ) -> _ModeComponent:
        """Parse one real-valued CRYSTAL mode component."""
        frequencies: list[float] = []
        vector_blocks: list[np.ndarray] = []
        atom_symbols: tuple[str, ...] | None = None
        index = 0

        while index < len(lines) and len(frequencies) < nmodes:
            match = patterns.PHONON_VECTOR_FREQUENCY_RE.search(lines[index])
            if match is None:
                index += 1
                continue

            block_frequencies = _numeric_values(match.group("values"))
            block_modes = int(block_frequencies.size)
            if block_modes <= 0:
                raise ValueError("CRYSTAL phonon eigenvector block has no modes")
            if len(frequencies) + block_modes > nmodes:
                raise ValueError(
                    "CRYSTAL phonon eigenvector blocks contain too many modes"
                )

            index += 1
            block = np.zeros((block_modes, natoms, 3), dtype=np.float64)
            symbols: list[str] = []
            for atom in range(1, natoms + 1):
                index = self._next_nonempty(lines, index)
                if index >= len(lines):
                    raise ValueError("CRYSTAL phonon eigenvector block is truncated")
                atom_match = patterns.PHONON_VECTOR_ATOM_RE.match(lines[index])
                if atom_match is None:
                    raise ValueError(
                        "unable to parse CRYSTAL phonon X-component row: "
                        f"{lines[index]!r}"
                    )
                if int(atom_match.group("atom")) != atom:
                    raise ValueError(
                        "CRYSTAL phonon eigenvector atom indices are out of order"
                    )
                if atom_match.group("axis").upper() != "X":
                    raise ValueError("CRYSTAL phonon atom row must begin with X")
                symbols.append(atom_match.group("symbol").capitalize())
                block[:, atom - 1, 0] = self._component_values(
                    atom_match.group("values"),
                    block_modes,
                )
                index += 1

                for coordinate, axis in enumerate(("Y", "Z"), start=1):
                    index = self._next_nonempty(lines, index)
                    if index >= len(lines):
                        raise ValueError(
                            "CRYSTAL phonon eigenvector block is truncated"
                        )
                    continuation = patterns.PHONON_VECTOR_CONTINUATION_RE.match(
                        lines[index]
                    )
                    if (
                        continuation is None
                        or continuation.group("axis").upper() != axis
                    ):
                        raise ValueError(
                            f"CRYSTAL phonon atom {atom} is missing its {axis} row"
                        )
                    block[:, atom - 1, coordinate] = self._component_values(
                        continuation.group("values"),
                        block_modes,
                    )
                    index += 1

            block_symbols = tuple(symbols)
            if atom_symbols is None:
                atom_symbols = block_symbols
            elif block_symbols != atom_symbols:
                raise ValueError(
                    "CRYSTAL phonon eigenvector atom ordering changes between "
                    "frequency blocks"
                )
            frequencies.extend(block_frequencies.tolist())
            vector_blocks.append(block)

        if len(frequencies) != nmodes or not vector_blocks or atom_symbols is None:
            raise ValueError(
                f"CRYSTAL phonon eigenvectors contain {len(frequencies)} of "
                f"{nmodes} expected modes"
            )
        return _ModeComponent(
            frequencies=np.asarray(frequencies, dtype=np.float64),
            vectors=np.concatenate(vector_blocks, axis=0),
            atom_symbols=atom_symbols,
        )

    @staticmethod
    def _next_nonempty(lines: Sequence[str], index: int) -> int:
        """Advance to the next non-empty source line."""
        while index < len(lines) and not lines[index].strip():
            index += 1
        return index

    @staticmethod
    def _component_values(text: str, expected: int) -> np.ndarray:
        """Parse and validate one CRYSTAL eigenvector coordinate row."""
        values = _numeric_values(text)
        if values.size != expected:
            raise ValueError(
                "CRYSTAL phonon eigenvector row contains "
                f"{values.size} values; expected {expected}"
            )
        return values

    @staticmethod
    def _remove_classical_amplitude(
        eigenvectors: np.ndarray,
        atom_symbols: tuple[str, ...],
    ) -> np.ndarray:
        """Return unit-norm mass-weighted eigenvectors.

        CRYSTAL prints mass-unweighted displacements multiplied by a scalar
        classical amplitude for each mode.  Multiplication by ``sqrt(mass)``
        restores the mass weighting.  The classical amplitude and source length
        unit are both scalar factors for a given mode and cancel when each mode
        is normalized to unit norm.  The result is therefore equivalent to
        removing the classical amplitude before normalization, as done by
        CRYSTALpytools.
        """
        masses = np.asarray(
            [atomic_mass(symbol) for symbol in atom_symbols],
            dtype=np.float64,
        )
        weighted = np.asarray(eigenvectors, dtype=np.complex128) * np.sqrt(masses)[
            np.newaxis, np.newaxis, :, np.newaxis
        ]
        flattened = weighted.reshape(weighted.shape[0], weighted.shape[1], -1)
        norms = np.linalg.norm(flattened, axis=2)
        if np.any(~np.isfinite(norms)) or np.any(norms <= np.finfo(np.float64).eps):
            raise ValueError("CRYSTAL phonon eigenvector has zero or invalid norm")
        return weighted / norms[:, :, np.newaxis, np.newaxis]
