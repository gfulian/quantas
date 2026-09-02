# -*- coding: utf-8 -*-

"""
Input generator for Quantas HA/QHA YAML files.

This module creates Quantas harmonic and quasi-harmonic YAML input files from
quantum-mechanical phonon outputs. The implementation is frontend-neutral: it
contains no Click prompts, printing logic, or command-line state. Command-line
and graphical interfaces are expected to collect user options and call the API
provided here.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, ClassVar

import numpy as np
from numpy.typing import ArrayLike
import yaml

from quantas.core.events import Event, EventLevel, NullObserver, Observer
from quantas.core.numerics.phonon_tracking import (
    PhononModeTrackingResult,
    PhononModeTrackingStep,
    track_phonon_modes,
)
from quantas.models import ReportTable
from quantas.models.phonons import default_phonon_input_units
from quantas.models.structures import StructureVolumeSeries

_CrystalPhononReader: type[Any] | None
_CrystalQHAReader: type[Any] | None
_PhonopyReader: type[Any] | None

try:  # pragma: no cover - import availability depends on installed interfaces
    from quantas.interfaces.crystal.phonons import (
        CrystalPhononReader as _ImportedCrystalPhononReader,
    )

    _CrystalPhononReader = _ImportedCrystalPhononReader
except Exception:  # pragma: no cover
    _CrystalPhononReader = None

try:  # pragma: no cover - import availability depends on installed interfaces
    from quantas.interfaces.crystal.qha import (
        CrystalQHAReader as _ImportedCrystalQHAReader,
    )

    _CrystalQHAReader = _ImportedCrystalQHAReader
except Exception:  # pragma: no cover
    _CrystalQHAReader = None

try:  # pragma: no cover - import availability depends on installed interfaces
    from quantas.interfaces.phonopy.phonons import (
        PhonopyReader as _ImportedPhonopyReader,
    )

    _PhonopyReader = _ImportedPhonopyReader
except Exception:  # pragma: no cover
    _PhonopyReader = None


class HAInputCreator:
    """
    Create Quantas HA/QHA YAML input files from QM phonon outputs.

    Parameters
    ----------
    interface : str or None, optional
        Input interface name. Supported values are ``"crystal"``,
        ``"crystal-qha"``, and ``"phonopy"``.
    interface_filter : dict or None, optional
        Optional mapping used to override the default interface readers. This
        is mainly useful for tests and for external applications that provide
        compatible reader classes.
    observer : Observer or None, optional
        Frontend-neutral observer receiving input-generation events.

    Attributes
    ----------
    interface_flag : str
        Selected interface name.
    phondata : list or object
        Parsed phonon data. Multiple-file interfaces store a list of readers;
        the ``"crystal-qha"`` interface stores a single reader.
    files : list of pathlib.Path
        Files read by the creator.
    """

    interface_filter: ClassVar[dict[str, type[Any] | None]] = {
        "crystal": _CrystalPhononReader,
        "crystal-qha": _CrystalQHAReader,
        "phonopy": _PhonopyReader,
    }

    def __init__(
        self,
        interface: str | None = None,
        interface_filter: dict[str, type[Any] | None] | None = None,
        observer: Observer | None = None,
    ) -> None:
        """
        Initialize the input creator.

        Parameters
        ----------
        interface : str or None, optional
            Input interface name.
        interface_filter : dict or None, optional
            Optional reader mapping overriding the default interfaces.
        observer : Observer or None, optional
            Frontend-neutral observer receiving input-generation events.
        """
        self.interface_flag = interface or "crystal"
        self._interface_filter = (
            dict(self.interface_filter)
            if interface_filter is None
            else dict(interface_filter)
        )
        self.phondata: list[Any] | Any = []
        self.files: list[Path] = []
        self.observer = observer if observer is not None else NullObserver()
        self.last_tracking_result: PhononModeTrackingResult | None = None

    def emit(
        self,
        message: str,
        *,
        level: EventLevel = EventLevel.INFO,
        data: dict[str, Any] | None = None,
    ) -> None:
        """Emit one frontend-neutral input-generation event.

        Parameters
        ----------
        message : str
            Human-readable event message.
        level : EventLevel, optional
            Event severity or structured-result level.
        data : dict or None, optional
            Optional structured payload for frontend observers.
        """
        self.observer(Event(message=message, level=level, data=data or {}))

    def read(
        self,
        filename: str | Path,
        is_list: bool = False,
        reference: int = 0,
        use_symm: bool = False,
    ) -> tuple[bool, str | None]:
        """
        Read one or more QM output files through the selected interface.

        Parameters
        ----------
        filename : str or pathlib.Path
            QM output file or text file containing a list of QM output files.
        is_list : bool, optional
            If ``True``, ``filename`` is interpreted as a file list.
        reference : int, optional
            Reference file index for multiple-file input generation.
        use_symm : bool, optional
            Kept for compatibility with the historical API. It is not used by
            the current HA input generator.

        Returns
        -------
        tuple of bool and str or None
            ``(True, None)`` if all files were read successfully. Otherwise,
            ``False`` and a human-readable error message are returned.
        """
        del use_symm

        if self.interface_flag == "crystal-qha" and is_list:
            return False, "Only one CRYSTAL output file for --crystal-qha interface"

        if self.interface_flag not in self._interface_filter:
            return False, f"Unsupported interface: {self.interface_flag}"

        reader_cls = self._interface_filter[self.interface_flag]
        if reader_cls is None:
            return False, f"Interface not available: {self.interface_flag}"

        try:
            files = self.get_files_from_list(filename) if is_list else [Path(filename)]
        except OSError as exc:
            return False, f"Unable to read file list: {exc}"

        if not files:
            return False, "No input files provided"

        if reference < 0 or reference >= len(files):
            return False, "Invalid reference provided"

        self.emit(
            f"Reading {len(files)} phonon output file(s) with the "
            f"{self.interface_flag} interface",
            data={
                "kind": "phonon_input_sources",
                "interface": self.interface_flag,
                "source_count": len(files),
                "reference_index": int(reference),
            },
        )
        self.phondata = [] if self.interface_flag != "crystal-qha" else None
        self.files = []

        for file in files:
            file = Path(file)
            if not file.exists():
                return False, f"File {file} does not exists"

            try:
                data = self._read_interface_file(reader_cls, file)
            except UnicodeDecodeError:
                return False, f"File '{file}' is in binary format"
            except Exception as exc:  # noqa: BLE001 - frontend-neutral error path
                return False, str(exc)

            if getattr(data, "completed", False) is not True:
                return False, getattr(data, "error", None) or f"Unable to read {file}"

            if self.interface_flag == "crystal-qha":
                self.phondata = data
            else:
                if not isinstance(self.phondata, list):
                    raise RuntimeError("multiple-file interface state is not a list")
                self.phondata.append(data)
            self.files.append(file)
            self.emit(
                f"Parsed phonon source {file}",
                level=EventLevel.DEBUG,
                data={
                    "kind": "phonon_source_parsed",
                    "source": str(file),
                    "source_index": len(self.files) - 1,
                },
            )

        self.emit(
            "Phonon source parsing completed",
            level=EventLevel.RESULT,
            data={
                "kind": "phonon_sources_parsed",
                "interface": self.interface_flag,
                "source_count": len(self.files),
                "reference_index": int(reference),
            },
        )
        return True, None

    @staticmethod
    def get_files_from_list(file_list: str | Path) -> list[Path]:
        """
        Read a list of filenames from a text file.

        Parameters
        ----------
        file_list : str or pathlib.Path
            Text file containing one filename per line. Empty lines are ignored.

        Returns
        -------
        list of pathlib.Path
            Parsed filenames.
        """
        path = Path(file_list)
        files: list[Path] = []
        with path.open("r", encoding="utf-8") as stream:
            for line in stream:
                value = line.strip()
                if not value or value.startswith("#"):
                    continue
                item = Path(value)
                if not item.is_absolute():
                    item = path.parent / item
                files.append(item)
        return files

    @staticmethod
    def _read_interface_file(reader_cls: type[Any], filename: Path) -> Any:
        """
        Instantiate an interface reader and ensure that ``load`` is called.

        Parameters
        ----------
        reader_cls : type
            Reader class implementing the Quantas interface reader attributes.
        filename : pathlib.Path
            File to read.

        Returns
        -------
        object
            Loaded interface reader instance.
        """
        try:
            reader = reader_cls(filename)
        except TypeError:
            reader = reader_cls()

        if getattr(reader, "completed", False) is not True and hasattr(reader, "load"):
            loaded = reader.load(filename)
            if loaded is not None:
                reader = loaded
        return reader

    def to_dict(
        self,
        jobname: str,
        reference: int = 0,
        formula_units: int = 1,
    ) -> dict[str, Any]:
        """
        Build the Quantas HA/QHA YAML data structure.

        Parameters
        ----------
        jobname : str
            Short description stored in the ``job`` field.
        reference : int, optional
            Reference dataset used for lattice, supercell, q-point and mode
            metadata in multiple-file input generation.
        formula_units : int, optional
            Number of chemical formula units in the thermodynamic normalization cell.

        Returns
        -------
        dict
            Ordered mapping ready to be serialized as YAML.

        Raises
        ------
        ValueError
            If no interface data have been loaded or the reference index is
            invalid.
        """
        if formula_units <= 0:
            raise ValueError("formula_units must be positive")
        if self.interface_flag == "crystal-qha":
            if self.phondata is None:
                raise ValueError("no phonon data have been loaded")
            data = self._dict_from_single_reader(
                self.phondata, jobname, formula_units
            )
            data["provenance"] = _provenance_dict(
                self.interface_flag, self.files, reference_index=0
            )
            self.last_tracking_result = None
            self.emit(
                "Using phonon-mode continuity reported by the source QHA workflow",
                level=EventLevel.RESULT,
                data={
                    "kind": "mode_continuity_source",
                    "status": data.get("mode_continuity", "unknown"),
                    "metadata": dict(data.get("mode_continuity_metadata", {})),
                },
            )
            self._emit_input_summary(data)
            return data

        if not isinstance(self.phondata, list) or not self.phondata:
            raise ValueError("no phonon data have been loaded")
        if reference < 0 or reference >= len(self.phondata):
            raise ValueError("invalid reference provided")
        data = self._dict_from_multiple_readers(
            self.phondata, jobname, reference, formula_units
        )
        data["provenance"] = _provenance_dict(
            self.interface_flag, self.files, reference_index=reference
        )
        self._emit_input_summary(data)
        return data

    def _emit_input_summary(self, data: dict[str, Any]) -> None:
        """Emit a compact structured summary of the generated input data.

        Parameters
        ----------
        data : dict
            Normalized YAML mapping before serialization.
        """
        table = phonon_input_summary_table(
            data,
            interface=self.interface_flag,
            source_count=len(self.files),
            eigenvectors_available=self.last_tracking_result is not None,
        )
        self.emit(
            "Phonon input data assembled",
            level=EventLevel.RESULT,
            data={
                "kind": "phonon_input_summary",
                "table": table,
            },
        )

    def to_yaml_lines(
        self,
        jobname: str,
        reference: int = 0,
        formula_units: int = 1,
    ) -> list[str]:
        """
        Build Quantas HA/QHA YAML lines from loaded interface data.

        This method is retained for compatibility with previous refactoring
        steps. Internally, data are first represented as a dictionary and then
        serialized with PyYAML.

        Parameters
        ----------
        jobname : str
            Short description stored in the ``job`` field.
        reference : int, optional
            Reference dataset for multiple-file input generation.
        formula_units : int, optional
            Number of chemical formula units in the thermodynamic normalization cell.

        Returns
        -------
        list of str
            Lines of the generated YAML file.
        """
        text = format_quantas_yaml(
            self.to_dict(
                jobname=jobname,
                reference=reference,
                formula_units=formula_units,
            )
        )
        return text.rstrip().splitlines()

    def write(
        self,
        outfile: str | Path,
        jobname: str,
        ref: int = 0,
        formula_units: int = 1,
    ) -> Path:
        """
        Write a Quantas HA/QHA YAML input file.

        Parameters
        ----------
        outfile : str or pathlib.Path
            Destination YAML file.
        jobname : str
            Short description stored in the ``job`` field.
        ref : int, optional
            Reference index for multiple-file input generation.
        formula_units : int, optional
            Number of chemical formula units in the thermodynamic normalization cell.

        Returns
        -------
        pathlib.Path
            Path of the written file.

        Raises
        ------
        ValueError
            If no data have been loaded or the reference index is invalid.
        OSError
            If the output file cannot be written.
        """
        outfile = Path(outfile)
        text = format_quantas_yaml(
            self.to_dict(
                jobname=jobname,
                reference=ref,
                formula_units=formula_units,
            )
        )
        outfile.write_text(text, encoding="utf-8")
        return outfile

    def _dict_from_multiple_readers(
        self,
        phondata: list[Any],
        jobname: str,
        reference: int,
        formula_units: int,
    ) -> dict[str, Any]:
        """
        Generate YAML data from several single-volume phonon readers.

        Parameters
        ----------
        phondata : list
            Loaded interface reader objects.
        jobname : str
            Input description.
        reference : int
            Reference reader index.
        formula_units : int
            Number of formula units in the thermodynamic normalization cell.

        Returns
        -------
        dict
            Quantas HA/QHA input data.
        """
        nvol = len(phondata)
        ref = phondata[reference]
        _validate_multiple_reader_qmeshes(phondata, reference)
        _validate_multiple_reader_units(phondata, reference)
        data = _header_dict(jobname, ref, formula_units)
        data.update(_q_position_metadata(ref))
        data["volume"] = [_as_float(getattr(item, "volume")) for item in phondata]
        data["energy"] = [_as_float(getattr(item, "energy")) for item in phondata]
        structure = _combine_structure_series(phondata, reference)
        if structure is not None:
            data["structure"] = structure.as_dict(include_source=False)
        tracking = _track_multiple_reader_modes(phondata, reference) if nvol > 1 else None
        self.last_tracking_result = tracking
        if nvol > 1 and tracking is None:
            data["mode_continuity"] = "unknown"
            data["mode_continuity_metadata"] = {
                "method": "none",
                "reason": "phonon_eigenvectors_unavailable",
            }
            tracked_frequencies = None
            self.emit(
                "Phonon-mode continuity could not be verified because one or "
                "more sources do not expose eigenvectors",
                level=EventLevel.WARNING,
                data={"kind": "mode_tracking_unavailable"},
            )
        elif tracking is not None:
            data["mode_continuity"] = tracking.status
            data["mode_continuity_metadata"] = _mode_tracking_metadata(tracking)
            tracked_frequencies = tracking.frequencies
            volumes = np.asarray(data["volume"], dtype=np.float64)
            qcoords = _fractional_qcoords(ref, int(getattr(ref, "qpoints")))
            frequency_unit = _reader_units(ref)["frequency"]
            self.emit(
                "Phonon-mode continuity analysis completed",
                level=EventLevel.RESULT,
                data={
                    "kind": "mode_tracking_summary",
                    "table": mode_tracking_summary_table(tracking, volumes),
                },
            )
            for step in tracking.steps:
                self.emit(
                    "Detailed phonon-mode continuity diagnostics",
                    level=EventLevel.DEBUG,
                    data={
                        "kind": "mode_tracking_detail",
                        "table": mode_tracking_detail_table(
                            tracking,
                            step,
                            volumes,
                            qcoords=qcoords,
                            frequency_unit=frequency_unit,
                        ),
                    },
                )
            if tracking.unresolved_assignments:
                self.emit(
                    f"{tracking.unresolved_assignments} phonon-mode assignment(s) "
                    "remain unresolved; frequency-based QHA should not be used",
                    level=EventLevel.WARNING,
                    data={
                        "kind": "mode_tracking_unresolved",
                        "count": tracking.unresolved_assignments,
                    },
                )
            elif tracking.caution_assignments:
                self.emit(
                    f"Mode continuity verified with {tracking.caution_assignments} "
                    "assignment caution(s)",
                    level=EventLevel.WARNING,
                    data={
                        "kind": "mode_tracking_cautions",
                        "count": tracking.caution_assignments,
                    },
                )
        else:
            tracked_frequencies = None
        data["phonon"] = _phonon_entries(
            ref,
            nvol,
            phondata=phondata,
            tracked_frequencies=tracked_frequencies,
        )
        return data

    @staticmethod
    def _dict_from_single_reader(
        reader: Any,
        jobname: str,
        formula_units: int,
    ) -> dict[str, Any]:
        """
        Generate YAML data from a single CRYSTAL-QHA reader.

        Parameters
        ----------
        reader : object
            Loaded interface reader with QHA arrays.
        jobname : str
            Input description.
        formula_units : int
            Number of formula units in the thermodynamic normalization cell.

        Returns
        -------
        dict
            Quantas HA/QHA input data.
        """
        nvol = int(getattr(reader, "points"))
        data = _header_dict(jobname, reader, formula_units)
        data.update(_q_position_metadata(reader))
        data["volume"] = _as_float_list(getattr(reader, "volume"))
        data["energy"] = _as_float_list(getattr(reader, "energy"))
        structure = getattr(reader, "structure_series", None)
        if structure is not None:
            data["structure"] = structure.as_dict(include_source=False)
        data["mode_continuity"] = str(
            getattr(reader, "mode_continuity", "unknown")
        )
        continuity_metadata = getattr(reader, "mode_continuity_metadata", {})
        data["mode_continuity_metadata"] = dict(continuity_metadata)
        data["phonon"] = _phonon_entries(reader, nvol, single_reader=reader)
        return data


class _FlowSequence(list):
    """YAML sequence rendered in flow style without changing its values."""


class _QuantasYamlDumper(yaml.SafeDumper):
    """Safe YAML dumper with presentation-only sequence controls."""


def _represent_flow_sequence(
    dumper: _QuantasYamlDumper,
    data: _FlowSequence,
) -> yaml.nodes.SequenceNode:
    """Represent one sequence as ``[a, b, c]`` in generated YAML."""
    node = dumper.represent_list(list(data))
    node.flow_style = True
    return node


_QuantasYamlDumper.add_representer(_FlowSequence, _represent_flow_sequence)


def _yaml_presentation_data(value: Any, *, path: tuple[str, ...] = ()) -> Any:
    """Return YAML-safe data with presentation-only sequence annotations.

    The returned object preserves all scalar values and mapping keys.  Selected
    vectors are rendered on one line and selected matrices keep one flow-style
    vector per physical row.  This function deliberately changes no numerical
    precision and is used only by :func:`format_quantas_yaml`.

    Parameters
    ----------
    value : object
        Plain or NumPy-rich data to prepare for YAML serialization.
    path : tuple of str, optional
        Mapping path used to identify vectors and matrices whose layout should
        be compacted.

    Returns
    -------
    object
        YAML-safe dictionaries/lists with selected lists tagged for flow style.
    """
    plain = _plain_data(value)
    if isinstance(plain, dict):
        return {
            str(key): _yaml_presentation_data(item, path=(*path, str(key)))
            for key, item in plain.items()
        }
    if not isinstance(plain, list):
        return plain

    key = path[-1] if path else ""
    inline_vector = key in {
        "atomic_numbers",
        "volume_order",
        "compression",
        "expansion",
        "equivalent_atoms",
        "origin_shift",
    } or path[-2:] == ("volume_series", "volume")
    row_vector_container = key in {
        "lattice",
        "fractional_positions",
        "expansion_matrix",
        "transformation_matrix",
        "primitive_to_crystallographic",
    }

    if inline_vector and _is_scalar_sequence(plain):
        return _FlowSequence(plain)
    if row_vector_container:
        return _compact_yaml_rows(plain)
    return [
        _yaml_presentation_data(item, path=path)
        if isinstance(item, (dict, list))
        else item
        for item in plain
    ]


def _compact_yaml_rows(value: list[Any]) -> list[Any] | _FlowSequence:
    """Render the innermost numeric vectors of one array-like value inline."""
    if _is_scalar_sequence(value):
        return _FlowSequence(value)
    return [
        _compact_yaml_rows(item) if isinstance(item, list) else item
        for item in value
    ]


def _is_scalar_sequence(value: list[Any]) -> bool:
    """Return whether one list contains no nested mappings or sequences."""
    return all(not isinstance(item, (dict, list, tuple)) for item in value)


def _dump_yaml_section(key: str, value: Any) -> str:
    """Serialize one mapping section with Quantas presentation conventions."""
    prepared = _yaml_presentation_data({key: value})
    return yaml.dump(
        prepared,
        Dumper=_QuantasYamlDumper,
        sort_keys=False,
        default_flow_style=False,
        width=4096,
    ).rstrip()

def format_quantas_yaml(data: dict[str, Any]) -> str:
    """
    Serialize Quantas HA/QHA input data using the reference readable layout.

    Parameters
    ----------
    data : dict
        Quantas HA/QHA input mapping generated by :class:`HAInputCreator`.

    Returns
    -------
    str
        YAML text preserving the compact matrix, volume, energy, q-point, and
        frequency-array layout used by Quantas input files.
    """
    lines: list[str] = []
    lines.append(f"job: {data.get('job', 'Quantas HA input')}")
    lines.append(f"natom: {int(data['natom']):3d}")
    lines.append(f"formula_units: {int(data.get('formula_units', 1))}")
    units_text = _dump_yaml_section(
        "units", data.get("units", default_phonon_input_units())
    )
    lines.extend(units_text.splitlines())
    if "provenance" in data:
        provenance_text = _dump_yaml_section("provenance", data["provenance"])
        lines.extend(provenance_text.splitlines())
    if np.asarray(data.get("volume", [])).size > 1:
        continuity = str(data.get("mode_continuity", "unknown"))
        lines.append(f"mode_continuity: {continuity}")
        continuity_metadata = data.get("mode_continuity_metadata")
        if continuity_metadata:
            metadata_text = _dump_yaml_section(
                "mode_continuity_metadata", continuity_metadata
            )
            lines.extend(metadata_text.splitlines())
    lines.append("supercell:")
    for row in data["supercell"]:
        values = ", ".join(f"{int(value):7d}" for value in row)
        lines.append(f"- [ {values} ]")
    if "structure" in data:
        structure_text = _dump_yaml_section("structure", data["structure"])
        lines.extend(structure_text.splitlines())
    if "q_position_source" in data:
        lines.append(f"q_position_source: {data['q_position_source']}")
    if "q_position_convention" in data:
        lines.append(f"q_position_convention: {data['q_position_convention']}")
    if "q_position_note" in data:
        note_text = _dump_yaml_section(
            "q_position_note", str(data["q_position_note"])
        )
        lines.extend(note_text.splitlines())
    lines.append(f"qpoints: {int(data['qpoints'])}")
    lines.append(f"volume: {_format_float_sequence(data['volume'], precision=8)}")
    lines.append(f"energy: {_format_energy_sequence(data['energy'])}")
    lines.append("phonon:")

    for qpoint in data["phonon"]:
        q_position = qpoint.get("q-position")
        if q_position is None:
            lines.append("- q-position: null")
        else:
            qpos = _format_float_sequence(q_position, precision=7)
            lines.append(f"- q-position: {qpos}")
        lines.append(f"  weight: {_format_weight(qpoint['weight'])}")
        lines.append("  band:")
        for index, band in enumerate(qpoint["band"], start=1):
            lines.append(f"  - # {index}")
            lines.append(
                "    frequency: "
                + _format_float_sequence(band["frequency"], precision=10)
            )
    return "\n".join(lines) + "\n"


def _format_float_sequence(values: Any, precision: int = 8) -> str:
    """
    Format a numeric sequence as a compact YAML flow-style list.

    Parameters
    ----------
    values : object
        Numeric scalar or sequence.
    precision : int, optional
        Number of decimal places.

    Returns
    -------
    str
        Formatted flow-style list.
    """
    array = np.atleast_1d(np.asarray(values, dtype=np.float64))
    numbers = ", ".join(f"{float(value):16.{precision}f}" for value in array)
    return f"[ {numbers} ]"


def _format_energy_sequence(values: Any) -> str:
    """
    Format static energies as a compact YAML flow-style list.

    Parameters
    ----------
    values : object
        Numeric scalar or sequence.

    Returns
    -------
    str
        Formatted flow-style list using scientific notation.
    """
    array = np.atleast_1d(np.asarray(values, dtype=np.float64))
    numbers = ", ".join(f"{float(value):19.12E}" for value in array)
    return f"[ {numbers} ]"


def _format_weight(value: Any) -> str:
    """
    Format a q-point weight as a YAML scalar.

    Parameters
    ----------
    value : object
        Weight value.

    Returns
    -------
    str
        Integer-like weights are written without decimal places.
    """
    number = float(value)
    if number.is_integer():
        return f"{int(number)}"
    return f"{number:.10g}"


def create_ha_input(
    filename: str | Path,
    outfile: str | Path,
    *,
    interface: str = "crystal",
    is_list: bool = False,
    reference: int = 0,
    jobname: str = "Quantas HA input",
    formula_units: int = 1,
    interface_filter: dict[str, type[Any] | None] | None = None,
    observer: Observer | None = None,
) -> Path:
    """
    Create a Quantas HA/QHA YAML input file from QM output data.

    Parameters
    ----------
    filename : str or pathlib.Path
        QM output file or file list.
    outfile : str or pathlib.Path
        Destination YAML file.
    interface : str, optional
        Input interface name: ``"crystal"``, ``"crystal-qha"``, or
        ``"phonopy"``.
    is_list : bool, optional
        If ``True``, ``filename`` is interpreted as a text file listing QM
        output files.
    reference : int, optional
        Reference input index for multiple-file generation.
    jobname : str, optional
        Description written in the YAML ``job`` field.
    formula_units : int, optional
        Number of chemical formula units in the thermodynamic normalization cell.
    interface_filter : dict or None, optional
        Optional reader mapping overriding the default interface readers.
    observer : Observer or None, optional
        Frontend-neutral observer receiving input-generation events.

    Returns
    -------
    pathlib.Path
        Path of the generated YAML file.

    Raises
    ------
    ValueError
        If the selected interface cannot read the provided input data.
    OSError
        If files cannot be read or written.
    """
    creator = HAInputCreator(
        interface=interface,
        interface_filter=interface_filter,
        observer=observer,
    )
    completed, error = creator.read(filename, is_list=is_list, reference=reference)
    if not completed:
        raise ValueError(error or "Unable to create HA input")
    return creator.write(
        outfile,
        jobname=jobname,
        ref=reference,
        formula_units=formula_units,
    )


QHAInputCreator = HAInputCreator


def _combine_structure_series(
    readers: list[Any],
    reference_index: int,
) -> StructureVolumeSeries | None:
    """Combine one-volume reader structures into a volume-dependent path.

    Parameters
    ----------
    readers : list
        Loaded single-volume interface readers.
    reference_index : int
        Reader fixing the reference atom order, symmetry, and orientation.

    Returns
    -------
    StructureVolumeSeries or None
        Combined structural path, or ``None`` when readers do not expose
        structural information.

    Raises
    ------
    ValueError
        If only part of the reader series has structural data or compact
        structures are inconsistent.
    """
    series = [getattr(reader, "structure_series", None) for reader in readers]
    if all(item is None for item in series):
        return None
    if any(item is None for item in series):
        raise ValueError(
            "structural data are available for only part of the phonon series"
        )
    values = [item for item in series if item is not None]
    reference_series = values[reference_index]
    reference = reference_series.reference
    for index, item in enumerate(values):
        if item.nvol != 1:
            raise ValueError(
                f"single-volume reader {index} exposes {item.nvol} structures"
            )
        if not np.array_equal(item.reference.atomic_numbers, reference.atomic_numbers):
            raise ValueError(
                "primitive atomic species/order changes across phonon inputs"
            )
        if item.normalization.basis != reference_series.normalization.basis:
            raise ValueError(
                "phonon inputs use inconsistent thermodynamic cell normalization"
            )
    diagnostics = tuple(item.diagnostics[0] for item in values if item.diagnostics)
    if diagnostics and len(diagnostics) != len(values):
        raise ValueError("structural diagnostics are incomplete across input files")
    return StructureVolumeSeries(
        reference=reference,
        lattices=np.asarray([item.lattices[0] for item in values], dtype=np.float64),
        fractional_positions=np.asarray(
            [item.fractional_positions[0] for item in values],
            dtype=np.float64,
        ),
        volumes=np.asarray([item.volumes[0] for item in values], dtype=np.float64),
        normalization=reference_series.normalization,
        symmetry=reference_series.symmetry,
        primitive_to_crystallographic=(reference_series.primitive_to_crystallographic),
        diagnostics=diagnostics,
        orientation=reference_series.orientation,
        reference_index=reference_index,
        metadata={
            "interface": "combined-volume-series",
            "source_count": len(values),
        },
    )


def _plain_data(value: Any) -> Any:
    """Convert NumPy-rich structures to YAML-safe Python containers.

    Parameters
    ----------
    value : object
        Recursively serializable value containing NumPy arrays or scalars.

    Returns
    -------
    object
        Plain dictionaries, lists, strings, integers, floats, and booleans.
    """
    if isinstance(value, dict):
        return {str(key): _plain_data(item) for key, item in value.items()}
    if isinstance(value, np.ndarray):
        return _plain_data(value.tolist())
    if isinstance(value, (list, tuple)):
        return [_plain_data(item) for item in value]
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    return value


def _header_dict(
    jobname: str,
    reader: Any,
    formula_units: int,
) -> dict[str, Any]:
    """
    Generate the common YAML header mapping.

    Parameters
    ----------
    jobname : str
        Input description.
    reader : object
        Reference interface reader.
    formula_units : int
        Number of formula units in the thermodynamic normalization cell.

    Returns
    -------
    dict
        Header data.
    """
    return {
        "job": jobname,
        "natom": int(getattr(reader, "natom")),
        "formula_units": int(formula_units),
        "units": _reader_units(reader),
        "supercell": _as_int_matrix(getattr(reader, "dim")),
        "qpoints": int(getattr(reader, "qpoints")),
    }


def _phonon_entries(
    reference: Any,
    nvol: int,
    *,
    phondata: list[Any] | None = None,
    single_reader: Any | None = None,
    tracked_frequencies: np.ndarray | None = None,
) -> list[dict[str, Any]]:
    """
    Generate the ``phonon`` section of the YAML input.

    Parameters
    ----------
    reference : object
        Reference reader providing q-point metadata.
    nvol : int
        Number of volumes.
    phondata : list or None, optional
        Multiple single-volume readers.
    single_reader : object or None, optional
        Single CRYSTAL-QHA reader.
    tracked_frequencies : ndarray or None, optional
        Mode-continuous frequencies with shape ``(nvol, qpoints, nphonon)``.

    Returns
    -------
    list of dict
        Q-point, weight and band entries for the YAML input.

    Raises
    ------
    ValueError
        If frequency arrays are inconsistent with the number of volumes.
    """
    qpoints = int(getattr(reference, "qpoints"))
    nphonon = int(getattr(reference, "nphonon"))
    qcoords = _fractional_qcoords(reference, qpoints)
    weights = _weights_array(getattr(reference, "weights"), qpoints)
    phonons = (
        getattr(single_reader, "phonons", None) if single_reader is not None else None
    )

    entries: list[dict[str, Any]] = []
    for i in range(qpoints):
        q_position = None if qcoords is None else qcoords[i].astype(float).tolist()
        bands: list[dict[str, list[float]]] = []
        for j in range(nphonon):
            if tracked_frequencies is not None:
                frequencies = np.asarray(
                    tracked_frequencies[:, i, j],
                    dtype=np.float64,
                )
            elif single_reader is not None:
                frequencies = _phonon_frequency_array(phonons, i, j)
            else:
                if phondata is None:
                    raise ValueError("missing phonon data")
                frequencies = np.asarray(
                    [_phonon_scalar_or_array(item.phonons, i, j) for item in phondata],
                    dtype=np.float64,
                )
            frequencies = np.atleast_1d(np.asarray(frequencies, dtype=np.float64))
            if frequencies.shape[0] != nvol:
                raise ValueError("phonon frequency arrays do not match volumes")
            bands.append({"frequency": frequencies.astype(float).tolist()})
        entries.append(
            {
                "q-position": q_position,
                "weight": _weight_value(weights[i]),
                "band": bands,
            }
        )
    return entries


def _track_multiple_reader_modes(
    readers: list[Any],
    reference_index: int,
) -> PhononModeTrackingResult | None:
    """Return mode-continuous frequencies when all readers expose eigenvectors.

    Parameters
    ----------
    readers : list
        Loaded single-volume phonon readers.
    reference_index : int
        Reader fixing phonon branch labels.

    Returns
    -------
    PhononModeTrackingResult or None
        Tracking result, or ``None`` when at least one interface reader does
        not expose phonon eigenvectors.

    Raises
    ------
    ValueError
        If eigenvector dimensions, atoms, q points, or frequency arrays are
        inconsistent across sampled volumes.
    """
    mode_data = [getattr(reader, "mode_data", None) for reader in readers]
    if any(item is None for item in mode_data):
        return None

    available = [item for item in mode_data if item is not None]
    reference = available[reference_index]
    for item in available:
        if item.atom_symbols != reference.atom_symbols:
            raise ValueError(
                "phonon eigenvector atom ordering changes between volume points"
            )
        if item.eigenvectors.shape != reference.eigenvectors.shape:
            raise ValueError(
                "phonon eigenvector dimensions change between volume points"
            )

    frequencies = np.asarray(
        [reader.phonons_array() for reader in readers],
        dtype=np.float64,
    )
    eigenvectors = np.asarray(
        [item.eigenvectors for item in available],
        dtype=np.complex128,
    )
    if frequencies.shape != eigenvectors.shape[:3]:
        raise ValueError(
            "phonon frequencies and eigenvectors use incompatible dimensions"
        )
    volumes = np.asarray(
        [_as_float(getattr(reader, "volume")) for reader in readers],
        dtype=np.float64,
    )
    return track_phonon_modes(
        frequencies,
        eigenvectors,
        volumes,
        reference_index=reference_index,
    )


def _mode_tracking_metadata(
    result: PhononModeTrackingResult,
) -> dict[str, Any]:
    """Return compact YAML metadata for one phonon-mode tracking result."""
    fit = result.fit
    fit_r2 = fit.r_squared[np.isfinite(fit.r_squared)]
    fit_rmse = fit.rmse[np.isfinite(fit.rmse)]
    fit_residual = fit.max_residual[np.isfinite(fit.max_residual)]
    volume_order = [int(index) for index in result.volume_order]
    reference_position = volume_order.index(int(result.reference_index))
    return {
        "method": "eigenvector_overlap",
        "reference_index": int(result.reference_index),
        "volume_order": volume_order,
        "traversal": {
            "compression": [
                int(result.reference_index),
                *reversed(volume_order[:reference_position]),
            ],
            "expansion": [
                int(result.reference_index),
                *volume_order[reference_position + 1 :],
            ],
        },
        "minimum_overlap": _finite_or_none(result.minimum_overlap),
        "minimum_subspace_singular_value": _finite_or_none(
            result.minimum_subspace_singular_value
        ),
        "reordered_assignments": int(result.reordered_assignments),
        "local_reordered_assignments": int(result.local_reordered_assignments),
        "ambiguous_assignments": int(result.ambiguous_assignments),
        "low_overlap_assignments": int(result.low_overlap_assignments),
        "caution_assignments": int(result.caution_assignments),
        "unresolved_assignments": int(result.unresolved_assignments),
        "degenerate_subspaces": int(result.degenerate_subspaces),
        "unresolved_degenerate_subspaces": int(
            result.unresolved_degenerate_subspaces
        ),
        "ambiguity_margin": 0.4,
        "minimum_required_overlap": 0.5,
        "degeneracy_atol_cm-1": 5.0e-2,
        "degeneracy_rtol": 1.0e-6,
        "minimum_subspace_overlap": 0.8,
        "frequency_fit": {
            "model": "polynomial",
            "maximum_degree": 3,
            "degree": fit.degree,
            "residual_degrees_of_freedom": int(
                fit.residual_degrees_of_freedom
            ),
            "minimum_r_squared": (
                float(np.min(fit_r2)) if fit_r2.size else None
            ),
            "maximum_rmse_cm-1": (
                float(np.max(fit_rmse)) if fit_rmse.size else None
            ),
            "maximum_residual_cm-1": (
                float(np.max(fit_residual)) if fit_residual.size else None
            ),
            "supported_branches": int(np.count_nonzero(fit.supported)),
            "total_branches": int(fit.supported.size),
            "support_min_r_squared": 0.98,
            "support_max_rmse_cm-1": 2.0,
            "leave_one_out_validation": {
                "degree": fit.predictive_degree,
                "residual_degrees_of_freedom": int(
                    fit.predictive_residual_degrees_of_freedom
                ),
                "absolute_residual_floor_cm-1": 2.0,
                "maximum_training_range_fraction": 0.10,
                "supported_low_overlap_assignments": int(
                    sum(
                        np.count_nonzero(
                            step.low_overlap_mask & ~step.unresolved_mask
                        )
                        for step in result.steps
                    )
                ),
                "unresolved_low_overlap_assignments": int(
                    sum(
                        np.count_nonzero(
                            step.low_overlap_mask & step.unresolved_mask
                        )
                        for step in result.steps
                    )
                ),
            },
        },
    }


def phonon_input_summary_table(
    data: dict[str, Any],
    *,
    interface: str,
    source_count: int,
    eigenvectors_available: bool,
) -> ReportTable:
    """Build a compact neutral table describing generated phonon input data.

    Parameters
    ----------
    data : dict
        Generated HA/QHA input mapping.
    interface : str
        Interface used to read source calculations.
    source_count : int
        Number of source files read.
    eigenvectors_available : bool
        Whether Quantas had eigenvectors available for its own mode tracking.

    Returns
    -------
    ReportTable
        Frontend-neutral input-generation summary.
    """
    volumes = np.atleast_1d(np.asarray(data.get("volume", []), dtype=np.float64))
    phonon = data.get("phonon", [])
    modes = 0
    if isinstance(phonon, list) and phonon and isinstance(phonon[0], dict):
        bands = phonon[0].get("band", [])
        if isinstance(bands, list):
            modes = len(bands)
    if interface == "crystal-qha":
        eigenvectors = "source-managed"
    else:
        eigenvectors = "available" if eigenvectors_available else "unavailable"

    rows: list[list[Any]] = [
        ["Interface", interface],
        ["Source files", int(source_count)],
        ["Atoms", int(data.get("natom", 0))],
        ["Volumes", int(volumes.size)],
        ["Q-points", int(data.get("qpoints", 0))],
        ["Modes per q-point", int(modes)],
        ["Eigenvectors", eigenvectors],
    ]
    if volumes.size:
        rows.append(
            [
                "Volume range (angstrom^3)",
                f"{float(np.min(volumes)):.6f} -> {float(np.max(volumes)):.6f}",
            ]
        )
    if volumes.size > 1:
        rows.append(["Mode continuity", str(data.get("mode_continuity", "unknown"))])
    return ReportTable(
        title="Phonon input generation",
        columns=["Property", "Value"],
        rows=rows,
    )


def mode_tracking_summary_table(
    result: PhononModeTrackingResult,
    volumes: ArrayLike,
) -> ReportTable:
    """Build the standard adjacent-volume mode-tracking summary table.

    Parameters
    ----------
    result : PhononModeTrackingResult
        Completed backend-neutral tracking result.
    volumes : array_like
        Original sampled volumes in cubic angstrom.

    Returns
    -------
    ReportTable
        Frontend-neutral summary with one row per adjacent volume interval.
    """
    volume = np.asarray(volumes, dtype=np.float64)
    grouped: dict[tuple[int, int], list[PhononModeTrackingStep]] = {}
    for step in result.steps:
        grouped.setdefault((step.predecessor_index, step.source_index), []).append(step)

    rows: list[list[Any]] = []
    for (lower, upper), steps in grouped.items():
        reordered = sum(
            np.count_nonzero(step.permutation != np.arange(step.permutation.size))
            for step in steps
        )
        degenerate = sum(len(step.degenerate_subspaces) for step in steps)
        low = sum(np.count_nonzero(step.low_overlap_mask) for step in steps)
        cautions = sum(np.count_nonzero(step.caution_mask) for step in steps)
        unresolved = sum(np.count_nonzero(step.unresolved_mask) for step in steps)
        nondegenerate = [
            step.overlaps[~step.degenerate_mask]
            for step in steps
            if np.any(~step.degenerate_mask)
        ]
        min_overlap = (
            float(np.min(np.concatenate(nondegenerate)))
            if nondegenerate
            else float("nan")
        )
        singular = [
            value
            for step in steps
            for value in step.subspace_min_singular_values
        ]
        min_singular = float(min(singular)) if singular else float("nan")
        rows.append(
            [
                float(volume[lower]),
                float(volume[upper]),
                int(reordered),
                int(degenerate),
                int(low),
                int(cautions),
                int(unresolved),
                min_overlap,
                min_singular,
            ]
        )

    fit = result.fit
    fit_text = "not available"
    if fit.degree is not None:
        fit_text = (
            f"degree {fit.degree}, residual dof {fit.residual_degrees_of_freedom}"
        )
    status = result.status
    if result.verified and result.caution_assignments:
        status = "verified with cautions"
    notes = [
        f"Overall mode continuity: {status}.",
        (
            f"Ambiguous assignments: {result.ambiguous_assignments}; "
            f"low-overlap assignments: {result.low_overlap_assignments}."
        ),
        (
            f"Final reordered assignments: {result.reordered_assignments}; "
            f"local reordered assignments: {result.local_reordered_assignments}."
        ),
        f"Global frequency-path diagnostic: polynomial fit ({fit_text}).",
        (
            "Low-overlap validation: symmetric leave-one-out prediction "
            f"(degree {fit.predictive_degree}, residual dof "
            f"{fit.predictive_residual_degrees_of_freedom})."
            if fit.predictive_degree is not None
            else "Low-overlap validation: leave-one-out prediction unavailable."
        ),
    ]
    return ReportTable(
        title="Phonon mode continuity",
        columns=[
            "Volume from",
            "Volume to",
            "Local reorder",
            "Deg. subspaces",
            "Low overlap",
            "Cautions",
            "Unresolved",
            "Min overlap",
            "Min sigma",
        ],
        rows=rows,
        metadata={
            "column_units": [
                "angstrom^3",
                "angstrom^3",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
            ],
            "column_formats": [
                ".6f",
                ".6f",
                "integer",
                "integer",
                "integer",
                "integer",
                "integer",
                ".4f",
                ".4f",
            ],
            "column_alignments": ["right"] * 9,
            "notes": notes,
        },
    )


def mode_tracking_detail_table(
    result: PhononModeTrackingResult,
    step: PhononModeTrackingStep,
    volumes: ArrayLike,
    *,
    qcoords: ArrayLike | None = None,
    frequency_unit: str = "cm^-1",
) -> ReportTable:
    """Build one detailed q-point table for debug mode-tracking output.

    Parameters
    ----------
    result : PhononModeTrackingResult
        Complete tracking result.
    step : PhononModeTrackingStep
        Adjacent-volume q-point assignment to describe.
    volumes : array_like
        Original sampled volumes in cubic angstrom.
    qcoords : array_like or None, optional
        Fractional q-point coordinates.
    frequency_unit : str, optional
        Frequency unit displayed by the renderer.

    Returns
    -------
    ReportTable
        Mode-by-mode diagnostics ordered by final branch label.
    """
    volume = np.asarray(volumes, dtype=np.float64)
    coordinate_text = "unavailable"
    if qcoords is not None:
        coordinates = np.asarray(qcoords, dtype=np.float64)
        q = coordinates[step.qpoint_index]
        coordinate_text = f"({q[0]:.6f}, {q[1]:.6f}, {q[2]:.6f})"
    title = (
        "Mode tracking: "
        f"{volume[step.predecessor_index]:.6f} -> "
        f"{volume[step.source_index]:.6f} (angstrom^3); "
        f"q #{step.qpoint_index + 1} = {coordinate_text}"
    )

    sigma_by_raw: dict[int, float] = {}
    for groups, singular_values in (
        (
            step.degenerate_subspaces,
            step.degenerate_subspace_min_singular_values,
        ),
        (
            step.unresolved_degenerate_subspaces,
            step.unresolved_degenerate_subspace_min_singular_values,
        ),
    ):
        for group, singular in zip(groups, singular_values, strict=True):
            for raw_mode in group:
                sigma_by_raw[int(raw_mode)] = float(singular)

    rows: list[list[Any]] = []
    order = np.argsort(step.branch_indices, kind="stable")
    for raw_from in order:
        branch = int(step.branch_indices[raw_from])
        raw_to = int(step.permutation[raw_from])
        if step.unresolved_mask[raw_from]:
            status = "unresolved"
        elif step.degenerate_mask[raw_from]:
            status = "degenerate"
        elif step.caution_mask[raw_from]:
            status = "caution"
        else:
            status = "matched"
        competitor: float | None = None
        gap: float | None = None
        if not step.degenerate_mask[raw_from]:
            competitor = float(step.competitor_overlaps[raw_from])
            gap = float(step.overlap_gaps[raw_from])
        fit_r2 = float(result.fit.r_squared[step.qpoint_index, branch])
        fit_rmse = float(result.fit.rmse[step.qpoint_index, branch])
        loo_residual: float | None = None
        loo_limit: float | None = None
        if step.low_overlap_mask[raw_from]:
            endpoint_diagnostics = []
            for state_index in (step.predecessor_index, step.source_index):
                residual = float(
                    result.fit.predictive_residuals[
                        state_index, step.qpoint_index, branch
                    ]
                )
                tolerance = float(
                    result.fit.predictive_tolerances[
                        state_index, step.qpoint_index, branch
                    ]
                )
                if np.isfinite(residual) and np.isfinite(tolerance):
                    ratio = residual / tolerance if tolerance > 0.0 else float("inf")
                    endpoint_diagnostics.append((ratio, residual, tolerance))
            if endpoint_diagnostics:
                _, loo_residual, loo_limit = max(endpoint_diagnostics)
        rows.append(
            [
                branch + 1,
                int(raw_from) + 1,
                raw_to + 1,
                float(result.frequencies[step.predecessor_index, step.qpoint_index, branch]),
                float(result.frequencies[step.source_index, step.qpoint_index, branch]),
                None if step.degenerate_mask[raw_from] else float(step.overlaps[raw_from]),
                competitor,
                gap,
                sigma_by_raw.get(int(raw_from)),
                status,
                loo_residual,
                loo_limit,
                fit_r2 if np.isfinite(fit_r2) else None,
                fit_rmse if np.isfinite(fit_rmse) else None,
            ]
        )

    return ReportTable(
        title=title,
        columns=[
            "Branch",
            "Raw mode from",
            "Raw mode to",
            "Frequency from",
            "Frequency to",
            "Overlap",
            "Competitor",
            "Gap",
            "Subspace sigma min",
            "Status",
            "LOO residual",
            "LOO limit",
            "Global fit R^2",
            "Global fit RMSE",
        ],
        rows=rows,
        metadata={
            "column_units": [
                "",
                "",
                "",
                frequency_unit,
                frequency_unit,
                "",
                "",
                "",
                "",
                "",
                frequency_unit,
                frequency_unit,
                "",
                frequency_unit,
            ],
            "column_formats": [
                "integer",
                "integer",
                "integer",
                ".4f",
                ".4f",
                ".4f",
                ".4f",
                ".4f",
                ".4f",
                None,
                ".4f",
                ".4f",
                ".6f",
                ".4f",
            ],
            "column_alignments": [
                "right",
                "right",
                "right",
                "right",
                "right",
                "right",
                "right",
                "right",
                "right",
                "left",
                "right",
                "right",
                "right",
                "right",
            ],
        },
    )


def _finite_or_none(value: float) -> float | None:
    """Return one finite diagnostic value or ``None`` for YAML output."""
    return float(value) if np.isfinite(value) else None


def _reader_units(reader: Any) -> dict[str, str]:
    """Return validated physical units exposed by one interface reader.

    Readers predating explicit unit metadata retain the historical Quantas
    HA/QHA convention of Hartree, cubic angstrom, and wavenumbers.

    Parameters
    ----------
    reader : object
        Loaded interface reader.

    Returns
    -------
    dict
        Energy, volume, frequency, and structural length unit labels.

    Raises
    ------
    ValueError
        If a reader exposes an incomplete or empty unit mapping.
    """
    raw = getattr(reader, "units", None)
    if raw is None:
        return default_phonon_input_units()
    if not isinstance(raw, dict):
        raise ValueError("interface reader units must be a mapping")

    units = {str(key): str(value).strip() for key, value in raw.items()}
    missing = [key for key in default_phonon_input_units() if not units.get(key)]
    if missing:
        missing_text = ", ".join(missing)
        raise ValueError(f"interface reader units missing: {missing_text}")
    return {key: units[key] for key in default_phonon_input_units()}


def _provenance_dict(
    interface: str,
    files: list[Path],
    *,
    reference_index: int,
) -> dict[str, Any]:
    """Return compact source provenance for one generated HA/QHA input.

    Parameters
    ----------
    interface : str
        Interface name selected by the input generator.
    files : list of pathlib.Path
        Source files read by the interface.
    reference_index : int
        Source index used as the reference dataset.

    Returns
    -------
    dict
        Interface name, source identifiers, and reference index.
    """
    return {
        "interface": str(interface),
        "sources": [str(path) for path in files],
        "reference_index": int(reference_index),
    }


def _q_position_metadata(reader: Any) -> dict[str, str]:
    """Return YAML metadata describing q-point coordinate provenance.

    Parameters
    ----------
    reader : object
        Interface reader providing optional ``q_position_source`` metadata.

    Returns
    -------
    dict of str
        Top-level YAML metadata for q-point coordinates.
    """
    source = str(getattr(reader, "q_position_source", "reader-provided"))
    metadata = {"q_position_source": source}
    if source.startswith("unavailable"):
        metadata["q_position_convention"] = "unavailable"
        metadata["q_position_note"] = (
            "Native CRYSTAL QHA follows supercell Gamma eigenmodes with volume "
            "but does not print a reliable mapping from each stored mode block "
            "to primitive-cell q-point coordinates. Equal weights preserve the "
            "thermodynamic normalization; q-position labels are therefore null."
        )
    else:
        metadata["q_position_convention"] = (
            "fractional primitive reciprocal coordinates"
        )
    return metadata


def _fractional_qcoords(reader: Any, qpoints: int) -> np.ndarray | None:
    """Return normalized fractional q-point coordinates from a reader.

    Parameters
    ----------
    reader : object
        Interface reader exposing q-point coordinates.
    qpoints : int
        Number of q-points expected in the result.

    Returns
    -------
    ndarray or None
        Fractional coordinates with shape ``(qpoints, 3)``, or ``None`` when
        coordinates are scientifically unavailable.

    Raises
    ------
    ValueError
        If available coordinates or shrinking factors have invalid shapes.
    """
    source = str(getattr(reader, "q_position_source", "reader-provided"))
    if source.startswith("unavailable"):
        return None

    fractional = getattr(reader, "qcoords_fractional", None)
    if callable(fractional):
        fractional = fractional()
    if fractional is not None:
        array = np.asarray(fractional, dtype=np.float64)
        if array.shape != (qpoints, 3):
            raise ValueError("fractional q-point coordinates have invalid shape")
        return array

    raw = _qcoords_array(getattr(reader, "qcoords"), qpoints)
    shrinkf = np.asarray(
        getattr(reader, "shrinkf", np.ones(3)),
        dtype=np.float64,
    )
    if shrinkf.shape != (3,) or np.any(shrinkf <= 0.0):
        raise ValueError("q-point shrinking factors must contain three positives")
    return raw / shrinkf[np.newaxis, :]


def _validate_multiple_reader_units(
    readers: list[Any],
    reference: int,
) -> None:
    """Require one consistent physical-unit contract across a file series.

    Parameters
    ----------
    readers : list
        Loaded single-volume interface readers.
    reference : int
        Reference reader index.

    Raises
    ------
    ValueError
        If one source reports units different from the reference source.
    """
    expected = _reader_units(readers[reference])
    for index, reader in enumerate(readers):
        if _reader_units(reader) != expected:
            raise ValueError(
                f"input file {index} uses physical units inconsistent "
                "with the reference source"
            )


def _validate_multiple_reader_qmeshes(
    readers: list[Any],
    reference_index: int,
) -> None:
    """Require identical q-point meshes across a multi-volume input series.

    Parameters
    ----------
    readers : list
        Loaded single-volume phonon readers.
    reference_index : int
        Reader used as q-point and mode-order reference.

    Raises
    ------
    ValueError
        If q-point counts, coordinates, weights, supercell matrices, or mode
        counts differ.  Frequencies from different volumes can only be stacked
        safely when all of these metadata use the same ordering.
    """
    reference = readers[reference_index]
    ref_qpoints = int(getattr(reference, "qpoints"))
    ref_modes = int(getattr(reference, "nphonon"))
    ref_coordinates = _fractional_qcoords(reference, ref_qpoints)
    ref_weights = _weights_array(getattr(reference, "weights"), ref_qpoints)
    ref_supercell = np.asarray(getattr(reference, "dim"), dtype=np.int64)

    for index, reader in enumerate(readers):
        qpoints = int(getattr(reader, "qpoints"))
        modes = int(getattr(reader, "nphonon"))
        if qpoints != ref_qpoints:
            raise ValueError(
                f"phonon input {index} has {qpoints} q-points; expected {ref_qpoints}"
            )
        if modes != ref_modes:
            raise ValueError(
                f"phonon input {index} has {modes} modes per q-point; "
                f"expected {ref_modes}"
            )
        supercell = np.asarray(getattr(reader, "dim"), dtype=np.int64)
        if not np.array_equal(supercell, ref_supercell):
            raise ValueError(f"phonon input {index} uses a different supercell matrix")
        coordinates = _fractional_qcoords(reader, qpoints)
        if (coordinates is None) != (ref_coordinates is None):
            raise ValueError(f"phonon input {index} has inconsistent q-point metadata")
        if coordinates is not None and ref_coordinates is not None:
            if not np.allclose(coordinates, ref_coordinates, rtol=0.0, atol=1.0e-12):
                raise ValueError(
                    f"phonon input {index} uses different q-point coordinates "
                    "or ordering"
                )
        weights = _weights_array(getattr(reader, "weights"), qpoints)
        if not np.allclose(weights, ref_weights, rtol=0.0, atol=1.0e-12):
            raise ValueError(f"phonon input {index} uses different q-point weights")


def _qcoords_array(values: Any, qpoints: int) -> np.ndarray:
    """
    Normalize q-point coordinates to a two-dimensional array.

    Parameters
    ----------
    values : object
        Coordinates as array-like object or dictionary indexed by q-point.
    qpoints : int
        Number of q-points.

    Returns
    -------
    numpy.ndarray
        Array with shape ``(qpoints, 3)``.
    """
    if isinstance(values, dict):
        return np.asarray([values[i] for i in range(qpoints)], dtype=np.float64)
    return np.asarray(values, dtype=np.float64).reshape(qpoints, 3)


def _weights_array(values: Any, qpoints: int) -> np.ndarray:
    """
    Normalize q-point weights to a one-dimensional array.

    Parameters
    ----------
    values : object
        Weights as array-like object or dictionary indexed by q-point.
    qpoints : int
        Number of q-points.

    Returns
    -------
    numpy.ndarray
        Array with shape ``(qpoints,)``.
    """
    if isinstance(values, dict):
        return np.asarray([values[i] for i in range(qpoints)], dtype=np.float64)
    return np.asarray(values, dtype=np.float64).reshape(qpoints)


def _phonon_frequency_array(phonons: Any, qpoint: int, band: int) -> np.ndarray:
    """
    Return all volume frequencies for one q-point and one phonon band.

    Parameters
    ----------
    phonons : object
        Phonon data as dictionary or array-like object.
    qpoint : int
        Q-point index.
    band : int
        Phonon band index.

    Returns
    -------
    numpy.ndarray
        Frequencies for all volumes.
    """
    if isinstance(phonons, dict):
        return np.asarray(phonons[qpoint][band], dtype=np.float64)
    return np.asarray(phonons[qpoint][band], dtype=np.float64)


def _phonon_scalar_or_array(phonons: Any, qpoint: int, band: int) -> Any:
    """
    Return the phonon entry for one single-volume reader.

    Parameters
    ----------
    phonons : object
        Phonon data as dictionary or array-like object.
    qpoint : int
        Q-point index.
    band : int
        Phonon band index.

    Returns
    -------
    object
        Scalar or array-like phonon frequency value.
    """
    if isinstance(phonons, dict):
        return phonons[qpoint][band]
    return phonons[qpoint][band]


def _weight_value(value: Any) -> int | float:
    """
    Convert q-point weights to stable YAML scalar values.

    Parameters
    ----------
    value : object
        Scalar weight value convertible to ``float``.

    Returns
    -------
    int or float
        Integer weight when possible, otherwise float.
    """
    number = float(value)
    if number.is_integer():
        return int(number)
    return number


def _as_float_list(values: Any) -> list[float]:
    """
    Convert an array-like object to a plain list of floats.

    Parameters
    ----------
    values : object
        Numeric sequence.

    Returns
    -------
    list of float
        Converted values.
    """
    return np.atleast_1d(np.asarray(values, dtype=np.float64)).astype(float).tolist()


def _as_int_matrix(values: Any) -> list[list[int]]:
    """
    Convert a matrix-like object to a nested list of integers.

    Parameters
    ----------
    values : object
        Matrix-like object.

    Returns
    -------
    list of list of int
        Converted matrix.
    """
    return np.asarray(values, dtype=int).tolist()


def _as_float(value: Any) -> float:
    """
    Convert a scalar-like value to float.

    Parameters
    ----------
    value : object
        Scalar value.

    Returns
    -------
    float
        Converted value.
    """
    return float(np.asarray(value, dtype=np.float64))
