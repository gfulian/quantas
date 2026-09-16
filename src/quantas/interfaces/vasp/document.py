# -*- coding: utf-8 -*-

"""Source resolution and XML document access for one VASP calculation."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
from pathlib import Path
from typing import Any
import xml.etree.ElementTree as ET

import numpy as np
from numpy.typing import NDArray

from quantas.core.chemistry.symbols import symbol2number


_VASPRUN_NAME = "vasprun.xml"
_OUTCAR_NAME = "OUTCAR"


@dataclass(frozen=True, slots=True)
class VaspRunSource:
    """Resolved files belonging to one VASP calculation.

    Parameters
    ----------
    directory : pathlib.Path
        Calculation directory.
    vasprun_xml : pathlib.Path
        Primary structured VASP output.
    outcar : pathlib.Path or None
        Sibling human-readable output when available.
    """

    directory: Path
    vasprun_xml: Path
    outcar: Path | None = None

    def __post_init__(self) -> None:
        """Normalize stored paths."""
        object.__setattr__(self, "directory", Path(self.directory))
        object.__setattr__(self, "vasprun_xml", Path(self.vasprun_xml))
        if self.outcar is not None:
            object.__setattr__(self, "outcar", Path(self.outcar))




def resolve_vasp_run_source(source: str | Path) -> VaspRunSource:
    """Resolve a VASP calculation directory or one of its primary outputs.

    Parameters
    ----------
    source : str or pathlib.Path
        Calculation directory, ``vasprun.xml``, or ``OUTCAR`` path.

    Returns
    -------
    VaspRunSource
        Resolved calculation directory and available primary outputs.

    Raises
    ------
    FileNotFoundError
        If the requested path or required sibling ``vasprun.xml`` is absent.
    ValueError
        If a file other than ``vasprun.xml`` or ``OUTCAR`` is supplied.
    """
    path = Path(source)
    if not path.exists():
        raise FileNotFoundError(f"VASP run source does not exist: {path}")

    if path.is_dir():
        directory = path
        vasprun = directory / _VASPRUN_NAME
        outcar = directory / _OUTCAR_NAME
    else:
        if path.name.casefold() == _VASPRUN_NAME.casefold():
            directory = path.parent
            vasprun = path
            outcar = directory / _OUTCAR_NAME
        elif path.name == _OUTCAR_NAME:
            directory = path.parent
            vasprun = directory / _VASPRUN_NAME
            outcar = path
        else:
            raise ValueError(
                "VASP run source must be a calculation directory, vasprun.xml, "
                "or OUTCAR"
            )

    if not vasprun.is_file():
        raise FileNotFoundError(
            f"VASP run source requires '{_VASPRUN_NAME}': {vasprun}"
        )
    return VaspRunSource(
        directory=directory,
        vasprun_xml=vasprun,
        outcar=outcar if outcar.is_file() else None,
    )


class VaspRunDocument:
    """Parsed primary documents for one VASP calculation directory.

    The class owns code-specific XML/text parsing and exposes stable access to
    source metadata and ionic-step containers.  Scientific selection policy is
    intentionally left to
    :class:`quantas.interfaces.vasp.output.VaspOutputParser` and later module
    adapters.

    Parameters
    ----------
    source : str or pathlib.Path
        VASP calculation directory, ``vasprun.xml``, or sibling ``OUTCAR``.

    Raises
    ------
    FileNotFoundError
        If the run source cannot be resolved to ``vasprun.xml``.
    ValueError
        If the XML document is malformed or is not a VASP ``modeling`` record.
    """

    def __init__(self, source: str | Path) -> None:
        self.source = resolve_vasp_run_source(source)
        try:
            self.tree = ET.parse(self.source.vasprun_xml)
        except ET.ParseError as exc:
            raise ValueError(
                f"malformed VASP XML document: {self.source.vasprun_xml}"
            ) from exc
        self.root = self.tree.getroot()
        if self.root.tag != "modeling":
            raise ValueError(
                "VASP vasprun.xml root element must be 'modeling', "
                f"observed {self.root.tag!r}"
            )
        self.outcar_text = (
            self.source.outcar.read_text(encoding="utf-8", errors="replace")
            if self.source.outcar is not None
            else None
        )

    def generator(self) -> dict[str, str]:
        """Return VASP generator metadata from the XML document.

        Returns
        -------
        dict
            Named string fields such as program, version, build, and platform.
        """
        node = self.root.find("generator")
        if node is None:
            return {}
        return {
            str(item.get("name")): _text(item)
            for item in node.findall("i")
            if item.get("name")
        }

    def incar_parameters(self) -> dict[str, Any]:
        """Return explicitly recorded INCAR scalar values.

        Returns
        -------
        dict
            Scalar INCAR values with XML integer/logical/string types retained
            where VASP declares them.
        """
        node = self.root.find("incar")
        if node is None:
            return {}
        return {
            str(item.get("name")): _typed_xml_value(item)
            for item in node.findall("i")
            if item.get("name")
        }

    def parameter(self, name: str, default: Any = None) -> Any:
        """Return one effective VASP scalar parameter when available.

        ``parameters`` contains VASP's effective values, including defaults;
        the explicitly supplied ``incar`` value is used as a fallback.

        Parameters
        ----------
        name : str
            Exact VASP parameter name, for example ``"ISMEAR"`` or ``"ISIF"``.
        default : Any, optional
            Value returned when neither XML section contains the parameter.

        Returns
        -------
        Any
            Parsed scalar parameter or ``default``.
        """
        parameters = self.root.find("parameters")
        if parameters is not None:
            for item in parameters.iter("i"):
                if item.get("name") == name:
                    return _typed_xml_value(item)
        incar = self.root.find("incar")
        if incar is not None:
            item = incar.find(f"i[@name='{name}']")
            if item is not None:
                return _typed_xml_value(item)
        return default

    def parameter_vector(self, name: str) -> tuple[Any, ...] | None:
        """Return one effective VASP vector parameter when available.

        Parameters
        ----------
        name : str
            Exact VASP vector parameter name.

        Returns
        -------
        tuple or None
            Parsed vector values, or ``None`` when absent.
        """
        return _parameter_vector(self.root, name)

    def kpoint_signature(self) -> tuple[str, ...]:
        """Return a stable description of the Brillouin-zone sampling.

        Returns
        -------
        tuple of str
            Generated-mesh fields or a digest of an explicit k-point list.
        """
        return _kpoint_signature(self.root)

    def atom_symbols(self) -> tuple[str, ...]:
        """Return atom symbols in the exact VASP atom order.

        Returns
        -------
        tuple of str
            Element symbols for all atoms.

        Raises
        ------
        ValueError
            If the ``atominfo`` table is absent or malformed.
        """
        atominfo = self.root.find("atominfo")
        array = None if atominfo is None else atominfo.find("array[@name='atoms']")
        rows = None if array is None else array.find("set")
        if rows is None:
            raise ValueError("VASP atominfo/atoms table is unavailable")
        symbols: list[str] = []
        for row in rows:
            cells = row.findall("c")
            if not cells:
                continue
            symbol = _text(cells[0])
            if not symbol:
                raise ValueError("VASP atominfo contains an empty element symbol")
            try:
                symbol2number(symbol)
            except ValueError as exc:
                raise ValueError(
                    f"VASP atominfo contains an unsupported element symbol: {symbol!r}"
                ) from exc
            symbols.append(symbol.strip().capitalize())
        if not symbols:
            raise ValueError("VASP atominfo/atoms table contains no atoms")
        return tuple(symbols)

    def pseudopotential_labels(self) -> tuple[str, ...]:
        """Return VASP pseudopotential labels recorded in ``atominfo``.

        Returns
        -------
        tuple of str
            One label per VASP atom type, for example ``"PAW_PBE Mg_pv
            13Apr2007"`` as recorded by ``vasprun.xml``.  Empty tuple is
            returned when the atom-type table does not expose labels.
        """
        atominfo = self.root.find("atominfo")
        array = None if atominfo is None else atominfo.find("array[@name='atomtypes']")
        rows = None if array is None else array.find("set")
        if rows is None:
            return ()
        labels: list[str] = []
        for row in rows.findall("rc"):
            cells = row.findall("c")
            if len(cells) < 5:
                continue
            label = _text(cells[4])
            if label:
                labels.append(label)
        return tuple(labels)

    def atomic_numbers(self) -> NDArray[np.int64]:
        """Return atomic numbers in VASP atom order.

        Returns
        -------
        numpy.ndarray
            One-dimensional ``int64`` atomic-number array.
        """
        return np.asarray(
            [symbol2number(symbol) for symbol in self.atom_symbols()],
            dtype=np.int64,
        )

    def structure_node(self, name: str) -> ET.Element:
        """Return a named top-level VASP structure node.

        Parameters
        ----------
        name : str
            VASP structure name, normally ``"initialpos"`` or ``"finalpos"``.

        Returns
        -------
        xml.etree.ElementTree.Element
            Matching XML node.

        Raises
        ------
        ValueError
            If the requested named structure is absent.
        """
        node = self.root.find(f"structure[@name='{name}']")
        if node is None:
            raise ValueError(f"VASP structure {name!r} is unavailable")
        return node

    def ionic_step_nodes(self) -> tuple[ET.Element, ...]:
        """Return normalized XML containers for ionic states.

        Older VASP XML, including the validated 5.4.4 fixture, encloses each
        ionic state in ``<calculation>``.  Current VASP documentation describes
        flat unnamed ``<structure>`` blocks under ``<modeling>``.  This method
        hides that layout difference from the semantic parser.

        Returns
        -------
        tuple of xml.etree.ElementTree.Element
            One synthetic or original container per ionic state.
        """
        calculations = tuple(self.root.findall("calculation"))
        if calculations:
            return calculations

        children = list(self.root)
        structures = [
            index
            for index, node in enumerate(children)
            if node.tag == "structure" and not node.get("name")
        ]
        containers: list[ET.Element] = []
        for position, start in enumerate(structures):
            stop = structures[position + 1] if position + 1 < len(structures) else len(children)
            container = ET.Element("ionic_step")
            for node in children[start:stop]:
                if node.tag == "structure" and node.get("name"):
                    break
                container.append(node)
            containers.append(container)
        return tuple(containers)



def _parameter_vector(root: ET.Element, name: str) -> tuple[Any, ...] | None:
    """Return one effective vector parameter from parameters or INCAR."""
    parameters = root.find("parameters")
    if parameters is not None:
        for item in parameters.iter("v"):
            if item.get("name") == name:
                return _typed_xml_vector(item)
    incar = root.find("incar")
    if incar is not None:
        item = incar.find(f"v[@name='{name}']")
        if item is not None:
            return _typed_xml_vector(item)
    return None


def _kpoint_signature(root: ET.Element) -> tuple[str, ...]:
    """Return generated-mesh fields or an explicit-list digest."""
    node = root.find("kpoints")
    if node is None:
        return ()
    generation = node.find("generation")
    if generation is not None:
        fields = [f"mode={generation.get('param', '').strip() or 'unknown'}"]
        for item in generation.findall("v"):
            name = item.get("name")
            if not name:
                continue
            rendered = ",".join(
                _stable_value(value) for value in _typed_xml_vector(item)
            )
            fields.append(f"{name}={rendered}")
        return tuple(fields)

    points = node.find("varray[@name='kpointlist']")
    weights = node.find("varray[@name='weights']")
    if points is None:
        return ()
    point_rows = tuple(_numeric_row(item) for item in points.findall("v"))
    weight_rows = (
        tuple(_numeric_row(item) for item in weights.findall("v"))
        if weights is not None
        else ()
    )
    payload = repr((point_rows, weight_rows)).encode("ascii")
    digest = hashlib.sha256(payload).hexdigest()
    return ("mode=explicit", f"nkpoints={len(point_rows)}", f"sha256={digest}")


def _text(node: ET.Element) -> str:
    """Return stripped XML element text."""
    return "" if node.text is None else node.text.strip()


def _as_float(value: str) -> float:
    """Parse a VASP floating-point token including Fortran ``D`` exponents."""
    return float(value.strip().replace("D", "E").replace("d", "e"))


def _typed_xml_value(node: ET.Element) -> Any:
    """Parse one scalar VASP XML value according to its declared type."""
    value = _text(node)
    kind = (node.get("type") or "").casefold()
    if kind == "string":
        return value
    if kind == "int":
        return int(value)
    if kind == "logical":
        return value.casefold() in {"t", ".true.", "true"}
    try:
        return _as_float(value)
    except ValueError:
        return value


def _typed_xml_vector(node: ET.Element) -> tuple[Any, ...]:
    """Parse one VASP XML vector according to its declared scalar type."""
    tokens = _text(node).split()
    kind = (node.get("type") or "").casefold()
    if kind == "int":
        return tuple(int(token) for token in tokens)
    if kind == "logical":
        return tuple(token.casefold() in {"t", ".true.", "true"} for token in tokens)
    if kind == "string":
        return tuple(tokens)
    values: list[Any] = []
    for token in tokens:
        try:
            values.append(_as_float(token))
        except ValueError:
            values.append(token)
    return tuple(values)


def _stable_value(value: Any) -> str:
    """Render one XML value deterministically for provenance signatures."""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (float, np.floating)):
        return format(float(value), ".15g")
    return str(value).strip()


def _numeric_row(node: ET.Element) -> tuple[float, ...]:
    """Return one numerical VASP varray row for deterministic hashing."""
    return tuple(_as_float(token) for token in _text(node).split())
