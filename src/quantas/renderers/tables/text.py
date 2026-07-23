# -*- coding: utf-8 -*-
"""Plain-text renderer for neutral Quantas report tables."""

from __future__ import annotations

from collections.abc import Sequence

from quantas.models import ReportTable
from quantas.renderers.tables.formatting import format_cell


def render_table(table: ReportTable) -> str:
    """Render one neutral report table as aligned plain text.

    Multiline cells are rendered as physical table rows and column widths are
    calculated from the longest visible line, not from the total string length.

    Parameters
    ----------
    table : ReportTable
        Neutral table specification.

    Returns
    -------
    str
        Plain-text table suitable for terminals and report files.
    """
    columns = [format_cell(column) for column in table.columns]
    column_formats = table.metadata.get("column_formats", [])
    if len(column_formats) != len(columns):
        column_formats = [None] * len(columns)
    rows = [
        [
            format_cell(
                row[index] if index < len(row) else None,
                column_formats[index],
            )
            for index in range(len(columns))
        ]
        for row in table.rows
    ]
    units = [
        (f"({format_cell(unit)})" if str(unit).strip() else "")
        for unit in table.metadata.get("column_units", [])
    ]
    if units and len(units) != len(columns):
        units = []

    widths = _column_widths(columns, rows, units=units)
    alignments = table.metadata.get("column_alignments", ["left"] * len(columns))
    if len(alignments) != len(columns):
        alignments = ["left"] * len(columns)

    title = str(table.title)
    lines = [title, "-" * len(title), ""]
    if columns:
        lines.extend(_format_row(columns, widths, ["center"] * len(columns)))
        if units:
            lines.extend(_format_row(units, widths, ["center"] * len(columns)))
        lines.append(_format_separator(widths))
    for row in rows:
        lines.extend(_format_row(row, widths, alignments))

    notes = table.metadata.get("notes", [])
    if notes:
        lines.append("")
        lines.extend(f"  Note: {note}" for note in notes)
    lines.append("")
    return "\n".join(lines)


def render_tables(tables: Sequence[ReportTable]) -> str:
    """Render an ordered collection of neutral report tables.

    Parameters
    ----------
    tables : sequence of ReportTable
        Tables to render.

    Returns
    -------
    str
        Concatenated plain-text report.
    """
    return "\n".join(render_table(table) for table in tables)


def _column_widths(
    columns: list[str],
    rows: list[list[str]],
    *,
    units: list[str] | None = None,
) -> list[int]:
    """Return visible display widths for all table columns."""
    widths = [_visible_width(column) for column in columns]
    if units:
        for index, unit in enumerate(units):
            widths[index] = max(widths[index], _visible_width(unit))
    for row in rows:
        for index, cell in enumerate(row):
            widths[index] = max(widths[index], _visible_width(cell))
    return widths


def _format_row(
    row: list[str],
    widths: list[int],
    alignments: list[str] | None = None,
) -> list[str]:
    """Format one logical row, including multiline cells."""
    if alignments is None:
        alignments = ["left"] * len(row)
    split_cells = [_cell_lines(cell) for cell in row]
    height = max((len(lines) for lines in split_cells), default=1)
    rendered: list[str] = []
    for line_index in range(height):
        cells: list[str] = []
        for index, lines in enumerate(split_cells):
            line = lines[line_index] if line_index < len(lines) else ""
            alignment = alignments[index]
            if alignment == "right":
                cells.append(line.rjust(widths[index]))
            elif alignment == "center":
                cells.append(line.center(widths[index]))
            else:
                cells.append(line.ljust(widths[index]))
        rendered.append(("  " + " | ".join(cells)).rstrip())
    return rendered


def _format_separator(widths: list[int]) -> str:
    """Format the table header separator."""
    return "  " + "-+-".join("-" * width for width in widths)


def _cell_lines(value: str) -> list[str]:
    """Split a rendered cell into visible lines without losing empty cells."""
    lines = str(value).splitlines()
    return lines if lines else [""]


def _visible_width(value: str) -> int:
    """Return the width of the longest visible line in a rendered cell."""
    return max((len(line) for line in _cell_lines(value)), default=0)
