# -*- coding: utf-8 -*-

"""Rich terminal renderer for neutral Quantas report tables."""

from __future__ import annotations

from typing import Any, Literal

from rich import box
from rich.console import Group, RenderableType
from rich.table import Table
from rich.text import Text

from quantas.models import ReportTable
from quantas.renderers.tables.formatting import format_cell


def build_rich_table(table: ReportTable) -> Table:
    """Build one Rich table from a neutral Quantas report table.

    Parameters
    ----------
    table : ReportTable
        Frontend-neutral table specification.

    Returns
    -------
    rich.table.Table
        Terminal table using the terminal's default foreground color.  Only
        titles and headers are bold; no scientific value is color-coded.
    """
    column_formats = _normalized_metadata_list(
        table.metadata.get("column_formats"),
        len(table.columns),
        None,
    )
    column_units = _normalized_metadata_list(
        table.metadata.get("column_units"),
        len(table.columns),
        "",
    )
    alignments = _normalized_metadata_list(
        table.metadata.get("column_alignments"),
        len(table.columns),
        "left",
    )

    rendered = Table(
        title=Text(str(table.title), style="bold"),
        title_justify="left",
        box=box.SIMPLE_HEAVY,
        header_style="bold",
        show_edge=False,
        pad_edge=False,
        collapse_padding=True,
        expand=bool(table.metadata.get("expand", False)),
    )
    for index, column in enumerate(table.columns):
        header = str(column)
        unit = str(column_units[index]).strip()
        if unit:
            header = f"{header}\n({unit})"
        rendered.add_column(
            header,
            justify=_rich_justify(str(alignments[index])),
            vertical="middle",
            no_wrap=bool(table.metadata.get("no_wrap", False)),
        )

    for row in table.rows:
        cells = []
        for index in range(len(table.columns)):
            value = row[index] if index < len(row) else None
            cells.append(format_cell(value, column_formats[index]))
        rendered.add_row(*cells)
    return rendered


def build_rich_renderable(table: ReportTable) -> RenderableType:
    """Return a Rich table together with optional neutral notes."""
    renderables: list[RenderableType] = [build_rich_table(table)]
    for note in table.metadata.get("notes", []):
        renderables.append(Text(f"Note: {note}", style="dim"))
    return Group(*renderables)


def _normalized_metadata_list(
    value: Any,
    size: int,
    default: Any,
) -> list[Any]:
    """Return a metadata list with one entry per table column."""
    if not isinstance(value, (list, tuple)) or len(value) != size:
        return [default] * size
    return list(value)


def _rich_justify(
    value: str,
) -> Literal["default", "left", "center", "right", "full"]:
    """Map neutral alignment names to Rich column justification."""
    if value == "left":
        return "left"
    if value == "center":
        return "center"
    if value == "right":
        return "right"
    if value == "full":
        return "full"
    return "left"


__all__ = [
    "build_rich_renderable",
    "build_rich_table",
]
