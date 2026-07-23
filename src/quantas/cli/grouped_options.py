# -*- coding: utf-8 -*-

"""Shared Click option-group presentation helpers.

The grouping logic is frontend-only. Scientific modules never import this
module; it is shared by CLI command families that need readable help output.
"""

from __future__ import annotations

from collections import OrderedDict
from collections.abc import Callable
from typing import Any

import click


class GroupedOption(click.Option):
    """Click option carrying a logical help-section label."""

    def __init__(self, *args: Any, help_group: str = "Options", **kwargs: Any) -> None:
        self.help_group = str(help_group)
        super().__init__(*args, **kwargs)


class GroupedCommand(click.Command):
    """Click command that renders options under stable logical headings."""

    def format_options(
        self, ctx: click.Context, formatter: click.HelpFormatter
    ) -> None:
        """Render command options under canonical semantic headings."""
        groups: OrderedDict[str, list[tuple[str, str]]] = OrderedDict()
        for parameter in self.get_params(ctx):
            record = parameter.get_help_record(ctx)
            if record is None:
                continue
            group = getattr(parameter, "help_group", None)
            if group is None:
                group = _infer_option_group(parameter, self.name or "")
            groups.setdefault(group, []).append(record)
        ordered_groups = sorted(
            groups.items(),
            key=lambda item: _help_group_priority(item[0]),
        )
        for title, records in ordered_groups:
            with formatter.section(title):
                formatter.write_dl(records)


def _infer_option_group(parameter: click.Parameter, command_name: str) -> str:
    """Infer a sensible group for ungrouped Click options.

    New or complex commands should still use :func:`grouped_option`
    explicitly.  This fallback keeps small utility commands consistent while
    they are migrated without embedding scientific behavior in the CLI layer.
    """
    if not isinstance(parameter, click.Option):
        return "Options"
    names = {name.lstrip("-").replace("-", "_") for name in parameter.opts}
    if names & {"help"}:
        return "Options"
    if names & {
        "output",
        "outfile",
        "outbase",
        "report",
        "quiet",
        "force",
        "show",
    }:
        return "Output and reporting"
    if any(
        name.endswith("unit") or name in {"eunit", "vunit", "funit", "tunit", "punit"}
        for name in names
    ):
        return "Units"
    if command_name == "plot":
        return "Plotting"
    if names & {
        "temperature",
        "pressure",
        "ntheta",
        "nphi",
        "2d",
        "3d",
        "depth",
    }:
        return "Calculation domain"
    if names & {"dpi", "format"}:
        return "Output and reporting"
    return "Scientific selection"


def _help_group_priority(title: str) -> int:
    """Return the canonical presentation priority for one option group.

    The classification intentionally uses semantic group names rather than
    command-specific registration order.  Stable sorting preserves the local
    order of groups that belong to the same broad category.
    """
    normalized = title.strip().lower()
    if normalized == "scientific model" or any(
        token in normalized
        for token in (
            "model",
            "coupling",
            "normalization",
            "adiabatic",
            "uncertainty propagation",
        )
    ):
        return 10
    if any(
        token in normalized
        for token in (
            "domain",
            "selection",
            "coordinate",
            "profile definition",
            "synthetic profile",
            "yaml evaluation",
        )
    ):
        return 20
    if any(
        token in normalized
        for token in (
            "numerical",
            "solver",
            "constraint",
            "effective variance",
            "odr",
        )
    ):
        return 30
    if any(
        token in normalized
        for token in (
            "validation",
            "diagnostic",
            "stability",
            "support",
            "batch execution",
        )
    ):
        return 40
    if "unit" in normalized:
        return 50
    if any(
        token in normalized
        for token in (
            "plot",
            "figure",
            "style",
            "line",
            "marker",
            "axes",
            "contour",
            "presentation",
            "annotation",
            "geometry",
            "typography",
        )
    ):
        return 60
    if any(
        token in normalized
        for token in (
            "output",
            "report",
            "persistence",
            "portable table",
        )
    ):
        return 70
    if "advanced" in normalized:
        return 80
    return 90


def grouped_option(*param_decls: str, group: str, **attrs: Any) -> Callable[[Any], Any]:
    """Return a Click option decorator assigned to one help group."""
    return click.option(*param_decls, cls=GroupedOption, help_group=group, **attrs)


__all__ = ["GroupedCommand", "GroupedOption", "grouped_option"]
