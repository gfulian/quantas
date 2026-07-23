# -*- coding: utf-8 -*-

"""Shared command-line contracts for Quantas frontends.

This module centralizes user-facing option names, help groups, report
verbosity, and default output paths.  It contains no scientific logic and is
never imported by numerical modules.
"""

from __future__ import annotations

from enum import Enum
from pathlib import Path
from typing import Any, Callable, TypeVar, cast

import click

from quantas.renderers.plots import FIGURE_PRESET_NAMES

from .grouped_options import grouped_option

_F = TypeVar("_F", bound=Callable[..., Any])


class ReportVerbosity(str, Enum):
    """Supported levels of deterministic scientific reporting.

    ``standard`` provides the normal scientific summary, ``extended`` adds
    complete result tables where available, and ``debug`` additionally emits
    solver and per-state diagnostics.  The selected level never changes the
    numerical calculation or persisted numerical results.
    """

    STANDARD = "standard"
    EXTENDED = "extended"
    DEBUG = "debug"

    @property
    def includes_extended(self) -> bool:
        """Return whether extended scientific tables are requested."""
        return self in {self.EXTENDED, self.DEBUG}

    @property
    def includes_debug(self) -> bool:
        """Return whether detailed numerical diagnostics are requested."""
        return self is self.DEBUG


SCIENTIFIC_GROUP = "Scientific model"
DOMAIN_GROUP = "Calculation domain"
NUMERICAL_GROUP = "Numerical method"
VALIDATION_GROUP = "Validation and diagnostics"
UNITS_GROUP = "Units"
OUTPUT_GROUP = "Output and reporting"
PLOTTING_GROUP = "Plotting"
ADVANCED_GROUP = "Advanced numerical controls"


def parse_verbosity(value: str | ReportVerbosity) -> ReportVerbosity:
    """Normalize one CLI verbosity value.

    Parameters
    ----------
    value : str or ReportVerbosity
        User-facing verbosity value.

    Returns
    -------
    ReportVerbosity
        Normalized enumeration value.
    """
    if isinstance(value, ReportVerbosity):
        return value
    return ReportVerbosity(str(value).lower())


def default_report_path(source: Path, report: Path | None) -> Path:
    """Return the explicit report path or ``source`` with a ``.log`` suffix.

    Parameters
    ----------
    source : pathlib.Path
        Primary scientific input file.
    report : pathlib.Path or None
        Explicit report destination.

    Returns
    -------
    pathlib.Path
        Deterministic report destination.
    """
    return source.with_suffix(".log") if report is None else report


def default_hdf5_path(
    source: Path,
    output: Path | None,
    *,
    suffix: str,
) -> Path:
    """Return an explicit HDF5 path or a module-specific default.

    Parameters
    ----------
    source : pathlib.Path
        Primary input file.
    output : pathlib.Path or None
        Explicit output destination.
    suffix : str
        Module suffix inserted after the source stem.

    Returns
    -------
    pathlib.Path
        HDF5 destination with a normalized ``.hdf5`` extension.
    """
    if output is None:
        return source.with_name(f"{source.stem}{suffix}.hdf5")
    return output.with_suffix(".hdf5")


def output_option(*, help: str) -> Callable[[_F], _F]:
    """Return the standard ``-o/--output`` path decorator."""
    return cast(
        Callable[[_F], _F],
        grouped_option(
            "-o",
            "--output",
            group=OUTPUT_GROUP,
            type=click.Path(dir_okay=False, path_type=Path),
            default=None,
            help=help,
        ),
    )


def report_option() -> Callable[[_F], _F]:
    """Return the standard automatic report option decorator."""
    return cast(
        Callable[[_F], _F],
        grouped_option(
            "-r",
            "--report",
            group=OUTPUT_GROUP,
            type=click.Path(dir_okay=False, path_type=Path),
            default=None,
            help=(
                "Deterministic plain-text report. Default: primary input "
                "file name with a '.log' extension."
            ),
        ),
    )


def verbosity_option() -> Callable[[_F], _F]:
    """Return the standard ``-v/--verbosity`` decorator."""
    return cast(
        Callable[[_F], _F],
        grouped_option(
            "-v",
            "--verbosity",
            group=OUTPUT_GROUP,
            type=click.Choice(
                [item.value for item in ReportVerbosity],
                case_sensitive=False,
            ),
            default=ReportVerbosity.STANDARD.value,
            show_default=True,
            help=(
                "Report detail: standard scientific summary, extended complete "
                "tables, or debug numerical diagnostics."
            ),
        ),
    )


def quiet_option() -> Callable[[_F], _F]:
    """Return the standard terminal-suppression decorator."""
    return cast(
        Callable[[_F], _F],
        grouped_option(
            "-q",
            "--quiet",
            group=OUTPUT_GROUP,
            is_flag=True,
            default=False,
            help=(
                "Suppress terminal presentation while still writing the report "
                "and structured result files."
            ),
        ),
    )


def progress_option() -> Callable[[_F], _F]:
    """Return the standard transient-progress decorator."""
    return cast(
        Callable[[_F], _F],
        grouped_option(
            "--progress/--no-progress",
            group=OUTPUT_GROUP,
            default=True,
            show_default=True,
            help="Display transient terminal progress; progress is never persisted.",
        ),
    )


def figure_preset_option(
    *,
    option_name: str = "--preset",
    parameter_name: str = "figure_preset",
    group: str = PLOTTING_GROUP,
) -> Callable[[_F], _F]:
    """Return the standard static-figure preset decorator.

    Parameters
    ----------
    option_name : str, optional
        Long option name. Commands that already use ``--preset`` for a
        scientific profile may request ``--plot-preset`` instead.
    parameter_name : str, optional
        Python parameter receiving the selected preset.
    group : str, optional
        Help group used by :class:`~quantas.cli.grouped_options.GroupedCommand`.

    Returns
    -------
    callable
        Click decorator exposing ``screen``, ``publication``, and
        ``monochrome``.
    """

    return cast(
        Callable[[_F], _F],
        grouped_option(
            option_name,
            parameter_name,
            group=group,
            type=click.Choice(FIGURE_PRESET_NAMES, case_sensitive=False),
            default="screen",
            show_default=True,
            help=(
                "Standard figure styling. Screen favors interactive output; "
                "publication uses print-ready geometry and resolution; "
                "monochrome adds grayscale-safe rendering."
            ),
        ),
    )


def force_option() -> Callable[[_F], _F]:
    """Return the standard overwrite-policy decorator."""
    return cast(
        Callable[[_F], _F],
        grouped_option(
            "-f",
            "--force",
            group=OUTPUT_GROUP,
            is_flag=True,
            default=False,
            help="Replace existing generated files without prompting.",
        ),
    )


__all__ = [
    "ADVANCED_GROUP",
    "DOMAIN_GROUP",
    "NUMERICAL_GROUP",
    "OUTPUT_GROUP",
    "PLOTTING_GROUP",
    "ReportVerbosity",
    "SCIENTIFIC_GROUP",
    "UNITS_GROUP",
    "VALIDATION_GROUP",
    "default_hdf5_path",
    "default_report_path",
    "figure_preset_option",
    "force_option",
    "output_option",
    "parse_verbosity",
    "progress_option",
    "quiet_option",
    "report_option",
    "verbosity_option",
]
