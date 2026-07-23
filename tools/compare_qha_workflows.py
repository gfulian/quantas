#!/usr/bin/env python3
"""Validate the four Quantas QHA workflow combinations on one input file."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from quantas.modules.qha.calculator import QHACalculator
from quantas.modules.qha.io.reader import read_qha_input
from quantas.modules.qha.models import QHAOptions, QHAResult
from quantas.modules.qha.validation import compare_qha_results, validate_qha_result


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed validation options.
    """
    parser = argparse.ArgumentParser(
        description="Run and compare freq/td x poly/eos QHA workflows."
    )
    parser.add_argument("input", type=Path, help="QHA YAML input file")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("qha_workflow_validation"),
        help="Output path without suffix",
    )
    parser.add_argument("--temperature-min", type=float, default=0.0)
    parser.add_argument("--temperature-max", type=float, default=2000.0)
    parser.add_argument("--temperature-step", type=float, default=10.0)
    parser.add_argument("--pressure-min", type=float, default=0.0)
    parser.add_argument("--pressure-max", type=float, default=30.0)
    parser.add_argument("--pressure-step", type=float, default=5.0)
    parser.add_argument("--eos", default="BM")
    parser.add_argument("-E", "--energy-degree", type=int, default=3)
    parser.add_argument("--free-energy-degree", type=int, default=3)
    parser.add_argument("-D", "--frequency-degree", type=int, default=3)
    return parser.parse_args()


def run_workflows(arguments: argparse.Namespace) -> tuple[Any, dict[str, QHAResult]]:
    """Run all four QHA workflow combinations.

    Parameters
    ----------
    arguments : argparse.Namespace
        Parsed validation settings.

    Returns
    -------
    tuple
        Input data and a mapping of workflow names to QHA results.
    """
    input_data = read_qha_input(arguments.input)
    results: dict[str, QHAResult] = {}
    for scheme in ("freq", "td"):
        for minimization in ("poly", "eos"):
            key = f"{scheme}_{minimization}"
            options = QHAOptions(
                temperature_min=arguments.temperature_min,
                temperature_max=arguments.temperature_max,
                temperature_step=arguments.temperature_step,
                pressure_min=arguments.pressure_min,
                pressure_max=arguments.pressure_max,
                pressure_step=arguments.pressure_step,
                scheme=scheme,
                minimization=minimization,
                eos=arguments.eos,
                energy_degree=arguments.energy_degree,
                free_energy_degree=arguments.free_energy_degree,
                frequency_degree=arguments.frequency_degree,
                polynomial_derivative_method="local_grid",
                estimate_uncertainties=minimization == "eos",
                uncertainty_method="covariance",
            )
            results[key] = QHACalculator(input_data, options).execute().results["qha"]
    return input_data, results


def build_validation_data(
    input_data: Any, results: dict[str, QHAResult]
) -> dict[str, Any]:
    """Build serializable validation summaries and pair comparisons.

    Parameters
    ----------
    input_data : QHAInput
        QHA input used by all workflows.
    results : dict
        Workflow result mapping.

    Returns
    -------
    dict
        Validation data suitable for JSON output.
    """
    summaries = {
        name: validate_qha_result(result, input_data).as_dict()
        for name, result in results.items()
    }
    pair_names = (
        ("freq_poly", "td_poly"),
        ("freq_eos", "td_eos"),
        ("freq_poly", "freq_eos"),
        ("td_poly", "td_eos"),
    )
    comparisons: dict[str, list[dict[str, Any]]] = {}
    for first, second in pair_names:
        comparisons[f"{first}__{second}"] = [
            {
                "property": item.property_name,
                "maximum_absolute": item.maximum_absolute,
                "maximum_relative": item.maximum_relative,
                "root_mean_square": item.root_mean_square,
                "compared_points": item.compared_points,
                "maximum_absolute_temperature": item.maximum_absolute_temperature,
                "maximum_absolute_pressure": item.maximum_absolute_pressure,
                "maximum_relative_temperature": item.maximum_relative_temperature,
                "maximum_relative_pressure": item.maximum_relative_pressure,
            }
            for item in compare_qha_results(results[first], results[second])
        ]
    return {"summaries": summaries, "comparisons": comparisons}


def render_markdown(data: dict[str, Any]) -> str:
    """Render validation data as a Markdown report.

    Parameters
    ----------
    data : dict
        Validation data returned by :func:`build_validation_data`.

    Returns
    -------
    str
        Markdown report.
    """
    lines = [
        "# QHA four-workflow validation",
        "",
        "## Workflow checks",
        "",
        "| Workflow | Valid | Finite | V(P) monotonic | K > 0 | 0 K identities | Cp >= Cv | Cv/DP | Outside V range |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for name, summary in data["summaries"].items():
        lines.append(
            "| {name} | {valid}/{total} | {finite} | {monotonic} | {positive} | "
            "{zero} | {cp} | {dp:.6f} | {outside} |".format(
                name=name,
                valid=summary["valid_points"],
                total=summary["total_points"],
                finite=_yes_no(summary["finite_properties"]),
                monotonic=_yes_no(summary["volume_decreases_with_pressure"]),
                positive=_yes_no(summary["positive_bulk_moduli"]),
                zero=_yes_no(summary["zero_kelvin_consistency"]),
                cp=_yes_no(summary["cp_not_below_cv"]),
                dp=float(summary["dulong_petit_ratio"] or float("nan")),
                outside=(
                    summary["volumes_below_sampled_range"]
                    + summary["volumes_above_sampled_range"]
                ),
            )
        )
        if summary["warnings"]:
            lines.append("")
            lines.append(f"Warnings for `{name}`:")
            lines.extend(f"- {warning}" for warning in summary["warnings"])
            lines.append("")

    for pair, differences in data["comparisons"].items():
        lines.extend(
            [
                "",
                f"## Comparison: `{pair.replace('__', ' vs ')}`",
                "",
                "| Property | Max abs | At (T, P) | Max relative | At (T, P) | RMS | Points |",
                "|---|---:|---:|---:|---:|---:|---:|",
            ]
        )
        for item in differences:
            lines.append(
                "| {property} | {maximum_absolute:.8E} | "
                "({maximum_absolute_temperature:.2f}, {maximum_absolute_pressure:.2f}) | "
                "{maximum_relative:.8E} | "
                "({maximum_relative_temperature:.2f}, {maximum_relative_pressure:.2f}) | "
                "{root_mean_square:.8E} | {compared_points} |".format(**item)
            )
    lines.append("")
    return "\n".join(lines)


def _yes_no(value: bool | None) -> str:
    """Return a compact Markdown label for a boolean value."""
    if value is None:
        return "n/a"
    return "yes" if value else "no"


def main() -> None:
    """Run validation and write JSON and Markdown reports."""
    arguments = parse_arguments()
    input_data, results = run_workflows(arguments)
    data = build_validation_data(input_data, results)
    base = arguments.output.with_suffix("")
    json_path = base.with_suffix(".json")
    markdown_path = base.with_suffix(".md")
    json_path.write_text(json.dumps(data, indent=2), encoding="utf8")
    markdown_path.write_text(render_markdown(data), encoding="utf8")
    print(markdown_path)
    print(json_path)


if __name__ == "__main__":
    main()
