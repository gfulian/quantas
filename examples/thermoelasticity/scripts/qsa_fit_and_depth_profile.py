"""Fit QSA elastic components and reconstruct C_ij(P,T) through the API.

The thermoelastic input is the readable YAML generated from CRYSTAL SOEC
outputs. The QHA source may be either a native Quantas HDF5 result or a QHA
YAML input. When a QHA YAML input is supplied, the QHA grid is calculated with
the options below before the thermoelastic fit is performed.

Examples
--------
Use a precomputed QHA result and add a linear geological profile::

    python scripts/qsa_fit_and_depth_profile.py dol_pbe0_thermoelastic.yaml \
        dol_qha.hdf5 --linear-profile --output dol_qsa.hdf5

Calculate QHA from YAML on a selected P-T grid and request an extended report::

    python scripts/qsa_fit_and_depth_profile.py dol_pbe0_thermoelastic.yaml \
        ../qha/crystal-phonons/dol_pbe0.yaml \
        --temperature-min 300 --temperature-max 1200 \
        --temperature-step 100 --pressure-min 0 --pressure-max 10 \
        --pressure-step 1 --report-level extended --output dol_qsa.hdf5

Use a tabulated profile with columns ``depth_km``, ``P_GPa``, and ``T_K``::

    python scripts/qsa_fit_and_depth_profile.py dol_pbe0_thermoelastic.yaml \
        dol_qha.hdf5 --profile-file mantle_profile.csv --output dol_qsa.hdf5
"""

from __future__ import annotations

import argparse
from pathlib import Path

from quantas.api import qha, rendering, thermoelasticity


def build_parser() -> argparse.ArgumentParser:
    """Return the command-line parser used by this API example."""
    parser = argparse.ArgumentParser(
        description=(
            "Fit independent cold-QSA elastic components and reconstruct "
            "isothermal C_ij(P,T)."
        )
    )
    parser.add_argument("thermoelastic_input", type=Path)
    parser.add_argument("qha_source", type=Path)
    parser.add_argument("--output", type=Path, default=Path("thermoelastic_qsa.hdf5"))
    parser.add_argument("--report", type=Path, default=Path("thermoelastic_qsa.dat"))
    parser.add_argument(
        "--report-level",
        choices=("standard", "extended", "debug"),
        default="standard",
    )

    qha = parser.add_argument_group("QHA grid used when qha_source is YAML")
    qha.add_argument("--temperature-min", type=float, default=300.0)
    qha.add_argument("--temperature-max", type=float, default=1200.0)
    qha.add_argument("--temperature-step", type=float, default=100.0)
    qha.add_argument("--pressure-min", type=float, default=0.0)
    qha.add_argument("--pressure-max", type=float, default=10.0)
    qha.add_argument("--pressure-step", type=float, default=1.0)

    profile = parser.add_argument_group("optional depth profile")
    profile.add_argument(
        "--profile-file",
        type=Path,
        help="CSV or whitespace table with depth_km, P_GPa, and T_K columns",
    )
    profile.add_argument(
        "--linear-profile",
        action="store_true",
        help="evaluate a profile defined by constant P and T gradients",
    )
    profile.add_argument("--profile-name", default="linear-geotherm")
    profile.add_argument("--depth-min", type=float, default=0.0)
    profile.add_argument("--depth-max", type=float, default=300.0)
    profile.add_argument("--depth-points", type=int, default=61)
    profile.add_argument("--pressure-at-depth-min", type=float, default=0.0)
    profile.add_argument("--pressure-gradient", type=float, default=0.03)
    profile.add_argument("--temperature-at-depth-min", type=float, default=300.0)
    profile.add_argument("--temperature-gradient", type=float, default=0.5)
    return parser


def profile_from_arguments(
    args: argparse.Namespace,
) -> list[thermoelasticity.DepthProfile]:
    """Build zero or one passive depth-profile contract from parsed arguments."""
    if args.profile_file is not None and args.linear_profile:
        raise ValueError("choose either --profile-file or --linear-profile")
    if args.profile_file is not None:
        return [thermoelasticity.read_depth_profile(args.profile_file)]
    if args.linear_profile:
        return [
            thermoelasticity.DepthProfile.linear(
                name=args.profile_name,
                depth_min=args.depth_min,
                depth_max=args.depth_max,
                npoints=args.depth_points,
                pressure_at_depth_min=args.pressure_at_depth_min,
                pressure_gradient=args.pressure_gradient,
                temperature_at_depth_min=args.temperature_at_depth_min,
                temperature_gradient=args.temperature_gradient,
            )
        ]
    return []


def main() -> None:
    """Run the API workflow, render the selected report, and write HDF5."""
    args = build_parser().parse_args()
    profiles = profile_from_arguments(args)

    qha_options = None
    if args.qha_source.suffix.lower() not in {".h5", ".hdf5"}:
        qha_options = qha.Options(
            temperature_min=args.temperature_min,
            temperature_max=args.temperature_max,
            temperature_step=args.temperature_step,
            pressure_min=args.pressure_min,
            pressure_max=args.pressure_max,
            pressure_step=args.pressure_step,
        )

    thermoelastic_options = thermoelasticity.Options(
        report_level=args.report_level,
        solver_debug=args.report_level == "debug",
    )
    result = thermoelasticity.run(
        args.thermoelastic_input,
        args.qha_source,
        options=thermoelastic_options,
        profiles=profiles,
        qha_options=qha_options,
    )

    report_text = rendering.render_tables(
        thermoelasticity.build_report(result, level=args.report_level)
    )
    print(report_text)
    args.report.write_text(report_text, encoding="utf-8")
    thermoelasticity.write_result(result, args.output, report_text=report_text)

    payload = thermoelasticity.get_result(result)
    print(f"HDF5 result: {args.output}")
    print(f"Text report: {args.report}")
    print(f"Independent components: {', '.join(payload.independent_labels)}")
    print(f"Completed: {payload.completed}")


if __name__ == "__main__":
    main()
