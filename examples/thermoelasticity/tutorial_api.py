"""Run the complete dolomite thermoelasticity tutorial through the public API.

The script expects the normalized thermoelastic input, the dolomite QHA YAML,
and the composed continental profile distributed with Quantas.  It performs
QHA, calibrates the quasi-static elastic model, evaluates a point, a P--T grid,
and a depth profile, writes HDF5/table products, and renders selected figures.

Examples
--------
Run from the project root::

    python examples/thermoelasticity/tutorial_api.py \
        --output-dir thermoelasticity_tutorial
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from quantas.api import qha, rendering, thermoelasticity


def parser() -> argparse.ArgumentParser:
    """Return the command-line parser for the standalone API tutorial."""
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument(
        "--thermoelastic-input",
        type=Path,
        default=Path("examples/thermoelasticity/dol_pbe0_thermoelastic.yaml"),
    )
    result.add_argument(
        "--qha-input",
        type=Path,
        default=Path("examples/qha/crystal-phonons/dol_pbe0.yaml"),
    )
    result.add_argument(
        "--profile-spec",
        type=Path,
        default=Path(
            "examples/thermoelasticity/dolomite_continental_profile.yaml"
        ),
    )
    result.add_argument(
        "--output-dir",
        type=Path,
        default=Path("thermoelasticity_tutorial"),
    )
    return result


def main() -> None:
    """Execute the complete API workflow and write tutorial products."""
    args = parser().parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    figures = args.output_dir / "figures"
    figures.mkdir(parents=True, exist_ok=True)

    qha_result = qha.run(
        args.qha_input,
        options=qha.Options(
            temperature_min=300.0,
            temperature_max=1800.0,
            temperature_step=100.0,
            pressure_min=0.0,
            pressure_max=10.0,
            pressure_step=2.0,
            scheme="freq",
            minimization="poly",
            energy_degree=3,
            frequency_degree=3,
            thermal_expansion_method="mixed_derivative",
        ),
    )
    qha_report = rendering.render_tables(qha.build_report(qha_result))
    qha.write_result(
        qha_result,
        args.output_dir / "dolomite_qha.hdf5",
        report_text=qha_report,
    )
    (args.output_dir / "dolomite_qha.log").write_text(
        qha_report,
        encoding="utf-8",
    )

    fit_result = thermoelasticity.run(
        args.thermoelastic_input,
        qha_result,
        options=thermoelasticity.Options(
            reference_eos="BM3",
            finite_strain_order=3,
            adiabatic_mode="auto",
            quality_policy="warn",
            report_level="standard",
        ),
    )
    fit_report = rendering.render_tables(
        thermoelasticity.build_report(fit_result, level="standard")
    )
    thermoelasticity.write_result(
        fit_result,
        args.output_dir / "dolomite_fit.hdf5",
        report_text=fit_report,
    )
    (args.output_dir / "dolomite_fit.log").write_text(
        fit_report,
        encoding="utf-8",
    )

    point_result = thermoelasticity.analyze_grid(
        fit_result,
        pressure=np.asarray([5.0], dtype=np.float64),
        temperature=np.asarray([800.0], dtype=np.float64),
        extrapolation_policy="fail",
    )
    point_payload = thermoelasticity.get_result(point_result)
    print(
        rendering.render_tables(
            [thermoelasticity.point_table(point_payload, tensor_condition="adiabatic")]
        )
    )
    thermoelasticity.write_state_input(
        point_payload,
        args.output_dir / "dolomite_5GPa_800K.dat",
        tensor_condition="adiabatic",
        overwrite=True,
    )

    grid_result = thermoelasticity.analyze_grid(
        fit_result,
        pressure=thermoelasticity.regular_grid(0.0, 10.0, 2.0),
        temperature=thermoelasticity.regular_grid(300.0, 1200.0, 100.0),
        extrapolation_policy="fail",
    )
    grid_payload = thermoelasticity.get_result(grid_result)
    thermoelasticity.write_result(
        grid_result,
        args.output_dir / "dolomite_grid.hdf5",
    )
    thermoelasticity.write_grid_table(
        grid_payload,
        args.output_dir / "dolomite_grid.csv",
        tensor_condition="both",
        include_uncertainties=True,
        overwrite=True,
    )

    profile = thermoelasticity.read_profile_spec(args.profile_spec)
    profile_result = thermoelasticity.analyze_profiles(
        fit_result,
        (profile,),
        extrapolation_policy="fail",
    )
    profile_payload = thermoelasticity.get_result(profile_result)
    thermoelasticity.write_result(
        profile_result,
        args.output_dir / "dolomite_profile.hdf5",
    )
    thermoelasticity.write_profile_table(
        profile_payload,
        args.output_dir / "dolomite_profile.csv",
        profile_name=profile.name,
        tensor_condition="both",
        include_uncertainties=True,
        overwrite=True,
    )

    plot_jobs = [
        (
            thermoelasticity.build_fit_plots(
                fit_result,
                components=("C11",),
                options=thermoelasticity.FitPlotOptions(
                    style=thermoelasticity.PlotStyleOptions(preset="publication")
                ),
            ),
            "fit_",
        ),
        (
            thermoelasticity.build_pt_plots(
                grid_result,
                components=("C11", "C33"),
                options=thermoelasticity.PTPlotOptions(
                    style=thermoelasticity.PlotStyleOptions(preset="publication"),
                    tensor_condition="isothermal",
                    layout="facets",
                ),
            ),
            "pt_",
        ),
        (
            thermoelasticity.build_compare_plots(
                grid_result,
                components=("C11", "C33"),
                options=thermoelasticity.ComparePlotOptions(
                    style=thermoelasticity.PlotStyleOptions(preset="publication"),
                    fixed_pressure=0.0,
                    layout="facets",
                ),
            ),
            "compare_",
        ),
        (
            thermoelasticity.build_profile_plots(
                profile_result,
                profile_name=profile.name,
                components=("C11", "C33", "C44"),
                options=thermoelasticity.ProfilePlotOptions(
                    style=thermoelasticity.PlotStyleOptions(preset="publication"),
                    tensor_condition="adiabatic",
                    mode="relative",
                    layout="overlay",
                    uncertainty="none",
                    color_by="component",
                ),
            ),
            "profile_",
        ),
        (
            thermoelasticity.build_domain_plot(
                profile_result,
                profile_names=(profile.name,),
                options=thermoelasticity.DomainPlotOptions(
                    style=thermoelasticity.PlotStyleOptions(preset="publication")
                ),
            ),
            "domain_",
        ),
    ]
    for collection, prefix in plot_jobs:
        rendering.render_plots(
            collection,
            output_dir=figures,
            filename_prefix=prefix,
            preset="publication",
            dpi=180,
            close=True,
        )

    fit_payload = thermoelasticity.get_result(fit_result)
    labels = {label: index for index, label in enumerate(fit_payload.independent_labels)}
    profile_values = profile_payload.profiles[profile.name]
    c11 = profile_values.stiffness_adiabatic[:, 0, 0]
    print(f"Independent components: {', '.join(fit_payload.independent_labels)}")
    print(f"Point C11_S: {point_payload.stiffness_adiabatic[0, 0, 0, 0]:.6f} GPa")
    print(f"Point density: {point_payload.density[0, 0]:.6f} kg m^-3")
    print(f"Grid C11 index: {labels['C11']}")
    print(
        "Profile minimum C11_S change: "
        f"{100.0 * (c11.min() / c11[0] - 1.0):.6f}%"
    )
    print(f"Output directory: {args.output_dir}")


if __name__ == "__main__":
    main()
