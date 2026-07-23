"""Reproduce the QHA MgO tutorial through the public Quantas API."""

from __future__ import annotations

import argparse
from pathlib import Path

from quantas.api import qha, rendering


def main() -> None:
    """Inspect and run QHA, then write reports, HDF5, and selected figures."""
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path("qha_tutorial"))
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    options = qha.Options(
        temperature_min=0.0,
        temperature_max=1000.0,
        temperature_step=100.0,
        pressure_min=0.0,
        pressure_max=10.0,
        pressure_step=2.0,
        scheme="freq",
        minimization="poly",
        eos="BM3",
        energy_degree=3,
        free_energy_degree=3,
        frequency_degree=3,
        structural_degree=3,
        thermal_expansion_method="mixed_derivative",
        calculate_gruneisen=True,
        calculate_mode_gruneisen=True,
    )

    preview = qha.inspect(
        args.input,
        options=options,
        include_polynomial=True,
        include_eos=True,
        polynomial_degree=3,
        eos="BM3",
    )
    print(
        "Static pressure coverage from BM3: "
        f"{preview.eos.pressure_min:.2f} to {preview.eos.pressure_max:.2f} "
        f"{preview.pressure_unit}"
    )

    result = qha.run(args.input, options=options)
    payload = qha.get_result(result)

    report = rendering.render_tables(qha.build_report(result))
    report_path = args.output_dir / "mgo_qha_api.log"
    report_path.write_text(report, encoding="utf-8")

    hdf5_path = qha.write_result(
        result,
        args.output_dir / "mgo_qha_api.hdf5",
        report_text=report,
    )

    volume_plots = qha.build_plots(
        result,
        properties=["VT"],
        options=qha.PlotOptions(),
    )
    rendering.render_plots(
        volume_plots,
        output_dir=args.output_dir,
        filename_prefix="mgo_qha_",
        preset="publication",
        dpi=180,
        close=True,
    )

    heat_capacity_plots = qha.build_plots(
        result,
        properties=["heat_capacities"],
        options=qha.PlotOptions(energy_unit="J/mol"),
    )
    rendering.render_plots(
        heat_capacity_plots,
        output_dir=args.output_dir,
        filename_prefix="mgo_qha_",
        preset="publication",
        dpi=180,
        close=True,
    )

    expansion_plots = qha.build_plots(
        result,
        properties=["alphaV"],
        options=qha.PlotOptions(
            include_contours=True,
            cmap="viridis",
            levels=12,
        ),
    )
    rendering.render_plots(
        expansion_plots,
        output_dir=args.output_dir,
        filename_prefix="mgo_qha_",
        preset="publication",
        dpi=180,
        close=True,
    )

    workflow = payload.metadata["workflow"]
    print(f"Grid shape: {payload.temperature.size} x {payload.pressure.size}")
    print(f"Successful points: {workflow['successful_points']}")
    print(f"Failed points: {workflow['failed_points']}")
    print(f"Written: {hdf5_path}")
    print(f"Written: {report_path}")


if __name__ == "__main__":
    main()
