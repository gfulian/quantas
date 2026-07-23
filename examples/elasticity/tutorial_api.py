"""Reproduce the Calcite Elasticity tutorial through the public Quantas API."""

from __future__ import annotations

import argparse
from pathlib import Path

from quantas.api import elasticity, rendering


def main() -> None:
    """Run the calculation, persist results, and render selected figures."""
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("elasticity_tutorial"),
    )
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    options = elasticity.Options(calculate_2d=True, ntheta=181)
    result = elasticity.run(args.input, options=options)
    payload = elasticity.get_result(result)

    report = rendering.render_tables(elasticity.build_report(result))
    report_path = args.output_dir / "calcite_elasticity_api.log"
    report_path.write_text(report, encoding="utf-8", newline="\n")
    hdf5_path = elasticity.write_result(
        result,
        args.output_dir / "calcite_elasticity_api.hdf5",
        report_text=report,
    )

    plots_2d = elasticity.build_2d_plots(
        result,
        properties=("young", "compressibility"),
    )
    rendering.render_plots(
        plots_2d,
        output_dir=args.output_dir,
        filename_prefix="calcite_",
        preset="publication",
        dpi=180,
        close=True,
    )

    plots_3d = elasticity.build_3d_plots(
        result,
        elasticity.SurfaceOptions(
            ntheta=61,
            nphi=121,
            properties=("young",),
        ),
        properties=("young",),
        geometry="physical",
        color_mode="property",
        show_mesh=True,
        mesh_color="black",
        mesh_line_width=0.35,
    )
    rendering.render_plots(
        plots_3d,
        output_dir=args.output_dir,
        filename_prefix="calcite_",
        preset="publication",
        dpi=180,
        close=True,
        axis_label_mode="crystal",
    )

    hill = payload.averages.hill
    young = payload.variations["young_modulus"]
    print(f"Crystal system: {payload.crystal_system}")
    print(f"Mechanically stable: {payload.stability.is_stable}")
    print(f"Hill bulk modulus / GPa: {hill.bulk_modulus:.6f}")
    print(f"Hill shear modulus / GPa: {hill.shear_modulus:.6f}")
    print(f"Young minimum / GPa: {young.minimum:.6f}")
    print(f"Young maximum / GPa: {young.maximum:.6f}")
    print(f"Written: {hdf5_path}")
    print(f"Written: {report_path}")


if __name__ == "__main__":
    main()
