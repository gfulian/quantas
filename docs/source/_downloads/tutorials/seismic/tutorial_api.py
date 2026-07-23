"""Reproduce the Hydroxylapatite SEISMIC tutorial with the public API."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from quantas.api import common, rendering, seismic


def main() -> None:
    """Run SEISMIC, persist fields, export CSV, and render selected figures."""
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("seismic_tutorial"),
    )
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    options = seismic.Options(
        ntheta=61,
        nphi=121,
        hemisphere="upper",
        level="enhancement",
        track_polarization_axes=True,
    )
    result = seismic.run(args.input, options=options)
    payload = seismic.get_result(result)

    report = rendering.render_tables(
        seismic.build_report(result, level="extended")
    )
    report_path = args.output_dir / "hydroxylapatite_seismic_api.log"
    report_path.write_text(report, encoding="utf-8", newline="\n")
    hdf5_path = seismic.write_result(
        result,
        args.output_dir / "hydroxylapatite_seismic_api.hdf5",
        report_text=report,
    )
    csv_path = seismic.write_csv(
        result,
        args.output_dir / "hydroxylapatite_seismic_api.csv",
    )

    summary_options = seismic.PlotOptions(
        include_polarizations=True,
        polarization_stride=8,
    )
    summary = seismic.build_summary(result, options=summary_options)
    rendering.render_plots(
        common.PlotCollection(plots=[summary]),
        output_dir=args.output_dir,
        filename_prefix="hydroxylapatite_",
        preset="publication",
        dpi=180,
        close=True,
        axis_label_mode="crystal",
    )

    vs1_map = seismic.build_plots(
        result,
        seismic.PlotOptions(
            properties=("phase_v_s1",),
            projection="equal_area",
            include_extrema_markers=True,
            include_polarizations=True,
            polarization_stride=8,
        ),
    )
    rendering.render_plots(
        vs1_map,
        output_dir=args.output_dir,
        filename_prefix="hydroxylapatite_",
        preset="publication",
        dpi=180,
        close=True,
        axis_label_mode="crystal",
    )

    vp_surface = seismic.build_surfaces(
        result,
        seismic.SurfaceOptions(
            properties=(),
            modes=("v_p",),
            surface_types=("phase",),
            geometry="physical",
            include_polarizations=False,
        ),
    )
    rendering.render_plots(
        vp_surface,
        output_dir=args.output_dir,
        filename_prefix="hydroxylapatite_",
        preset="publication",
        dpi=180,
        close=True,
        axis_label_mode="crystal",
    )

    phase = payload.field.phase
    vp = phase.for_mode("v_p")
    vs1 = phase.for_mode("v_s1")
    vs2 = phase.for_mode("v_s2")
    print(f"Sampled directions: {phase.n_points}")
    print(f"Isotropic V_S / km s^-1: {payload.isotropic_velocities.shear:.8f}")
    print(
        "Isotropic V_P / km s^-1: "
        f"{payload.isotropic_velocities.compressional:.8f}"
    )
    print(f"V_P range / km s^-1: {np.nanmin(vp.phase_speeds):.8f} "
          f"{np.nanmax(vp.phase_speeds):.8f}")
    print(f"V_S1 range / km s^-1: {np.nanmin(vs1.phase_speeds):.8f} "
          f"{np.nanmax(vs1.phase_speeds):.8f}")
    print(f"V_S2 range / km s^-1: {np.nanmin(vs2.phase_speeds):.8f} "
          f"{np.nanmax(vs2.phase_speeds):.8f}")
    print(f"Written: {hdf5_path}")
    print(f"Written: {report_path}")
    print(f"Written: {csv_path}")


if __name__ == "__main__":
    main()
