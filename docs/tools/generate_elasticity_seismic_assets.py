"""Regenerate Elasticity and SEISMIC tutorial figures from public APIs."""

from __future__ import annotations

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
ELASTICITY_OUTPUT = ROOT / "docs" / "source" / "_static" / "tutorials" / "elasticity"
SEISMIC_OUTPUT = ROOT / "docs" / "source" / "_static" / "tutorials" / "seismic"


def main() -> None:
    """Run curated tutorial calculations and write deterministic figures."""
    source_path = str(ROOT / "src")
    if source_path not in sys.path:
        sys.path.insert(0, source_path)

    from quantas.api import common, elasticity, rendering, seismic

    ELASTICITY_OUTPUT.mkdir(parents=True, exist_ok=True)
    SEISMIC_OUTPUT.mkdir(parents=True, exist_ok=True)

    elasticity_result = elasticity.run(
        ROOT / "examples" / "elasticity" / "calcite.dat",
        options=elasticity.Options(calculate_2d=True, ntheta=181),
    )
    rendering.render_plots(
        elasticity.build_2d_plots(elasticity_result, properties=("young",)),
        output_dir=ELASTICITY_OUTPUT,
        filename_prefix="calcite_young_modulus_",
        preset="publication",
        dpi=180,
        close=True,
    )
    rendering.render_plots(
        elasticity.build_3d_plots(
            elasticity_result,
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
        ),
        output_dir=ELASTICITY_OUTPUT,
        filename_prefix="calcite_young_modulus_",
        preset="publication",
        dpi=180,
        close=True,
        axis_label_mode="crystal",
    )

    seismic_result = seismic.run(
        ROOT / "examples" / "elasticity" / "hydroxylapatite.dat",
        options=seismic.Options(
            ntheta=61,
            nphi=121,
            hemisphere="upper",
            level="enhancement",
            track_polarization_axes=True,
        ),
    )
    summary = seismic.build_summary(
        seismic_result,
        options=seismic.PlotOptions(
            include_polarizations=True,
            polarization_stride=8,
        ),
    )
    rendering.render_plots(
        common.PlotCollection(plots=[summary]),
        output_dir=SEISMIC_OUTPUT,
        filename_prefix="hydroxylapatite_",
        preset="publication",
        dpi=180,
        close=True,
        axis_label_mode="crystal",
    )
    rendering.render_plots(
        seismic.build_plots(
            seismic_result,
            seismic.PlotOptions(
                properties=("phase_v_s1",),
                projection="equal_area",
                include_extrema_markers=True,
                include_polarizations=True,
                polarization_stride=8,
            ),
        ),
        output_dir=SEISMIC_OUTPUT,
        filename_prefix="hydroxylapatite_",
        preset="publication",
        dpi=180,
        close=True,
        axis_label_mode="crystal",
    )
    rendering.render_plots(
        seismic.build_surfaces(
            seismic_result,
            seismic.SurfaceOptions(
                properties=(),
                modes=("v_p",),
                surface_types=("phase",),
                geometry="physical",
                include_polarizations=False,
            ),
        ),
        output_dir=SEISMIC_OUTPUT,
        filename_prefix="hydroxylapatite_",
        preset="publication",
        dpi=180,
        close=True,
        axis_label_mode="crystal",
    )

    # Normalize renderer-generated names to the stable documentation paths.
    renames = {
        ELASTICITY_OUTPUT / "calcite_young_modulus_2d_young.png":
            ELASTICITY_OUTPUT / "calcite_young_modulus_2d.png",
        ELASTICITY_OUTPUT / "calcite_young_modulus_3d_young.png":
            ELASTICITY_OUTPUT / "calcite_young_modulus_3d.png",
        SEISMIC_OUTPUT / "hydroxylapatite_seismic_summary_polarization.png":
            SEISMIC_OUTPUT / "hydroxylapatite_summary.png",
        SEISMIC_OUTPUT / "hydroxylapatite_seismic_2d_phase_v_s1_polarization.png":
            SEISMIC_OUTPUT / "hydroxylapatite_vs1_polarization.png",
        SEISMIC_OUTPUT / "hydroxylapatite_seismic_3d_phase_v_p_physical.png":
            SEISMIC_OUTPUT / "hydroxylapatite_vp_phase_surface.png",
    }
    for source, destination in renames.items():
        if not source.exists():
            raise RuntimeError(f"Expected tutorial figure was not generated: {source}")
        source.replace(destination)


if __name__ == "__main__":
    main()
