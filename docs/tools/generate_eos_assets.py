"""Regenerate EOS tutorial archives, tables, reports, and figures."""

from __future__ import annotations

import csv
from pathlib import Path
import shutil
import sys
import tempfile

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
STATIC_OUTPUT = ROOT / "docs" / "source" / "_static" / "tutorials" / "eos"
DOWNLOAD_OUTPUT = ROOT / "docs" / "source" / "_downloads" / "tutorials" / "eos"
EXAMPLE_ROOT = ROOT / "examples" / "eos"

CASES = {
    "pv": (
        EXAMPLE_ROOT / "PV_quartz.dat",
        EXAMPLE_ROOT / "specs" / "quartz_pv_tutorial.spec",
        "pv/volume",
    ),
    "vt": (
        EXAMPLE_ROOT / "VT_rutile.dat",
        EXAMPLE_ROOT / "specs" / "rutile_vt_tutorial.spec",
        "vt/volume",
    ),
    "pvt": (
        EXAMPLE_ROOT / "PVT_NaF.dat",
        EXAMPLE_ROOT / "specs" / "naf_pvt_tutorial.spec",
        "pvt/volume",
    ),
}


def _write_rows(
    path: Path,
    fieldnames: list[str],
    rows: list[dict[str, object]],
) -> None:
    """Write one deterministic CSV table."""
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _parameter_error(record: object, name: str) -> float | None:
    """Return one named parameter standard error from an archive record."""
    fit = record.result.fit
    if fit.errors is None or name not in fit.parameter_names:
        return None
    return float(fit.errors[fit.parameter_names.index(name)])


def _run_case(kind: str, work: Path) -> Path:
    """Run one declarative EOS batch and return its archive path."""
    from quantas.api import eos, rendering

    data_path, spec_path, _ = CASES[kind]
    dataset = eos.read_input(data_path)
    resolved = eos.resolve_spec(eos.read_spec(spec_path), dataset)
    for job in resolved.plan.jobs:
        eos.validate_request(dataset, job.request)

    archive = work / f"{kind}.hdf5"
    result = eos.run_batch(
        dataset,
        resolved.plan,
        archive,
        overwrite=True,
        creator="docs/tools/generate_eos_assets.py",
    )
    report = rendering.render_tables(
        (
            *eos.build_batch_preamble(
                dataset,
                resolved.plan,
                archive,
                resolved.report_options,
            ),
            *eos.build_batch_report(result, resolved.report_options),
        )
    )
    (work / f"{kind}.log").write_text(report, encoding="utf-8")
    return archive


def _write_comparisons(archives: dict[str, Path]) -> None:
    """Write solver, model, and coupling comparison tables."""
    from quantas.api import eos

    with eos.open_archive(archives["pv"]) as archive:
        rows = []
        for record_id in archive.record_ids:
            record = archive.record(record_id)
            fit = record.result.fit
            diagnostics = fit.diagnostics
            values = record.result.parameter_values
            rows.append(
                {
                    "record_id": record_id,
                    "job": record.request.request_id,
                    "model": record.request.model.tag,
                    "solver": fit.method.value if fit.method is not None else "",
                    "V0": f"{values['V0']:.9f}",
                    "sigma_V0": f"{_parameter_error(record, 'V0'):.9f}",
                    "K0_GPa": f"{values['K0']:.9f}",
                    "sigma_K0_GPa": f"{_parameter_error(record, 'K0'):.9f}",
                    "KP": f"{values['KP']:.9f}",
                    "sigma_KP": f"{_parameter_error(record, 'KP'):.9f}",
                    "RMSE_GPa": f"{fit.rmse:.9f}",
                    "reduced_chi_square": ""
                    if diagnostics is None or diagnostics.reduced_chi_square is None
                    else f"{diagnostics.reduced_chi_square:.9f}",
                }
            )
        _write_rows(
            DOWNLOAD_OUTPUT / "quartz_pv_comparison.csv",
            list(rows[0]),
            rows,
        )

    with eos.open_archive(archives["vt"]) as archive:
        rows = []
        for record_id in archive.record_ids:
            record = archive.record(record_id)
            if record.request.target != "volume":
                continue
            fit = record.result.fit
            diagnostics = fit.diagnostics
            values = record.result.parameter_values
            calculated = eos.calculate(
                archives["vt"],
                record_id=record_id,
                temperature=[300.0, 1000.0],
                propagate_uncertainty=False,
            )
            alpha = np.asarray(calculated.columns["expansion_coefficient"])
            rows.append(
                {
                    "record_id": record_id,
                    "job": record.request.request_id,
                    "model": record.request.model.tag,
                    "solver": fit.method.value if fit.method is not None else "",
                    "V0": f"{values['V0']:.10f}",
                    "alpha_300_K": f"{alpha[0]:.10e}",
                    "alpha_1000_K": f"{alpha[1]:.10e}",
                    "RMSE": f"{fit.rmse:.10e}",
                    "reduced_chi_square": ""
                    if diagnostics is None or diagnostics.reduced_chi_square is None
                    else f"{diagnostics.reduced_chi_square:.9f}",
                }
            )
        _write_rows(
            DOWNLOAD_OUTPUT / "rutile_vt_comparison.csv",
            list(rows[0]),
            rows,
        )

    with eos.open_archive(archives["pvt"]) as archive:
        rows = []
        for record_id in archive.record_ids:
            record = archive.record(record_id)
            fit = record.result.fit
            diagnostics = fit.diagnostics
            values = record.result.parameter_values
            calculated = eos.calculate(
                archives["pvt"],
                record_id=record_id,
                pressure=[1.5],
                temperature=[230.0],
                propagate_uncertainty=False,
            )
            rows.append(
                {
                    "record_id": record_id,
                    "job": record.request.request_id,
                    "model": record.request.model.tag,
                    "K0_GPa": f"{values['K0']:.9f}",
                    "KP": f"{values['KP']:.9f}",
                    "V0": f"{values['V0']:.9f}",
                    "V_1.5GPa_230K": f"{float(calculated.columns['volume'][0]):.9f}",
                    "K_1.5GPa_230K_GPa": (
                        f"{float(calculated.columns['bulk_modulus'][0]):.9f}"
                    ),
                    "alpha_1.5GPa_230K": (
                        f"{float(calculated.columns['expansion_coefficient'][0]):.10e}"
                    ),
                    "RMSE_GPa": f"{fit.rmse:.9f}",
                    "reduced_chi_square": ""
                    if diagnostics is None or diagnostics.reduced_chi_square is None
                    else f"{diagnostics.reduced_chi_square:.9f}",
                }
            )
        _write_rows(
            DOWNLOAD_OUTPUT / "naf_pvt_coupling_comparison.csv",
            list(rows[0]),
            rows,
        )


def _render_plots(archives: dict[str, Path]) -> None:
    """Render the compact set of figures embedded by the tutorial."""
    from quantas.api import eos, rendering

    specifications = {
        "pv": (
            ("fit", "residuals", "normalized-pressure"),
            eos.PlotOptions(
                show_uncertainties=True,
                point_size=5.0,
                curve_width=1.8,
            ),
        ),
        "vt": (
            ("fit", "residuals"),
            eos.PlotOptions(
                show_uncertainties=True,
                point_size=5.0,
                curve_width=1.8,
            ),
        ),
        "pvt": (
            ("coverage", "isotherms", "isobars", "residuals"),
            eos.PlotOptions(
                show_uncertainties=False,
                point_size=5.0,
                curve_width=1.8,
                curve_points=60,
                isotherms=(140.0, 230.0, 350.0),
                isobars=(0.5, 1.5, 2.5),
            ),
        ),
    }
    for kind, (plot_types, options) in specifications.items():
        _, _, slot = CASES[kind]
        collection = eos.build_plots(
            archives[kind],
            plot_types,
            slot=slot,
            options=options,
        )
        rendered = rendering.render_plots(
            collection,
            output_dir=STATIC_OUTPUT,
            filename_prefix=f"{kind}_",
            preset="publication",
            dpi=180,
            close=True,
        )
        for warning in rendered.warnings:
            print(f"Warning while rendering {kind}: {warning}")

    generated_names = {
        "pv_record_000003_pv_volume_fit.png": "quartz_pv_fit.png",
        "pv_record_000003_pv_volume_residuals_vs_pressure.png": (
            "quartz_pv_residuals.png"
        ),
        "pv_record_000003_pv_volume_normalized_pressure.png": (
            "quartz_pv_normalized_pressure.png"
        ),
        "vt_record_000003_vt_volume_fit.png": "rutile_vt_fit.png",
        "vt_record_000003_vt_volume_residuals_vs_temperature.png": (
            "rutile_vt_residuals.png"
        ),
        "pvt_record_000004_pvt_volume_coverage.png": "naf_pvt_coverage.png",
        "pvt_record_000004_pvt_volume_isotherms.png": "naf_pvt_isotherms.png",
        "pvt_record_000004_pvt_volume_isobars.png": "naf_pvt_isobars.png",
        "pvt_record_000004_pvt_volume_residuals_vs_pressure.png": (
            "naf_pvt_residuals_pressure.png"
        ),
        "pvt_record_000004_pvt_volume_residuals_vs_temperature.png": (
            "naf_pvt_residuals_temperature.png"
        ),
    }
    renames = {
        STATIC_OUTPUT / source: STATIC_OUTPUT / destination
        for source, destination in generated_names.items()
    }
    for source, destination in renames.items():
        if not source.is_file():
            raise RuntimeError(
                f"Expected EOS tutorial figure was not generated: {source}"
            )
        source.replace(destination)


def main() -> None:
    """Regenerate all EOS tutorial assets from public APIs."""
    source = str(ROOT / "src")
    if source not in sys.path:
        sys.path.insert(0, source)
    STATIC_OUTPUT.mkdir(parents=True, exist_ok=True)
    DOWNLOAD_OUTPUT.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="quantas-eos-docs-") as tmp:
        work = Path(tmp)
        archives = {kind: _run_case(kind, work) for kind in CASES}
        _write_comparisons(archives)
        _render_plots(archives)
        for kind in CASES:
            shutil.copyfile(
                work / f"{kind}.log",
                DOWNLOAD_OUTPUT / f"{kind}_report.txt",
            )

    for source_path in (
        EXAMPLE_ROOT / "specs" / "quartz_pv_tutorial.spec",
        EXAMPLE_ROOT / "specs" / "rutile_vt_tutorial.spec",
        EXAMPLE_ROOT / "specs" / "naf_pvt_tutorial.spec",
        EXAMPLE_ROOT / "scripts" / "run_spec_api.py",
        EXAMPLE_ROOT / "scripts" / "pv_fit_api.py",
        EXAMPLE_ROOT / "scripts" / "vt_fit_api.py",
        EXAMPLE_ROOT / "scripts" / "pvt_fit_api.py",
    ):
        shutil.copyfile(source_path, DOWNLOAD_OUTPUT / source_path.name)


if __name__ == "__main__":
    main()
