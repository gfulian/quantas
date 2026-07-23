"""Reproduce the HA MgO tutorial through the public Quantas API."""

from __future__ import annotations

import argparse
from pathlib import Path

from quantas.api import ha, rendering


def main() -> None:
    """Run HA, write HDF5/report files, and render selected figures."""
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path("ha_tutorial"))
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    options = ha.Options(
        temperature_min=0.0,
        temperature_max=1000.0,
        temperature_step=100.0,
        energy_unit="Ha",
        volume_unit="A",
        frequency_unit="cm^-1",
        temperature_unit="K",
    )

    result = ha.run(args.input, options=options)
    payload = ha.get_result(result)

    report = rendering.render_tables(ha.build_report(result))
    report_path = args.output_dir / "mgo_ha_api.log"
    report_path.write_text(report, encoding="utf-8")

    hdf5_path = ha.write_result(
        result,
        args.output_dir / "mgo_ha_api.hdf5",
        report_text=report,
    )

    plot_requests = (
        ("Cv", "J/mol", "mgo_ha_cv_"),
        ("Fvib", "kJ/mol", "mgo_ha_fvib_"),
    )
    for property_name, unit, prefix in plot_requests:
        collection = ha.build_plots(
            result,
            properties=property_name,
            unit=unit,
        )
        rendering.render_plots(
            collection,
            output_dir=args.output_dir,
            filename_prefix=prefix,
            preset="publication",
            dpi=180,
            close=True,
        )

    print(f"Temperatures: {payload.temperature.size}")
    print(f"Volumes: {payload.volume.size}")
    print(f"Cv shape: {payload.isochoric_heat_capacity.shape}")
    print(f"Written: {hdf5_path}")
    print(f"Written: {report_path}")


if __name__ == "__main__":
    main()
