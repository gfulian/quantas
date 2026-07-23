"""Run the first Quantas elasticity calculation through the public API."""

from __future__ import annotations

from pathlib import Path

from quantas.api import elasticity, rendering


INPUT = Path("calcite.dat")
RESULT = Path("calcite_elasticity_api.hdf5")
REPORT = Path("calcite_elasticity_api.log")


def main() -> None:
    """Calculate, report, and persist the calcite elasticity example."""
    result = elasticity.run(INPUT)
    payload = elasticity.get_result(result)

    report_text = rendering.render_tables(elasticity.build_report(result))
    REPORT.write_text(report_text, encoding="utf-8")
    elasticity.write_result(result, RESULT, report_text=report_text)

    young = payload.variations["young_modulus"]
    print(f"Crystal system: {payload.crystal_system}")
    print(f"Mechanically stable: {payload.stability.is_stable}")
    print(f"Hill bulk modulus / GPa: {payload.averages.hill.bulk_modulus:.6f}")
    print(f"Hill shear modulus / GPa: {payload.averages.hill.shear_modulus:.6f}")
    print(f"Young minimum / GPa: {young.minimum:.6f}")
    print(f"Young maximum / GPa: {young.maximum:.6f}")


if __name__ == "__main__":
    main()
