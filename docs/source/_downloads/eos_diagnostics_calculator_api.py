"""Inspect one fitted EOS record and calculate representative properties."""

from __future__ import annotations

import argparse
from pathlib import Path

from quantas.api import eos


def parse_args() -> argparse.Namespace:
    """Return command-line arguments for the example."""
    parser = argparse.ArgumentParser()
    parser.add_argument("archive", type=Path, help="Quantas EOS HDF5 archive")
    parser.add_argument(
        "--slot",
        default=None,
        help="Accepted DOMAIN/TARGET slot when the archive contains several results",
    )
    parser.add_argument(
        "--record-id",
        type=int,
        default=None,
        help="Explicit immutable record instead of the accepted result",
    )
    parser.add_argument(
        "--output-directory",
        type=Path,
        default=Path("eos_postfit"),
        help="Directory receiving diagnostic and calculator CSV files",
    )
    return parser.parse_args()


def main() -> None:
    """Evaluate diagnostics and a small domain-appropriate state grid."""
    args = parse_args()
    args.output_directory.mkdir(parents=True, exist_ok=True)

    diagnostics = eos.diagnose(
        args.archive,
        slot=args.slot,
        record_id=args.record_id,
    )
    diagnostic_path = eos.write_diagnostic_csv(
        diagnostics,
        args.output_directory / "diagnostics.csv",
        overwrite=True,
    )

    domain = eos.record_domain(
        args.archive,
        slot=args.slot,
        record_id=args.record_id,
    ).value
    calculation_kwargs: dict[str, list[float]]
    if domain == "pv":
        calculation_kwargs = {"pressure": [0.0, 5.0, 10.0]}
    elif domain == "vt":
        calculation_kwargs = {"temperature": [300.0, 600.0, 900.0]}
    elif domain == "pvt":
        calculation_kwargs = {
            "pressure": [0.0, 5.0, 10.0],
            "temperature": [300.0, 600.0, 900.0],
        }
    else:  # pragma: no cover - guarded by the public API
        raise RuntimeError(f"unsupported EOS domain: {domain}")
    calculation = eos.calculate(
        args.archive,
        slot=args.slot,
        record_id=args.record_id,
        **calculation_kwargs,
    )

    calculation_path = eos.write_calculation_csv(
        calculation,
        args.output_directory / "calculated_properties.csv",
        overwrite=True,
    )

    print("Quantas EOS post-fit analysis")
    print("=============================")
    print(f"Record                  : {calculation.record_id}")
    print(f"Slot                    : {calculation.slot.key}")
    print(f"Diagnostic observations : {diagnostics.nrows}")
    print(f"Calculated states       : {calculation.nrows}")
    print(f"Diagnostics CSV         : {diagnostic_path}")
    print(f"Properties CSV          : {calculation_path}")
    if calculation.warnings:
        print("Warnings:")
        for warning in calculation.warnings:
            print(f"- {warning}")


if __name__ == "__main__":
    main()
