"""Run one QUANTAS EOS SPEC 1 batch through the public Python API."""

from __future__ import annotations

import argparse
from pathlib import Path

from quantas.api import eos, rendering


def parse_args() -> argparse.Namespace:
    """Return command-line arguments for the reusable batch example."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data", type=Path, help="EOS data table")
    parser.add_argument("spec", type=Path, help="QUANTAS EOS SPEC 1 file")
    parser.add_argument("archive", type=Path, help="Output EOS HDF5 archive")
    parser.add_argument(
        "--report",
        type=Path,
        default=None,
        help="Optional deterministic plain-text report",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Resolve and validate every request without fitting",
    )
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def main() -> None:
    """Resolve the declarative plan, validate it, and optionally execute it."""
    args = parse_args()
    dataset = eos.read_input(args.data)
    document = eos.read_spec(args.spec)
    resolved = eos.resolve_spec(document, dataset)

    for job in resolved.plan.jobs:
        eos.validate_request(dataset, job.request)

    preamble = eos.build_batch_preamble(
        dataset,
        resolved.plan,
        args.archive,
        resolved.report_options,
    )
    print(rendering.render_tables(preamble))

    if args.dry_run:
        print("Dry run completed: every expanded fit request is valid.")
        return

    result = eos.run_batch(
        dataset,
        resolved.plan,
        args.archive,
        overwrite=args.force,
        creator="examples/eos/scripts/run_spec_api.py",
    )
    report_tables = (
        *preamble,
        *eos.build_batch_report(result, resolved.report_options),
    )
    report_text = rendering.render_tables(report_tables)
    print(report_text)
    if args.report is not None:
        if args.report.exists() and not args.force:
            raise FileExistsError(args.report)
        args.report.write_text(report_text, encoding="utf-8")


if __name__ == "__main__":
    main()
