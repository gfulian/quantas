"""Generate standard EOS plots from an existing Quantas HDF5 archive."""

from __future__ import annotations

import argparse
from pathlib import Path

from quantas.api import eos, rendering


def parse_args() -> argparse.Namespace:
    """Return command-line arguments for the tutorial script."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path, help="Quantas EOS HDF5 archive")
    parser.add_argument(
        "--slot",
        help="Accepted EOS slot, for example pv/volume or pvt/volume",
    )
    parser.add_argument(
        "--record-id",
        type=int,
        help="Explicit immutable record identifier",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("eos_plots"),
        help="Directory receiving plot files",
    )
    parser.add_argument("--dpi", type=int, default=150)
    parser.add_argument("--point-size", type=float, default=5.0)
    parser.add_argument("--curve-width", type=float, default=1.8)
    parser.add_argument(
        "--no-uncertainties",
        action="store_true",
        help="Suppress input one-sigma error bars",
    )
    return parser.parse_args()


def main() -> None:
    """Build neutral plot specifications and render them with Matplotlib."""
    args = parse_args()
    inventory = eos.describe_plots(
        args.archive,
        slot=args.slot,
        record_id=args.record_id,
    )
    for warning in inventory.warnings:
        print(f"Warning: {warning}")
    if inventory.selected_plots is None:
        print("No unique EOS fit record is selected. Available slots:")
        for slot in inventory.slots:
            print(
                f"  {slot.key}: status={slot.status.value}; "
                f"accepted={slot.accepted_record_id}; "
                f"plottable={slot.plottable_record_ids}"
            )
        return
    available = tuple(
        item.key for item in inventory.selected_plots.representations
    )
    print("Available plot types:", ", ".join(available))

    collection = eos.build_plots(
        args.archive,
        available,
        slot=args.slot,
        record_id=args.record_id,
        options=eos.PlotOptions(
            show_uncertainties=not args.no_uncertainties,
            point_size=args.point_size,
            curve_width=args.curve_width,
        )
    )
    result = rendering.render_plots(
        collection,
        output_dir=args.output,
        image_format="png",
        dpi=args.dpi,
        close=True,
    )
    for warning in result.warnings:
        print(f"Warning: {warning}")
    for artifact in result.artifacts:
        if artifact.path is not None:
            print(f"Written {artifact.path}")


if __name__ == "__main__":
    main()
