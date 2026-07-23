"""Generate static three-dimensional elasticity surfaces through the API."""

from __future__ import annotations

import argparse
from pathlib import Path

from quantas.api import elasticity, rendering


def main() -> None:
    """Read an HDF5 tensor and render transient 3D plot specifications."""
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    parser.add_argument("--ntheta", type=int, default=61)
    parser.add_argument("--nphi", type=int, default=121)
    args = parser.parse_args()

    result = elasticity.read_result(args.input)
    plots = elasticity.build_3d_plots(
        result,
        elasticity.SurfaceOptions(ntheta=args.ntheta, nphi=args.nphi),
    )
    rendered = rendering.render_plots(
        plots,
        output_dir=args.output_dir,
        close=True,
    )
    print(f"Generated {len(rendered.artifacts)} static figures")


if __name__ == "__main__":
    main()
