"""Run an elasticity calculation through the public Quantas Python API."""

from __future__ import annotations

import argparse
from pathlib import Path

from quantas.api import elasticity, rendering


def main() -> None:
    """Read input, run elasticity, write HDF5, and render a neutral report."""
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path, default=Path("elasticity_result.hdf5"))
    parser.add_argument("--2d", dest="calculate_2d", action="store_true")
    args = parser.parse_args()

    input_data = elasticity.read_input(args.input)
    options = elasticity.Options(calculate_2d=args.calculate_2d)
    result = elasticity.run(input_data, options=options)
    output = elasticity.write_result(result, args.output)

    print(rendering.render_tables(elasticity.build_report(result)))
    print(f"Written: {output}")


if __name__ == "__main__":
    main()
