"""Run a harmonic calculation through the public Quantas Python API."""

from __future__ import annotations

import argparse
from pathlib import Path

from quantas.api import ha, rendering


def main() -> None:
    """Read input, run HA, write HDF5, and render a neutral report."""
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path, default=Path("ha_result.hdf5"))
    args = parser.parse_args()

    input_data = ha.read_input(args.input)
    options = ha.Options()
    result = ha.run(input_data, options=options)
    output = ha.write_result(result, args.output)

    print(rendering.render_tables(ha.build_report(result)))
    print(f"Written: {output}")


if __name__ == "__main__":
    main()
