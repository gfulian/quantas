"""Run a quasi-harmonic calculation through the public Quantas Python API."""

from __future__ import annotations

import argparse
from pathlib import Path

from quantas.api import qha, rendering


def main() -> None:
    """Read input, run QHA, write HDF5, and render a neutral report."""
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path, default=Path("qha_result.hdf5"))
    args = parser.parse_args()

    input_data = qha.read_input(args.input)
    options = qha.Options()
    result = qha.run(input_data, options=options)
    output = qha.write_result(result, args.output)

    print(rendering.render_tables(qha.build_report(result)))
    print(f"Written: {output}")


if __name__ == "__main__":
    main()
