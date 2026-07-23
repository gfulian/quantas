"""Generate and validate quasi-static thermoelastic input through the API.

Examples
--------
Generate the elastic YAML from a list file and pair it with an existing QHA
HDF5 record::

    python qsa_input_api.py --soec-list dol_soec.txt \
        --output dol_thermoelastic.yaml --qha dol_qha.hdf5
"""

from __future__ import annotations

import argparse
from pathlib import Path

from quantas.api import qha, thermoelasticity


def main() -> None:
    """Generate elastic input and validate its coupling to a QHA source."""
    parser = argparse.ArgumentParser()
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--soec-list", type=Path)
    source.add_argument("--soec", type=Path, nargs="+")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--qha", type=Path, required=True)
    parser.add_argument("--jobname", default="Quasi-static thermoelastic example")
    args = parser.parse_args()

    sources = args.soec_list if args.soec_list is not None else args.soec
    output = thermoelasticity.create_input(
        sources,
        args.output,
        is_list=args.soec_list is not None,
        jobname=args.jobname,
    )

    qha_options = None
    if args.qha.suffix.lower() not in {".h5", ".hdf5"}:
        qha_options = qha.Options(
            temperature_min=300.0,
            temperature_max=1000.0,
            temperature_step=100.0,
            pressure_min=0.0,
            pressure_max=5.0,
            pressure_step=1.0,
        )

    context = thermoelasticity.prepare_context(
        output,
        args.qha,
        qha_options=qha_options,
    )
    series = context.input_data.elastic_series
    print(f"Thermoelastic input: {output}")
    print(f"Elastic points: {series.npoints}")
    print(
        "Elastic volume interval: "
        f"{series.volume_bounds[0]:.6f} -- {series.volume_bounds[1]:.6f} A^3"
    )
    print(
        "QHA volume interval: "
        f"{context.metadata['qha_volume_min']:.6f} -- "
        f"{context.metadata['qha_volume_max']:.6f} A^3"
    )
    print(
        "QHA points outside elastic coverage: "
        f"{context.metadata['extrapolated_points']} / "
        f"{context.metadata['total_points']}"
    )
    print(f"Complete cold-QSA inputs: {context.has_complete_quasistatic_inputs}")
    print(
        "Complete anisotropic adiabatic inputs: "
        f"{context.has_complete_adiabatic_inputs}"
    )
    if context.missing_qha_fields:
        print("Missing QHA fields: " + ", ".join(context.missing_qha_fields))


if __name__ == "__main__":
    main()
