"""Run the supported Quantas cross-module workflow through the public API.

The example couples a dolomite QHA result to the quasi-static thermoelastic
workflow, reconstructs one adiabatic state, and sends the same state to both
Elasticity and SEISMIC without importing private implementation modules.

Run from the project root::

    python examples/interoperability/workflow_api.py \
        --output-dir interoperability_example
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from quantas.api import (
    elasticity,
    interop,
    qha,
    rendering,
    seismic,
    thermoelasticity,
)


def parser() -> argparse.ArgumentParser:
    """Return the standalone-example command-line parser."""
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument(
        "--qha-input",
        type=Path,
        default=Path("examples/qha/crystal-phonons/dol_pbe0.yaml"),
    )
    result.add_argument(
        "--thermoelastic-input",
        type=Path,
        default=Path("examples/thermoelasticity/dol_pbe0_thermoelastic.yaml"),
    )
    result.add_argument("--pressure", type=float, default=5.0)
    result.add_argument("--temperature", type=float, default=800.0)
    result.add_argument(
        "--output-dir",
        type=Path,
        default=Path("interoperability_example"),
    )
    return result


def main() -> None:
    """Execute QHA -> Thermoelasticity -> Elasticity/SEISMIC."""
    args = parser().parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    qha_result = qha.run(
        args.qha_input,
        options=qha.Options(
            temperature_min=300.0,
            temperature_max=1800.0,
            temperature_step=100.0,
            pressure_min=0.0,
            pressure_max=10.0,
            pressure_step=2.0,
            scheme="freq",
            minimization="poly",
            energy_degree=3,
            frequency_degree=3,
            thermal_expansion_method="mixed_derivative",
        ),
    )
    qha_path = qha.write_result(
        qha_result,
        args.output_dir / "dolomite_qha.hdf5",
    )

    context = interop.qha_to_thermoelastic_context(
        args.thermoelastic_input,
        qha_result,
    )
    fit_result = thermoelasticity.run_context(
        context,
        options=thermoelasticity.Options(
            reference_eos="BM3",
            finite_strain_order=3,
            adiabatic_mode="auto",
            quality_policy="warn",
            report_level="standard",
        ),
    )
    fit_path = thermoelasticity.write_result(
        fit_result,
        args.output_dir / "dolomite_fit.hdf5",
    )

    seismic_input = interop.thermoelastic_to_seismic(
        fit_result,
        pressure=args.pressure,
        temperature=args.temperature,
        tensor_condition="adiabatic",
        extrapolation_policy="fail",
    )

    elasticity_input = elasticity.Input(
        jobname=seismic_input.jobname,
        stiffness=np.asarray(seismic_input.stiffness, dtype=np.float64),
        source=seismic_input.source,
    )
    elasticity_result = elasticity.run(elasticity_input)
    elasticity_path = elasticity.write_result(
        elasticity_result,
        args.output_dir / "dolomite_state_elasticity.hdf5",
    )

    seismic_result = seismic.run(
        seismic_input,
        options=seismic.Options(
            ntheta=31,
            nphi=61,
            hemisphere="upper",
            level="phase",
            track_polarization_axes=True,
        ),
    )
    seismic_path = seismic.write_result(
        seismic_result,
        args.output_dir / "dolomite_state_seismic.hdf5",
    )

    state_path = args.output_dir / "dolomite_5GPa_800K.dat"
    point_result = thermoelasticity.analyze_grid(
        fit_result,
        pressure=np.asarray([args.pressure], dtype=np.float64),
        temperature=np.asarray([args.temperature], dtype=np.float64),
        extrapolation_policy="fail",
    )
    thermoelasticity.write_state_input(
        thermoelasticity.get_result(point_result),
        state_path,
        tensor_condition="adiabatic",
        overwrite=True,
    )

    elastic = elasticity.get_result(elasticity_result)
    seismic_payload = seismic.get_result(seismic_result)
    vp = seismic_payload.field.phase.for_mode("v_p").phase_speeds

    print(f"QHA archive: {qha_path}")
    print(f"Thermoelastic fit archive: {fit_path}")
    print(f"Shared state input: {state_path}")
    print(f"Elasticity archive: {elasticity_path}")
    print(f"SEISMIC archive: {seismic_path}")
    print(f"Context extrapolated states: {context.metadata['extrapolated_points']}")
    print(f"State density / kg m^-3: {seismic_input.density:.6f}")
    print(f"State C11_S / GPa: {seismic_input.stiffness[0, 0]:.6f}")
    print(f"Hill bulk modulus / GPa: {elastic.averages.hill.bulk_modulus:.6f}")
    print(f"Sampled V_P range / km s^-1: {np.nanmin(vp):.6f} {np.nanmax(vp):.6f}")

    report = rendering.render_tables(elasticity.build_report(elasticity_result))
    (args.output_dir / "dolomite_state_elasticity.log").write_text(
        report,
        encoding="utf-8",
        newline="\n",
    )


if __name__ == "__main__":
    main()
