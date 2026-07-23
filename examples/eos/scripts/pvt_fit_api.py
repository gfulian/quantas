"""Fit the NaF P-V-T example with a full Mie-Gruneisen-Debye model."""

from __future__ import annotations

from pathlib import Path

from quantas.api import eos

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "PVT_NaF.dat"


def main() -> None:
    """Run the public volume-only MGD fit using molar-volume normalization."""
    dataset = eos.read_input(DATA)
    normalization = eos.MGDNormalization.molar_formula_unit(formula="NaF")
    model = eos.PVTModel(
        pressure_model="BM4",
        coupling="thermal-pressure",
        thermal_pressure_model="mgd",
        mgd_normalization=normalization,
    )
    constraints = (
        eos.ParameterConstraint.fixed("temperature_ref", 295.0),
        eos.ParameterConstraint.free(
            "theta_d0", 459.0, lower_bound=1.0, upper_bound=5000.0
        ),
        eos.ParameterConstraint.free(
            "gamma0", 1.5, lower_bound=0.01, upper_bound=10.0
        ),
        eos.ParameterConstraint.free(
            "q", 1.0, lower_bound=-10.0, upper_bound=10.0
        ),
    )
    request = eos.FitRequest(
        model=model,
        domain="pvt",
        target="volume",
        constraints=constraints,
        options=eos.FitOptions(
            solver_options=eos.EffectiveVarianceOptions(
                max_iterations=50,
                inner_max_iterations=10000,
            )
        ),
        request_id="naf-bm4-mgd-effective-variance",
    )
    result = eos.fit(dataset, request)
    if not result.fit.success:
        raise RuntimeError(result.fit.message)

    print("NaF BM4 + full-MGD effective-variance fit")
    print("==========================================")
    for name in ("V0", "K0", "KP", "KPP", "theta_d0", "gamma0", "q"):
        print(f"{name:10s} = {result.parameter_values[name]:.10g}")
    print(f"RMSE       = {result.fit.rmse:.10g} GPa")


if __name__ == "__main__":
    main()
