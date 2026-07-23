"""Fit the rutile V-T example with the quadratic Berman model."""

from __future__ import annotations

from pathlib import Path

from quantas.api import eos

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "VT_rutile.dat"


def main() -> None:
    """Run one public V-T request and calculate representative states."""
    dataset = eos.read_input(DATA)
    request = eos.FitRequest(
        model="berman:quadratic",
        domain="vt",
        target="volume",
        options=eos.FitOptions(
            solver_options=eos.EffectiveVarianceOptions(
                max_iterations=30,
                inner_max_iterations=5000,
            )
        ),
        request_id="rutile-berman-effective-variance",
    )
    result = eos.fit(dataset, request)
    if not result.fit.success:
        raise RuntimeError(result.fit.message)

    print("Rutile quadratic-Berman effective-variance fit")
    print("================================================")
    for name in ("V0", "temperature_ref", "alpha0", "alpha1"):
        print(f"{name:16s} = {result.parameter_values[name]:.10g}")
    print(f"RMSE             = {result.fit.rmse:.10g}")


if __name__ == "__main__":
    main()
