"""Fit the quartz P-V example with BM3 effective variance."""

from __future__ import annotations

from pathlib import Path

from quantas.api import eos

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "PV_quartz.dat"


def main() -> None:
    """Run one public P-V request and print the principal results."""
    dataset = eos.read_input(DATA)
    request = eos.FitRequest(
        model="BM3",
        domain="pv",
        target="volume",
        options=eos.FitOptions(
            solver_options=eos.EffectiveVarianceOptions(
                max_iterations=50,
                inner_max_iterations=10000,
            )
        ),
        request_id="quartz-bm3-effective-variance",
    )
    result = eos.fit(dataset, request)
    if not result.fit.success:
        raise RuntimeError(result.fit.message)

    print("Quartz BM3 effective-variance fit")
    print("=================================")
    for name in ("V0", "K0", "KP", "KPP"):
        print(f"{name:4s} = {result.parameter_values[name]:.9g}")
    print(f"RMSE = {result.fit.rmse:.8f} GPa")
    diagnostics = result.fit.diagnostics
    if diagnostics is not None:
        print(f"Reduced chi-square = {diagnostics.reduced_chi_square:.8f}")


if __name__ == "__main__":
    main()
