"""Fit the MgO static E-V example with a third-order Birch-Murnaghan EOS."""

from __future__ import annotations

from pathlib import Path

from quantas.api import eos

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "EV_mgo_pbe.dat"


def main() -> None:
    """Run one public E-V request and print the equilibrium parameters."""
    dataset = eos.read_input(DATA)
    request = eos.FitRequest(
        model="BM3",
        domain="ev",
        target="energy",
        options=eos.FitOptions(solver_options=eos.OLSOptions()),
        request_id="mgo-pbe-bm3-energy",
    )
    result = eos.fit(dataset, request)
    if not result.fit.success:
        raise RuntimeError(result.fit.message)

    print("MgO PBE BM3 energy-volume fit")
    print("=============================")
    for name in ("E0", "V0", "K0", "KP", "KPP"):
        print(f"{name:4s} = {result.parameter_values[name]:.12g}")
    print(f"RMSE = {result.fit.rmse:.8e} Ha")
    print(f"Pressure relation = {result.metadata['pressure_relation']}")


if __name__ == "__main__":
    main()
