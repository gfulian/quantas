"""Derive MgO axial response from an Energy EOS and structural path."""

from __future__ import annotations

from pathlib import Path

from quantas.api import eos

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "EV_mgo_pbe_crystallographic.dat"


def main() -> None:
    """Fit SJEOS, derive the primary axial response, and add a BM3 axial fit."""
    dataset = eos.read_input(DATA)
    request = eos.FitRequest(
        model="SJ",
        axial_model="BM3",
        domain="ev",
        target="energy",
        options=eos.FitOptions(solver_options=eos.OLSOptions()),
        request_id="mgo-pbe-sj-structural",
    )
    result = eos.fit(dataset, request)
    if not result.fit.success:
        raise RuntimeError(result.fit.message)

    print("MgO theoretical axial response")
    print("==============================")
    print(f"a0    = {result.derived['a0']:.9f} angstrom")
    print(f"eta_a = {result.derived['eta_a']:.9f}")
    print(f"M_a   = {result.derived['M_a']:.6f} GPa")

    secondary = result.metadata["secondary_axial_fits"]["fits"]["a"]
    values = secondary["parameter_values"]
    print("\nSecondary BM3 P(a^3) parameterization")
    for name in ("M0", "MP", "MPP", "L0"):
        print(f"{name:4s} = {values[name]:.12g}")


if __name__ == "__main__":
    main()
