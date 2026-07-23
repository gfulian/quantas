"""Compare QHA interpolation and volume-minimization choices at one P,T state."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
from scipy.constants import Avogadro, physical_constants

from quantas.api import qha

_HARTREE_TO_J_MOL = physical_constants["Hartree energy"][0] * Avogadro


def scalar(value: np.ndarray) -> float:
    """Return the only numerical value in a one-state result array."""
    return float(np.asarray(value, dtype=np.float64).reshape(-1)[0])


def run_case(source: Path, scheme: str, minimization: str, eos: str) -> dict[str, object]:
    """Run one QHA method combination and return selected observables."""
    options = qha.Options(
        temperature_min=300.0,
        temperature_max=300.0,
        temperature_step=1.0,
        pressure_min=0.0,
        pressure_max=0.0,
        pressure_step=1.0,
        scheme=scheme,  # type: ignore[arg-type]
        minimization=minimization,  # type: ignore[arg-type]
        eos=eos,
        energy_degree=3,
        free_energy_degree=3,
        frequency_degree=3,
        thermal_expansion_method="mixed_derivative",
        calculate_gruneisen=True,
        calculate_mode_gruneisen=scheme == "freq",
        estimate_uncertainties=False,
    )
    payload = qha.get_result(qha.run(source, options=options))
    return {
        "scheme": scheme,
        "minimization": minimization,
        "eos": eos if minimization == "eos" else "-",
        "V_A3": scalar(payload.equilibrium_volume),
        "KT_GPa": scalar(payload.isothermal_bulk_modulus),
        "KS_GPa": scalar(payload.adiabatic_bulk_modulus),
        "alphaV_1e5_K-1": 1.0e5 * scalar(payload.thermal_expansion),
        "Cv_J_mol-1_K-1": _HARTREE_TO_J_MOL * scalar(payload.isochoric_heat_capacity),
        "Cp_J_mol-1_K-1": _HARTREE_TO_J_MOL * scalar(payload.isobaric_heat_capacity),
        "gamma": scalar(payload.gruneisen),
    }


def main() -> None:
    """Run the four core combinations and write a comparison CSV file."""
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path, default=Path("qha_method_comparison.csv"))
    args = parser.parse_args()

    cases = (
        ("freq", "poly", "BM3"),
        ("td", "poly", "BM3"),
        ("freq", "eos", "BM3"),
        ("td", "eos", "BM3"),
    )
    rows = [run_case(args.input, *case) for case in cases]
    with args.output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
                stream, 
                fieldnames=list(rows[0]),
                lineterminator="\n",
                )
        writer.writeheader()
        writer.writerows(rows)

    for row in rows:
        print(
            f"{row['scheme']:>4s} {row['minimization']:>4s} "
            f"V={row['V_A3']:.6f} A^3  KT={row['KT_GPa']:.4f} GPa  "
            f"alphaV={row['alphaV_1e5_K-1']:.6f} x10^-5 K^-1"
        )
    print(f"Written: {args.output}")


if __name__ == "__main__":
    main()
