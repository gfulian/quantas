"""Generate the frozen Quantas 0.9 SEISMIC characterization baseline."""

from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np

PROJECT_ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from tests.reference.seismic_reference import SeismicFormulaReference  # noqa: E402


def read_fixture(path: Path) -> tuple[np.ndarray, float]:
    """Read the triangular stiffness fixture and density.

    Parameters
    ----------
    path : pathlib.Path
        Fixture path.

    Returns
    -------
    tuple
        Stiffness matrix in GPa and density in kg m-3.
    """
    lines = path.read_text(encoding="utf-8").splitlines()
    stiffness = np.zeros((6, 6), dtype=float)
    for row, line in enumerate(lines[1:7]):
        values = np.asarray([float(value) for value in line.split()], dtype=float)
        stiffness[row, row:] = values
    stiffness += np.triu(stiffness, 1).T
    return stiffness, float(lines[7])


def main() -> None:
    """Generate NPZ arrays and their JSON provenance manifest."""
    fixture = (
        PROJECT_ROOT
        / "tests"
        / "physics"
        / "seismic"
        / "data"
        / "hydroxylapatite.dat"
    )
    stiffness, density = read_fixture(fixture)
    directions = np.asarray(
        [
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 1.0, 0.0],
            [1.0, 2.0, 3.0],
            [2.0, -1.0, 3.0],
            [-2.0, 3.0, 1.0],
            [1.0, -3.0, 2.0],
            [-3.0, -2.0, 1.0],
        ],
        dtype=float,
    )
    directions /= np.linalg.norm(directions, axis=1)[:, None]
    reference = SeismicFormulaReference(stiffness, density)
    results = [reference.solve(direction) for direction in directions]

    arrays = {
        "directions": directions,
        "christoffel_hessian": reference.christoffel_hessian,
        "christoffel": np.stack([result.christoffel for result in results]),
        "christoffel_gradient": np.stack(
            [result.christoffel_gradient for result in results]
        ),
        "eigenvalues": np.stack([result.eigenvalues for result in results]),
        "polarizations": np.stack([result.polarizations for result in results]),
        "phase_speeds": np.stack([result.phase_speeds for result in results]),
        "eigenvalue_gradients": np.stack(
            [result.eigenvalue_gradients for result in results]
        ),
        "eigenvalue_hessians": np.stack(
            [result.eigenvalue_hessians for result in results]
        ),
        "group_velocities": np.stack([result.group_velocities for result in results]),
        "group_speeds": np.stack([result.group_speeds for result in results]),
        "group_directions": np.stack([result.group_directions for result in results]),
        "power_flow_angles": np.stack([result.power_flow_angles for result in results]),
        "enhancement": np.stack([result.enhancement for result in results]),
        "log10_enhancement": np.log10(
            np.stack([result.enhancement for result in results])
        ),
    }
    baseline_dir = PROJECT_ROOT / "tests" / "baselines"
    np.savez(  # type: ignore[arg-type]
        baseline_dir / "seismic_reference.npz",
        **arrays,  # type: ignore[arg-type]
    )

    metadata = {
        "name": "Quantas 0.9 SEISMIC formula reference",
        "source": {
            "repository": "https://github.com/gfulian/quantas",
            "implementation": "quantas/seismic/utils/seismic_obj.py",
            "method": "Jaeken-Cottenier analytical Christoffel derivatives",
        },
        "material": {
            "name": "Hydroxylapatite",
            "stiffness_unit": "GPa",
            "density": density,
            "density_unit": "kg m^-3",
        },
        "mode_order": ["v_s2", "v_s1", "v_p"],
        "mode_aliases": {
            "v_s2": "slow_secondary",
            "v_s1": "fast_secondary",
            "v_p": "primary",
        },
        "units": {
            "eigenvalues": "km^2 s^-2",
            "phase_speeds": "km s^-1",
            "group_velocities": "km s^-1",
            "group_speeds": "km s^-1",
            "power_flow_angles": "rad",
            "enhancement": "dimensionless",
            "log10_enhancement": "dimensionless",
        },
        "density_scaling": {
            "reduced_stiffness": "C_GPa * 1000 / rho",
            "reason": "10^9 Pa per GPa divided by 10^6 m^2 per km^2",
        },
        "known_reference_limitations": [
            "Cartesian azimuth uses atan2(x, y) instead of atan2(y, x).",
            "The HDF5 header labels column 14 as log10(A), but raw A is stored.",
            "The original Quantas 0.9 2D polarization plot reads group-vector columns.",
            "No explicit degeneracy mask or polarization sign tracking is present.",
        ],
        "tolerance": {"rtol": 1.0e-11, "atol": 1.0e-11},
    }
    (baseline_dir / "seismic_reference.json").write_text(
        json.dumps(metadata, indent=2) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
