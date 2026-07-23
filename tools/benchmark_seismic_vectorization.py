"""Benchmark scalar reference and batched SEISMIC kernels."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from time import perf_counter

import numpy as np

PROJECT_ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from tests.reference.seismic_reference import SeismicFormulaReference  # noqa: E402
from tests.reference.seismic_vectorized import solve_batched  # noqa: E402


def read_fixture(path: Path) -> tuple[np.ndarray, float]:
    """Read a triangular Quantas elasticity fixture.

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


def create_grid(ntheta: int, nphi: int) -> np.ndarray:
    """Create an upper-hemisphere grid without a duplicated azimuth seam.

    Parameters
    ----------
    ntheta : int
        Number of polar samples.
    nphi : int
        Number of azimuthal samples.

    Returns
    -------
    numpy.ndarray
        Cartesian directions with shape ``(ntheta * nphi, 3)``.
    """
    theta = np.linspace(0.0, 0.5 * np.pi, ntheta)
    phi = np.linspace(0.0, 2.0 * np.pi, nphi, endpoint=False)
    theta_grid, phi_grid = np.meshgrid(theta, phi, indexing="ij")
    return np.stack(
        [
            np.sin(theta_grid) * np.cos(phi_grid),
            np.sin(theta_grid) * np.sin(phi_grid),
            np.cos(theta_grid),
        ],
        axis=-1,
    ).reshape(-1, 3)


def elapsed(callable_) -> float:
    """Return elapsed wall time for a zero-argument callable.

    Parameters
    ----------
    callable_ : callable
        Operation to time.

    Returns
    -------
    float
        Elapsed seconds.
    """
    start = perf_counter()
    callable_()
    return perf_counter() - start


def main() -> None:
    """Run the reproducible scalar-versus-batched benchmark."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--ntheta", type=int, default=91)
    parser.add_argument("--nphi", type=int, default=181)
    parser.add_argument("--skip-scalar", action="store_true")
    arguments = parser.parse_args()

    fixture = (
        PROJECT_ROOT
        / "tests"
        / "physics"
        / "seismic"
        / "data"
        / "hydroxylapatite.dat"
    )
    stiffness, density = read_fixture(fixture)
    directions = create_grid(arguments.ntheta, arguments.nphi)
    reference = SeismicFormulaReference(stiffness, density)

    print(f"Grid points: {directions.shape[0]}")
    if not arguments.skip_scalar:
        scalar_time = elapsed(lambda: [reference.solve(q) for q in directions])
        print(f"Scalar formula reference: {scalar_time:.6f} s")
    else:
        scalar_time = None

    phase_group_time = elapsed(
        lambda: solve_batched(
            stiffness,
            density,
            directions,
            calculate_enhancement=False,
        )
    )
    complete_time = elapsed(
        lambda: solve_batched(
            stiffness,
            density,
            directions,
            calculate_enhancement=True,
        )
    )
    print(f"Batched phase/group:      {phase_group_time:.6f} s")
    print(f"Batched complete:         {complete_time:.6f} s")
    if scalar_time is not None:
        print(f"Complete speed-up:        {scalar_time / complete_time:.2f}x")


if __name__ == "__main__":
    main()
