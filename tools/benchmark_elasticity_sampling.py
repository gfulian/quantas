"""Benchmark pointwise optimization against exact batched elasticity fields."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from time import perf_counter

import numpy as np
from scipy import optimize

PROJECT_ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from quantas.core.physics.elasticity import (  # noqa: E402
    ElasticTensor,
    poisson_ratio,
    sample_elastic_directional_field,
    sample_elasticity_surfaces,
    shear_modulus,
)
from quantas.api import elasticity  # noqa: E402


def _reference_transverse_extrema(
    tensor: ElasticTensor,
    theta: np.ndarray,
    phi: np.ndarray,
) -> None:
    """Reproduce the former per-direction Powell sampling workload."""
    options = {"xtol": 0.01, "ftol": 0.001}
    for polar, azimuth in zip(theta, phi, strict=True):
        for function in (
            lambda chi: shear_modulus(tensor, (polar, azimuth, chi)),
            lambda chi: poisson_ratio(tensor, (polar, azimuth, chi)),
        ):
            optimize.minimize(
                function,
                np.pi / 2.0,
                method="Powell",
                options=options,
            )
            optimize.minimize(
                lambda chi, objective=function: -objective(chi),
                np.pi / 2.0,
                method="Powell",
                options=options,
            )


def _elapsed(operation) -> float:
    """Return elapsed wall time for a zero-argument operation."""
    start = perf_counter()
    operation()
    return perf_counter() - start


def main() -> None:
    """Run reproducible two- and three-dimensional sampling benchmarks."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--points-2d", type=int, default=361)
    parser.add_argument("--ntheta", type=int, default=61)
    parser.add_argument("--nphi", type=int, default=121)
    parser.add_argument("--skip-reference", action="store_true")
    arguments = parser.parse_args()

    fixture = (
        PROJECT_ROOT
        / "tests"
        / "physics"
        / "elasticity"
        / "data"
        / "hydroxylapatite.dat"
    )
    input_data = elasticity.read_input(fixture)
    tensor = ElasticTensor(input_data.stiffness)

    phi = np.linspace(0.0, 2.0 * np.pi, arguments.points_2d, endpoint=True)
    theta = np.full_like(phi, np.pi / 2.0)
    if not arguments.skip_reference:
        reference = _elapsed(lambda: _reference_transverse_extrema(tensor, theta, phi))
        print(f"Pointwise 2D reference: {reference:.6f} s")
    else:
        reference = None

    batch_2d = _elapsed(lambda: sample_elastic_directional_field(tensor, theta, phi))
    batch_3d = _elapsed(
        lambda: sample_elasticity_surfaces(
            tensor,
            ntheta=arguments.ntheta,
            nphi=arguments.nphi,
        )
    )
    print(f"Exact batched 2D field:          {batch_2d:.6f} s")
    print(f"Exact batched 3D surfaces:       {batch_3d:.6f} s")
    if reference is not None:
        print(f"2D speed-up:                     {reference / batch_2d:.2f}x")


if __name__ == "__main__":
    main()
