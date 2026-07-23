"""Frozen real-data validation matrix for thermoelastic QSA.

The compact fixture was extracted from user-supplied CRYSTAL calculations that
used ``PRESSURE`` and from native Quantas QHA HDF5 archives for cubic MgO and
trigonal dolomite.  The large raw calculations are intentionally not shipped.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from quantas.core.physics.elasticity import cold_finite_strain_component
from quantas.core.physics.units import energy_to_pressure


_FIXTURE = Path(__file__).with_name("data") / "thermoelastic_reference.json"


def _data() -> dict[str, object]:
    return json.loads(_FIXTURE.read_text(encoding="utf-8"))


@pytest.mark.parametrize("system", ("mgo", "dolomite"))
def test_cold_strain_fits_match_linear_reference(
    system: str,
) -> None:
    """A NumPy linear solve independently reproduces every real QSA fit."""
    block = _data()[system]["fit"]
    eos = block["reference_eos_resolved"]
    v0 = float(eos["V0"])
    k0 = float(eos["K0"]) * float(energy_to_pressure(1.0, "Ha", "A", "GPa"))
    kp = float(eos["KP"])
    for label, component in block["components"].items():
        volumes = np.asarray(component["volumes_A3"], dtype=np.float64)
        observed = np.asarray(component["observed_GPa"], dtype=np.float64)
        delta = float(component["wallace_delta"])
        base = cold_finite_strain_component(
            volumes,
            reference_volume=v0,
            bulk_modulus=k0,
            bulk_modulus_derivative=kp,
            reference_component=0.0,
            component_pressure_derivative=0.0,
            wallace_delta=delta,
            order=3,
        )
        column_c0 = (
            cold_finite_strain_component(
                volumes,
                reference_volume=v0,
                bulk_modulus=k0,
                bulk_modulus_derivative=kp,
                reference_component=1.0,
                component_pressure_derivative=0.0,
                wallace_delta=delta,
                order=3,
            )
            - base
        )
        column_cp = (
            cold_finite_strain_component(
                volumes,
                reference_volume=v0,
                bulk_modulus=k0,
                bulk_modulus_derivative=kp,
                reference_component=0.0,
                component_pressure_derivative=1.0,
                wallace_delta=delta,
                order=3,
            )
            - base
        )
        solved = np.linalg.lstsq(
            np.column_stack((column_c0, column_cp)),
            observed - base,
            rcond=None,
        )[0]
        np.testing.assert_allclose(
            solved,
            np.asarray(component["expected_parameters"], dtype=np.float64),
            rtol=2.0e-8,
            atol=2.0e-8,
            err_msg=f"real-data regression failed for {system} {label}",
        )


@pytest.mark.parametrize("system", ("mgo", "dolomite"))
def test_adiabatic_states_follow_rank_one_update(
    system: str,
) -> None:
    """Frozen states satisfy the independent anisotropic C_S-C_T identity."""
    for state in _data()[system]["pt"]["states"]:
        ct = np.asarray(state["stiffness_isothermal"], dtype=np.float64)
        cs = np.asarray(state["stiffness_adiabatic"], dtype=np.float64)
        correction = np.asarray(state["adiabatic_correction"], dtype=np.float64)
        temperature = float(state["temperature_K"])
        if temperature == 0.0:
            np.testing.assert_array_equal(cs, ct)
            np.testing.assert_array_equal(correction, np.zeros((6, 6)))
            continue
        alpha = np.asarray(state["thermal_expansion_tensor"], dtype=np.float64)
        alpha_voigt = np.asarray(
            [
                alpha[0, 0],
                alpha[1, 1],
                alpha[2, 2],
                2.0 * alpha[1, 2],
                2.0 * alpha[0, 2],
                2.0 * alpha[0, 1],
            ]
        )
        beta_pa = (ct * 1.0e9) @ alpha_voigt
        factor = (
            temperature
            * float(state["equilibrium_volume"])
            * 1.0e-30
            / float(state["isochoric_heat_capacity_cell"])
        )
        independently_calculated = factor * np.outer(beta_pa, beta_pa) / 1.0e9
        np.testing.assert_allclose(
            correction,
            independently_calculated,
            rtol=2.0e-13,
            atol=2.0e-12,
        )
        np.testing.assert_allclose(cs, ct + correction, rtol=2.0e-13, atol=2.0e-12)
        eigenvalues = np.linalg.eigvalsh(correction)
        assert eigenvalues.min() >= -1.0e-10
        assert np.count_nonzero(eigenvalues > 1.0e-8) <= 1


def test_real_mgo_adiabatic_update_preserves_cubic_symmetry() -> None:
    """Hydrostatic cubic MgO changes C11 and C12 equally and not C44."""
    for state in _data()["mgo"]["pt"]["states"]:
        correction = np.asarray(state["adiabatic_correction"], dtype=np.float64)
        np.testing.assert_allclose(correction[0, 0], correction[0, 1], atol=1.0e-12)
        np.testing.assert_allclose(correction[3:, 3:], 0.0, atol=1.0e-12)


def test_real_validation_domains_remain_mechanically_stable() -> None:
    """Both supplied P-T domains retain positive Wallace stability margins."""
    data = _data()
    assert data["mgo"]["pt"]["global"]["minimum_stability_eigenvalue_GPa"] > 100.0
    assert data["dolomite"]["pt"]["global"]["minimum_stability_eigenvalue_GPa"] > 25.0


def test_dolomite_frame_normalization_is_material() -> None:
    """The low-symmetry reference records a finite frame-migration effect."""
    frame = _data()["dolomite"]["frame_normalization"]
    assert frame["maximum_removed_rotation_degrees"] > 1.0
    assert frame["maximum_ordered_atom_displacement_A"] < 0.05
    assert frame["raw_vs_normalized_grid_max_abs_GPa"] > 1.0
