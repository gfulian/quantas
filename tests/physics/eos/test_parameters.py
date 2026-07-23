from __future__ import annotations

import numpy as np
import pytest

from quantas.core.physics.eos import implied_kpp, parse_eos_model
from quantas.core.physics.eos.parameters import (
    resolve_energy_parameters,
    resolve_pressure_parameters,
)

K0 = 160.0
KP = 4.6


def test_birch_murnaghan_implied_parameters_follow_truncation_rules():
    bm2 = resolve_energy_parameters("BM2", [-10.0, K0, 20.0])
    bm3 = resolve_energy_parameters("BM3", [-10.0, K0, KP, 20.0])
    bm4 = resolve_energy_parameters("BM4", [-10.0, K0, KP, -0.02, 20.0])

    assert bm2.KP == pytest.approx(4.0)
    assert bm2.KPP == pytest.approx(-35.0 / (9.0 * K0))
    assert bm3.KPP == pytest.approx(-(((3.0 - KP) * (4.0 - KP)) + 35.0 / 9.0) / K0)
    assert bm4.KPP == pytest.approx(-0.02)


def test_natural_strain_implied_parameters_follow_truncation_rules():
    pt2 = resolve_pressure_parameters("PT2", [K0, 20.0])
    pt3 = resolve_pressure_parameters("PT3", [K0, KP, 20.0])

    assert pt2.KP == pytest.approx(2.0)
    assert pt2.KPP == pytest.approx(-1.0 / K0)
    delta = KP - 2.0
    assert pt3.KPP == pytest.approx(-(1.0 + delta + delta**2) / K0)


def test_vinet_uses_jeanloz_relation_with_nineteen_over_thirty_six():
    expected = -((KP / 2.0) ** 2 + KP / 2.0 - 19.0 / 36.0) / K0

    assert implied_kpp(parse_eos_model("V3"), K0, KP) == pytest.approx(expected)
    assert resolve_pressure_parameters("V2", [K0, 20.0]).KP == pytest.approx(1.0)


def test_tait_second_third_and_fourth_order_parameter_rules():
    t2 = resolve_pressure_parameters("T2", [K0, 20.0])
    t3 = resolve_pressure_parameters("T3", [K0, KP, 20.0])
    t4 = resolve_pressure_parameters("T4", [K0, KP, -0.03, 20.0])

    assert t2.KP == pytest.approx(4.0)
    assert t2.KPP == pytest.approx(-4.0 / K0)
    assert t3.KPP == pytest.approx(-KP / K0)
    assert t4.KPP == pytest.approx(-0.03)


def test_bm2_resolved_covariance_propagates_implied_parameters() -> None:
    from quantas.core.physics.eos import (
        resolved_energy_parameter_covariance,
        resolved_energy_parameter_jacobian,
    )

    parameters = np.array([-100.0, 0.55, 72.0])
    covariance = np.diag([1.0e-8, 4.0e-6, 1.0e-4])

    jacobian = resolved_energy_parameter_jacobian("BM2", parameters)
    transformed = resolved_energy_parameter_covariance("BM2", parameters, covariance)

    assert jacobian.shape == (5, 3)
    assert transformed.shape == (5, 5)
    assert transformed[2, 2] == pytest.approx(0.0)
    assert transformed[3, 3] > 0.0
    np.testing.assert_allclose(transformed, jacobian @ covariance @ jacobian.T)


def test_bm3_resolved_kpp_jacobian_matches_finite_differences() -> None:
    from quantas.core.physics.eos import resolved_energy_parameter_jacobian
    from quantas.core.physics.eos.parameters import resolve_energy_parameters

    parameters = np.array([-100.0, 0.55, 4.2, 72.0])
    jacobian = resolved_energy_parameter_jacobian("BM3", parameters)
    step = 1.0e-6

    for index in (1, 2):
        plus = parameters.copy()
        minus = parameters.copy()
        plus[index] += step
        minus[index] -= step
        numerical = (
            resolve_energy_parameters("BM3", plus).KPP
            - resolve_energy_parameters("BM3", minus).KPP
        ) / (2.0 * step)
        assert jacobian[3, index] == pytest.approx(numerical, rel=1.0e-7)
