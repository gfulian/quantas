"""Shared elastic-state and Kieffer cutoff contracts."""

from __future__ import annotations

import numpy as np
import pytest

from quantas.models import (
    CutoffVolumeSource,
    ElasticState,
    ElasticStateSeries,
    ElasticTensorKind,
    KiefferCutoffState,
    KiefferVolumeSeries,
    PressureSource,
    PrestressProvenance,
    VolumeMatchPolicy,
    match_sampled_volumes,
)


def _prestress(kind: ElasticTensorKind) -> PrestressProvenance:
    """Return minimal valid provenance for a tensor convention."""
    return PrestressProvenance(
        tensor_kind=kind,
        pressure_gpa=1.0,
        pressure_source=PressureSource.APPLIED_PRESTRESS,
    )


def _state(volume: float, kind: ElasticTensorKind) -> ElasticState:
    """Return one elastic state with a cubic lattice."""
    edge = volume ** (1.0 / 3.0)
    return ElasticState(
        volume=volume,
        density=3000.0,
        stiffness=np.identity(6) * 100.0,
        prestress=_prestress(kind),
        energy=-10.0,
        energy_unit="Ha",
        lattice=np.identity(3) * edge,
        source="elastic.out",
    )


def test_incremental_provenance_prevents_double_correction() -> None:
    provenance = PrestressProvenance(
        tensor_kind=ElasticTensorKind.WALLACE_HYDROSTATIC,
        pressure_gpa=2.0,
        pressure_source=PressureSource.ENERGY_EOS,
        correction_method="Barron-Klein hydrostatic",
        correction_applied_by="quantas",
        source_tensor_kind=ElasticTensorKind.RAW_ENERGY_STRAIN,
    )
    provenance.require_incremental()

    with pytest.raises(ValueError, match="second time"):
        PrestressProvenance(
            tensor_kind=ElasticTensorKind.WALLACE_HYDROSTATIC,
            pressure_gpa=2.0,
            pressure_source=PressureSource.ENERGY_EOS,
            correction_method="duplicate",
            correction_applied_by="quantas",
            source_tensor_kind=ElasticTensorKind.WALLACE_HYDROSTATIC,
        )


def test_raw_tensor_is_retained_but_rejected_for_acoustics() -> None:
    raw = _state(100.0, ElasticTensorKind.RAW_ENERGY_STRAIN)
    series = ElasticStateSeries(states=(raw,), reference_index=0)
    with pytest.raises(ValueError, match="incremental"):
        series.require_incremental()


def test_elastic_series_exposes_float64_arrays() -> None:
    series = ElasticStateSeries(
        states=(
            _state(100.0, ElasticTensorKind.WALLACE_HYDROSTATIC),
            _state(110.0, ElasticTensorKind.WALLACE_HYDROSTATIC),
        ),
        reference_index=0,
    )
    series.require_incremental()
    np.testing.assert_allclose(series.volumes, [100.0, 110.0])
    assert series.densities.dtype == np.dtype("float64")
    assert series.stiffness.shape == (2, 6, 6)


def test_elastic_series_rejects_order_and_energy_mismatch() -> None:
    first = _state(100.0, ElasticTensorKind.WALLACE_HYDROSTATIC)
    second = _state(90.0, ElasticTensorKind.WALLACE_HYDROSTATIC)
    with pytest.raises(ValueError, match="increasing"):
        ElasticStateSeries(states=(first, second), reference_index=0)

    second = _state(110.0, ElasticTensorKind.WALLACE_HYDROSTATIC)
    second.energy_unit = "eV"
    with pytest.raises(ValueError, match="energy units"):
        ElasticStateSeries(states=(first, second), reference_index=0)


def test_cutoff_series_uses_thermodynamic_array_orientation() -> None:
    states = tuple(
        KiefferCutoffState(
            volume=volume,
            frequencies_hz=np.array([3.0e12, 4.0e12, 5.0e12]) * factor,
            effective_velocities_km_s=np.array([4.0, 5.0, 8.0]) * factor,
            source=CutoffVolumeSource.DIRECT,
            source_elastic_indices=(index,),
        )
        for index, (volume, factor) in enumerate(((100.0, 1.0), (110.0, 0.9)))
    )
    series = KiefferVolumeSeries(states=states)
    assert series.frequencies_hz.shape == (3, 2)
    np.testing.assert_allclose(series.volumes, [100.0, 110.0])
    assert not states[0].frequencies_hz.flags.writeable


def test_interpolated_cutoff_requires_complete_provenance() -> None:
    kwargs = {
        "volume": 105.0,
        "frequencies_hz": [3.0e12, 4.0e12, 5.0e12],
        "effective_velocities_km_s": [4.0, 5.0, 8.0],
        "source": CutoffVolumeSource.INTERPOLATED,
        "source_elastic_indices": (0, 1),
    }
    state = KiefferCutoffState(**kwargs)
    with pytest.raises(ValueError, match="interpolation method"):
        KiefferVolumeSeries(states=(state,))
    series = KiefferVolumeSeries(states=(state,), interpolation_method="pchip")
    assert series.interpolation_method == "pchip"


def test_volume_matching_records_differences_and_order() -> None:
    policy = VolumeMatchPolicy(relative_tolerance=0.0, absolute_tolerance=1.0e-5)
    matches = match_sampled_volumes(
        [100.0, 110.0], [110.000001, 90.0, 100.000001], policy=policy
    )
    assert [match.source_index for match in matches] == [2, 0]
    assert matches[0].absolute_difference == pytest.approx(1.0e-6)
    assert matches[0].relative_difference == pytest.approx(1.0e-8)


def test_volume_matching_rejects_missing_and_ambiguous_sources() -> None:
    with pytest.raises(ValueError, match="no source match"):
        match_sampled_volumes([100.0], [101.0])
    with pytest.raises(ValueError, match="multiple source matches"):
        match_sampled_volumes(
            [100.0],
            [99.9999999, 100.0000001],
            policy=VolumeMatchPolicy(absolute_tolerance=1.0e-5),
        )
