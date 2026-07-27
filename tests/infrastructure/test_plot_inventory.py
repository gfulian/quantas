"""Shared public plot-inventory contracts and first module implementations."""

from __future__ import annotations

from dataclasses import FrozenInstanceError
from pathlib import Path

import numpy as np
import pytest

from quantas.api import elasticity, ha, plotting, registry, seismic
from quantas.api.registry import Capability
from quantas.core.geometry import Hemisphere
from quantas.core.physics.seismic import SamplingLevel
from quantas.models import ResultData, ResultMetadata


ELASTICITY_DATA = Path(__file__).parents[2] / "examples" / "elasticity" / "calcite.dat"
SEISMIC_DATA = (
    Path(__file__).parents[1]
    / "physics"
    / "seismic"
    / "data"
    / "hydroxylapatite.dat"
)


def _ha_result() -> ResultData:
    payload = ha.Result(
        jobname="Inventory",
        temperature=np.asarray([0.0, 100.0, 200.0], dtype=np.float64),
        volume=np.asarray([10.0, 11.0], dtype=np.float64),
        static_energy=np.asarray([-1.0, -0.9], dtype=np.float64),
        zero_point_energy=np.asarray([[0.1, 0.2]], dtype=np.float64),
        free_energy=np.asarray(
            [[-0.9, -0.7], [-0.8, -0.6], [-0.7, -0.5]],
            dtype=np.float64,
        ),
        entropy=np.asarray(
            [[0.0, 0.0], [0.01, 0.02], [0.02, 0.03]],
            dtype=np.float64,
        ),
        isochoric_heat_capacity=np.asarray(
            [[0.0, 0.0], [0.1, 0.2], [0.2, 0.3]],
            dtype=np.float64,
        ),
        metadata={
            "units": {
                "energy": "Ha",
                "entropy": "Ha cell^-1 K^-1",
                "heat_capacity": "Ha cell^-1 K^-1",
                "volume": "A^3",
                "temperature": "K",
            }
        },
    )
    return ResultData(
        metadata=ResultMetadata(module="ha", method="harmonic"),
        results={"ha": payload},
    )


def test_inventory_descriptors_are_frozen_and_cross_validated() -> None:
    """Descriptors are immutable and reject dangling compatibility keys."""
    property_descriptor = plotting.PlotPropertyDescriptor(
        key="volume",
        name="Volume",
        symbol_math="V",
        symbol_plain="V",
        unit="A^3",
        representations=("temperature_curves",),
    )
    representation = plotting.PlotRepresentationDescriptor(
        key="temperature_curves",
        name="Temperature curves",
        plot_kind="line",
        property_keys=("volume",),
        supported_contexts=("temperature_grid",),
    )
    context = plotting.PlotContextDescriptor(
        key="temperature_grid",
        name="Temperature grid",
        values=(0.0, 300.0),
        unit="K",
        selectable=False,
    )
    inventory = plotting.PlotInventory(
        module="test",
        properties=(property_descriptor,),
        representations=(representation,),
        contexts=(context,),
    )

    assert inventory.property_by_key("volume") is property_descriptor
    assert inventory.representation_by_key("temperature_curves") is representation
    assert inventory.context_by_key("temperature_grid") is context
    with pytest.raises(FrozenInstanceError):
        property_descriptor.name = "Changed"  # type: ignore[misc]
    with pytest.raises(ValueError, match="unknown representations"):
        plotting.PlotInventory(
            module="test",
            properties=(
                plotting.PlotPropertyDescriptor(
                    key="volume",
                    name="Volume",
                    symbol_math="V",
                    symbol_plain="V",
                    unit="A^3",
                    representations=("missing",),
                ),
            ),
            representations=(),
        )
    with pytest.raises(ValueError, match="plot kind"):
        plotting.PlotRepresentationDescriptor(
            key="invalid",
            name="Invalid",
            plot_kind="plotly_surface",  # type: ignore[arg-type]
        )


def test_registry_exposes_inventory_capability_incrementally() -> None:
    """The first three module namespaces advertise the common operation."""
    for module in ("elasticity", "seismic", "ha"):
        descriptor = registry.get(module)
        assert descriptor.has(Capability.PLOT_INVENTORY)
        operation = descriptor.operation(Capability.PLOT_INVENTORY)
        assert operation.__name__ == "describe_plots"

    for module in ("qha", "thermoelasticity", "eos"):
        assert not registry.get(module).has(Capability.PLOT_INVENTORY)


@pytest.mark.elasticity
def test_elasticity_inventory_matches_builders() -> None:
    """Every advertised elasticity property is accepted by its public builder."""
    result = elasticity.run(
        ELASTICITY_DATA,
        options=elasticity.Options(calculate_2d=True, ntheta=37),
    )
    payload = elasticity.get_result(result)
    stiffness_before = np.array(payload.stiffness, copy=True)
    assert payload.properties_3d is None

    inventory = elasticity.describe_plots(result)

    assert inventory.module == "elasticity"
    assert {item.key for item in inventory.properties} == {
        "young",
        "compressibility",
        "shear",
        "poisson",
    }
    assert {item.key for item in inventory.representations} == {
        "polar_2d",
        "surface_3d",
    }
    assert inventory.context_by_key("principal_plane").values == ("xy", "xz", "yz")
    assert inventory.property_by_key("compressibility").unit == "TPa^-1"
    assert inventory.property_by_key("poisson").unit is None

    for descriptor in inventory.properties:
        if "polar_2d" in descriptor.representations:
            collection = elasticity.build_2d_plots(
                result,
                properties=(descriptor.key,),  # type: ignore[arg-type]
            )
            assert collection.plots
        if "surface_3d" in descriptor.representations:
            collection = elasticity.build_3d_plots(
                result,
                options=elasticity.SurfaceOptions(
                    ntheta=5,
                    nphi=7,
                    properties=(descriptor.key,),  # type: ignore[arg-type]
                ),
                properties=(descriptor.key,),  # type: ignore[arg-type]
            )
            assert collection.plots

    np.testing.assert_array_equal(payload.stiffness, stiffness_before)
    assert payload.properties_3d is None


@pytest.mark.seismic
def test_seismic_inventory_matches_result_level() -> None:
    """Group and enhancement fields appear only when present in the result."""
    result = seismic.run(
        SEISMIC_DATA,
        options=seismic.Options(
            ntheta=5,
            nphi=9,
            hemisphere=Hemisphere.UPPER,
            level=SamplingLevel.ENHANCEMENT,
            batch_size=10,
            track_polarization_axes=True,
        ),
    )
    inventory = seismic.describe_plots(result)
    keys = tuple(item.key for item in inventory.properties)

    assert inventory.module == "seismic"
    assert len(keys) == 16
    assert "group_v_p" in keys
    assert "power_flow_v_s1" in keys
    assert "log10_enhancement_v_s2" in keys
    assert inventory.context_by_key("sampling_level").values == ("enhancement",)
    assert inventory.context_by_key("surface_type").values == (
        "phase",
        "slowness",
        "group",
    )
    assert inventory.context_by_key("polarization_overlay").values == (False, True)

    maps = seismic.build_plots(
        result,
        options=seismic.PlotOptions(properties=keys),
    )
    assert len(maps.plots) == len(keys)
    assert not maps.warnings


@pytest.mark.seismic
def test_phase_only_seismic_inventory_excludes_unavailable_fields() -> None:
    """Discovery never advertises group or enhancement data by assumption."""
    result = seismic.run(
        SEISMIC_DATA,
        options=seismic.Options(
            ntheta=4,
            nphi=7,
            hemisphere=Hemisphere.UPPER,
            level=SamplingLevel.PHASE,
            batch_size=8,
            track_polarization_axes=False,
        ),
    )
    inventory = seismic.describe_plots(result)
    keys = {item.key for item in inventory.properties}

    assert len(keys) == 7
    assert not any(key.startswith("group_") for key in keys)
    assert not any(key.startswith("power_flow_") for key in keys)
    assert not any(key.startswith("log10_enhancement_") for key in keys)
    assert inventory.context_by_key("surface_type").values == ("phase", "slowness")
    assert inventory.context_by_key("polarization_overlay").values == (False,)
    assert inventory.warnings


@pytest.mark.ha
def test_ha_inventory_matches_grids_and_builders() -> None:
    """HA discovery reports only arrays accepted by the public curve builder."""
    result = _ha_result()
    inventory = ha.describe_plots(result)

    assert inventory.module == "ha"
    assert tuple(item.key for item in inventory.properties) == (
        "static_energy",
        "zero_point_energy",
        "entropy",
        "free_energy",
        "isochoric_heat_capacity",
    )
    assert inventory.context_by_key("temperature_grid").values == (0.0, 100.0, 200.0)
    assert inventory.context_by_key("sampled_volume").values == (10.0, 11.0)
    assert inventory.context_by_key("sampled_volume").unit == "A^3"
    assert inventory.property_by_key("zero_point_energy").symbol_plain == "U_ZP"
    assert inventory.property_by_key("isochoric_heat_capacity").unit == (
        "Ha cell^-1 K^-1"
    )

    for descriptor in inventory.properties:
        collection = ha.build_plots(result, properties=(descriptor.key,))
        assert len(collection.plots) == 1
        assert collection.plots[0].metadata["module"] == "ha"


def test_inventory_metadata_contains_no_renderer_specific_contracts() -> None:
    """Shared discovery names remain independent of frontend implementations."""
    inventory = ha.describe_plots(_ha_result())
    rendered = repr(inventory).lower()

    for forbidden in ("matplotlib", "plotly", "dash", "rich", "camera", "hover"):
        assert forbidden not in rendered
