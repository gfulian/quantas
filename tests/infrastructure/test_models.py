import numpy as np
import pytest

from quantas.models import (
    HarmonicThermodynamicResult,
    InputData,
    LinePlotSpec,
    PhononInputData,
    PlotAxis,
    PlotCollection,
    PlotSeries,
    ReportTable,
    ResultData,
    ResultMetadata,
)


def test_input_data_container():
    input_data = InputData(source="input.txt", raw="test", data={"x": 1})

    assert input_data.source == "input.txt"
    assert input_data.raw == "test"
    assert input_data.data["x"] == 1


def test_result_data_container():
    metadata = ResultMetadata(module="dummy", method="test")
    result = ResultData(
        metadata=metadata,
        options={"degree": 2},
        results={"value": 10.0},
    )

    assert result.metadata.module == "dummy"
    assert result.options["degree"] == 2
    assert result.results["value"] == 10.0


def test_phonon_input_data_exposes_shared_dimensions_and_weights():
    input_data = PhononInputData(
        natoms=2,
        formula_units=2,
        supercell=np.diag([2.0, 2.0, 1.0]),
        qpoints=2,
        volume=np.array([10.0, 11.0]),
        frequencies=np.zeros((2, 6, 2)),
        weights=np.array([1.0, 3.0]),
        qcoords=np.array([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]]),
    )

    assert input_data.nvol == 2
    assert input_data.nmodes == 6
    assert input_data.kpoints == 4
    assert input_data.total_q_points == 4.0
    assert input_data.natoms_per_formula_unit == 1.0
    np.testing.assert_allclose(input_data.normalized_weights(), [0.25, 0.75])


def test_phonon_input_data_rejects_invalid_normalization():
    with pytest.raises(ValueError, match="formula_units"):
        _ = PhononInputData(natoms=2, formula_units=0).natoms_per_formula_unit

    with pytest.raises(ValueError, match="weights"):
        PhononInputData().normalized_weights()


def test_harmonic_result_exposes_only_available_properties():
    result = HarmonicThermodynamicResult(
        zero_point_energy=np.array([[1.0]]),
        entropy=np.array([[0.1]]),
    )

    assert result.has_thermodynamic_data() is True
    assert set(result.as_property_dict()) == {"Uzp", "S"}


def test_report_table_is_a_frontend_neutral_data_container():
    """ReportTable stores tabular content without renderer-specific objects."""
    table = ReportTable(
        title="Example",
        columns=["Property", "Value"],
        rows=[["Bulk modulus", 160.0]],
        metadata={"unit": "GPa"},
    )

    assert table.title == "Example"
    assert table.columns == ["Property", "Value"]
    assert table.rows == [["Bulk modulus", 160.0]]
    assert table.metadata == {"unit": "GPa"}


def test_scientific_reports_do_not_reexport_the_shared_report_table():
    """Scientific report modules expose builders, not shared model aliases."""
    from quantas.modules.ha import report as ha_report
    from quantas.modules.qha import report as qha_report
    from quantas.modules.elasticity import report as elasticity_report

    assert not hasattr(ha_report, "ReportTable")
    assert not hasattr(qha_report, "ReportTable")
    assert not hasattr(elasticity_report, "ReportTable")


def test_plot_models_are_frontend_neutral_data_containers():
    """Plot specifications store prepared arrays without renderer objects."""
    series = PlotSeries(
        key="curve",
        label="Curve",
        x=np.array([0.0, 1.0]),
        y=np.array([2.0, 3.0]),
    )
    spec = LinePlotSpec(
        key="example",
        title="Example",
        filename_stem="example",
        x_axis=PlotAxis(key="x", label="X"),
        y_axis=PlotAxis(key="y", label="Y"),
        series=[series],
    )
    collection = PlotCollection(plots=[spec])

    assert collection.plots == [spec]
    np.testing.assert_allclose(spec.series[0].y, [2.0, 3.0])
    assert "matplotlib" not in repr(spec).lower()


def test_structure_energy_series_keeps_units_and_provenance():
    from quantas.models.computation import (
        EnergyKind,
        EnergyRecord,
        SourceProvenance,
        StructureEnergyPoint,
        StructureEnergySeries,
    )
    from quantas.models.structures import CrystalStructure

    def structure(scale: float) -> CrystalStructure:
        return CrystalStructure(
            lattice=np.identity(3) * scale,
            fractional_positions=np.array([[0.0, 0.0, 0.0]]),
            atomic_numbers=np.array([8]),
        )

    points = tuple(
        StructureEnergyPoint(
            structure=structure(scale),
            energy=EnergyRecord(value=energy, unit="Ha", kind=EnergyKind.TOTAL),
            provenance=SourceProvenance(
                interface="crystal",
                source="eos.out",
                record_index=index,
            ),
        )
        for index, (scale, energy) in enumerate(((2.0, -10.0), (2.1, -10.2)))
    )
    series = StructureEnergySeries(points=points, reference_index=1)

    assert series.npoints == 2
    assert series.energy_unit == "Ha"
    assert series.energy_kind is EnergyKind.TOTAL
    assert series.volume_unit == "angstrom^3"
    np.testing.assert_allclose(series.volumes, [8.0, 2.1**3])
    np.testing.assert_allclose(series.energies, [-10.0, -10.2])
    assert series.points[1].provenance is not None
    assert series.points[1].provenance.record_index == 1


def test_structure_energy_series_rejects_mixed_units():
    from quantas.models.computation import (
        EnergyKind,
        EnergyRecord,
        StructureEnergyPoint,
        StructureEnergySeries,
    )
    from quantas.models.structures import CrystalStructure

    structure = CrystalStructure(
        lattice=np.identity(3),
        fractional_positions=np.array([[0.0, 0.0, 0.0]]),
        atomic_numbers=np.array([8]),
    )
    points = (
        StructureEnergyPoint(
            structure=structure,
            energy=EnergyRecord(value=-1.0, unit="Ha", kind=EnergyKind.DFT),
        ),
        StructureEnergyPoint(
            structure=structure,
            energy=EnergyRecord(value=-27.2, unit="eV", kind=EnergyKind.DFT),
        ),
    )

    with pytest.raises(ValueError, match="inconsistent energy units"):
        StructureEnergySeries(points=points)
