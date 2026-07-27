"""Session-aware public EOS plot inventory and CLI discovery."""

from __future__ import annotations

from pathlib import Path

from click.testing import CliRunner
import numpy as np

from quantas.api import eos, plotting
from quantas.cli.main import main
from quantas.core.physics.eos import PVTEOS


DATA = Path(__file__).with_name("data")


def _fit_request(target: str = "volume") -> eos.FitRequest:
    """Return a deterministic P-V request for inventory tests."""
    return eos.FitRequest(
        model="BM3",
        domain="pv",
        target=target,
        options=eos.FitOptions(solver_options=eos.OLSOptions()),
    )


def _quartz_archive(tmp_path: Path, *, accept: bool = True) -> Path:
    """Create one successful single-slot archive."""
    dataset = eos.read_input(DATA / "PV_quartz.dat")
    request = _fit_request()
    result = eos.fit(dataset, request)
    path = tmp_path / "quartz_inventory.hdf5"
    with eos.Archive.create(path, dataset=dataset) as archive:
        archive.append_fit(1, request, result, accept=accept)
    return path


def _topaz_archive(tmp_path: Path) -> Path:
    """Create an archive with two independently accepted P-V slots."""
    dataset = eos.read_input(DATA / "PV_topaz.dat")
    path = tmp_path / "topaz_inventory.hdf5"
    with eos.Archive.create(path, dataset=dataset) as archive:
        for target in ("volume", "a"):
            request = _fit_request(target)
            result = eos.fit(dataset, request)
            archive.append_fit(1, request, result, accept=True)
    return path


def _vt_archive(tmp_path: Path) -> Path:
    """Create one accepted V-T archive through the public API."""
    dataset = eos.read_input(DATA / "rutile.dat")
    request = eos.FitRequest(
        model="berman:quadratic",
        domain="vt",
        target="volume",
        options=eos.FitOptions(
            solver_options=eos.EffectiveVarianceOptions(
                max_iterations=30,
                inner_max_iterations=5000,
            )
        ),
    )
    result = eos.fit(dataset, request)
    path = tmp_path / "rutile_inventory.hdf5"
    with eos.Archive.create(path, dataset=dataset) as archive:
        archive.append_fit(1, request, result, accept=True)
    return path


def _pvt_archive(tmp_path: Path) -> Path:
    """Create one accepted synthetic P-V-T archive."""
    model = eos.PVTModel("BM3", "linear", "berman:quadratic")
    pressure_parameters = {"K0": 160.0, "KP": 4.2, "KPP": -0.02, "V0": 100.0}
    thermal_parameters = {
        "V0": 100.0,
        "temperature_ref": 300.0,
        "alpha0": 3.0e-5,
        "alpha1": 1.0e-8,
    }
    coupling_parameters = {"dK0_dT": -0.02}
    temperature = np.repeat([300.0, 600.0, 900.0], 6)
    pressure = np.tile(np.linspace(0.0, 10.0, 6), 3)
    volume = PVTEOS().volume(
        model,
        pressure_parameters,
        thermal_parameters,
        coupling_parameters,
        pressure,
        temperature,
    )
    dataset = eos.Dataset(
        jobname="Synthetic PVT inventory",
        columns={
            "pressure": pressure,
            "temperature": temperature,
            "volume": volume,
            "sigma_pressure": np.full(pressure.size, 0.01),
            "sigma_temperature": np.full(pressure.size, 1.0),
            "sigma_volume": np.full(pressure.size, 0.005),
        },
        units={
            "pressure": "GPa",
            "temperature": "K",
            "volume": "angstrom^3",
            "sigma_pressure": "GPa",
            "sigma_temperature": "K",
            "sigma_volume": "angstrom^3",
        },
    )
    request = eos.FitRequest(
        model=model,
        domain="pvt",
        constraints=(
            eos.ParameterConstraint.free("K0", 158.0),
            eos.ParameterConstraint.free("KP", 4.0),
            eos.ParameterConstraint.free("V0", 99.8),
            eos.ParameterConstraint.fixed("temperature_ref", 300.0),
            eos.ParameterConstraint.free("alpha0", 2.8e-5),
            eos.ParameterConstraint.free("alpha1", 0.8e-8),
            eos.ParameterConstraint.free("dK0_dT", -0.018),
        ),
        options=eos.FitOptions(
            solver_options=eos.OLSOptions(max_iterations=5000)
        ),
    )
    result = eos.fit(dataset, request)
    path = tmp_path / "pvt_inventory.hdf5"
    with eos.Archive.create(path, dataset=dataset) as archive:
        archive.append_fit(1, request, result, accept=True)
    return path


def test_eos_inventory_keeps_dataset_only_archive_browsable(tmp_path: Path) -> None:
    """An archive without fit attempts still exposes its embedded dataset."""
    dataset = eos.read_input(DATA / "PV_quartz.dat")
    path = tmp_path / "dataset_only_inventory.hdf5"
    with eos.Archive.create(path, dataset=dataset):
        pass

    inventory = eos.describe_plots(path)

    assert inventory.selected_record_id is None
    assert inventory.selected_plots is None
    assert inventory.records == ()
    assert len(inventory.slots) == 1
    assert inventory.slots[0].key == "pv/volume"
    assert inventory.slots[0].status is eos.SlotStatus.NOT_PROCESSED
    assert inventory.slots[0].record_ids == ()
    assert inventory.dataset_by_id(1).slot_keys == ("pv/volume",)
    assert any("no accepted record" in item for item in inventory.warnings)


def test_eos_inventory_describes_archive_and_selected_record(tmp_path: Path) -> None:
    """A unique accepted record exposes archive and common plot metadata."""
    path = _quartz_archive(tmp_path)

    before = path.stat().st_size
    inventory = eos.describe_plots(path)
    after = path.stat().st_size

    assert isinstance(inventory, eos.ArchivePlotInventory)
    assert inventory.path == path
    assert inventory.schema_version == "1.1"
    assert inventory.event_count == 2
    assert inventory.selected_record_id == 1
    assert before == after

    dataset = inventory.dataset_by_id(1)
    assert dataset.jobname == "PV_Quartz"
    assert dataset.npoints == 22
    assert dataset.selected_npoints == 22
    assert dataset.excluded_npoints == 0
    assert dataset.slot_keys == ("pv/volume",)
    assert dataset.unit_for("pressure") == "GPa"
    assert dataset.unit_for("volume") == "angstrom^3"

    slot = inventory.slot_by_key("pv/volume")
    assert slot.status is eos.SlotStatus.ACCEPTED
    assert slot.accepted_record_id == 1
    assert slot.record_ids == (1,)
    assert slot.plottable_record_ids == (1,)

    record = inventory.record_by_id(1)
    assert record.disposition is eos.RecordDisposition.ACCEPTED
    assert record.current_accepted
    assert record.plottable
    assert record.representation_keys == (
        "fit",
        "residuals",
        "normalized_pressure",
    )

    plots = inventory.selected_plots
    assert isinstance(plots, plotting.PlotInventory)
    assert plots.module == "eos"
    assert {item.key for item in plots.representations} == {
        "fit",
        "residuals",
        "normalized_pressure",
    }
    assert plots.property_by_key("pressure").unit == "GPa"
    assert plots.property_by_key("residual").unit == "GPa"
    assert plots.property_by_key("normalized_pressure").symbol_math == "F_E"
    assert plots.context_by_key("record_id").values == (1,)
    assert plots.context_by_key("result_slot").values == ("pv/volume",)


def test_every_advertised_eos_representation_is_buildable(tmp_path: Path) -> None:
    """Inventory keys are accepted directly by the existing public builder."""
    path = _quartz_archive(tmp_path)
    inventory = eos.describe_plots(path)
    assert inventory.selected_plots is not None

    for representation in inventory.selected_plots.representations:
        collection = eos.build_plots(path, (representation.key,))
        assert collection.plots


def test_vt_inventory_matches_builder_availability(tmp_path: Path) -> None:
    """V-T discovery exposes only the fitted curve and available diagnostics."""
    path = _vt_archive(tmp_path)

    inventory = eos.describe_plots(path)

    assert inventory.selected_plots is not None
    keys = tuple(item.key for item in inventory.selected_plots.representations)
    assert keys == tuple(
        item.replace("-", "_") for item in eos.available_plot_types(path)
    )
    assert "fit" in keys
    assert "residuals" in keys
    assert inventory.selected_plots.property_by_key("volume").unit == "dimensionless"
    for key in keys:
        assert eos.build_plots(path, (key,)).plots


def test_pvt_inventory_exposes_calculated_curve_contexts(tmp_path: Path) -> None:
    """P-V-T discovery preserves archive context and calculated-curve semantics."""
    path = _pvt_archive(tmp_path)

    inventory = eos.describe_plots(path)

    assert inventory.selected_plots is not None
    plots = inventory.selected_plots
    keys = {item.key for item in plots.representations}
    assert {"coverage", "isotherms", "isobars", "residuals"} <= keys
    assert plots.context_by_key("isotherm_temperature").unit == "K"
    assert plots.context_by_key("isobar_pressure").unit == "GPa"
    assert "not interpolated" in " ".join(
        plots.representation_by_key("isotherms").constraints
    )
    assert "not interpolated" in " ".join(
        plots.representation_by_key("isobars").constraints
    )
    for representation in plots.representations:
        collection = eos.build_plots(
            path,
            (representation.key,),
            options=eos.PlotOptions(isotherms=(300.0,), isobars=(0.0,)),
        )
        assert collection.plots


def test_eos_inventory_keeps_unaccepted_archive_browsable(tmp_path: Path) -> None:
    """Archive discovery does not require an accepted result."""
    path = _quartz_archive(tmp_path, accept=False)

    inventory = eos.describe_plots(path)

    assert inventory.selected_record_id is None
    assert inventory.selected_plots is None
    assert any("no accepted record" in warning for warning in inventory.warnings)
    record = inventory.record_by_id(1)
    assert record.disposition is eos.RecordDisposition.UNDECIDED
    assert record.plottable
    slot = inventory.slot_by_key("pv/volume")
    assert slot.status is eos.SlotStatus.ATTEMPTED
    assert slot.accepted_record_id is None
    assert slot.plottable_record_ids == (1,)

    selected = eos.describe_plots(path, record_id=1)
    assert selected.selected_record_id == 1
    assert selected.selected_plots is not None


def test_eos_inventory_requires_selection_for_multiple_accepted_slots(
    tmp_path: Path,
) -> None:
    """Multiple accepted slots remain visible without an arbitrary default."""
    path = _topaz_archive(tmp_path)

    inventory = eos.describe_plots(path)

    assert inventory.selected_plots is None
    assert any("multiple accepted records" in warning for warning in inventory.warnings)
    assert {item.key for item in inventory.slots if item.accepted_record_id} == {
        "pv/a",
        "pv/volume",
    }

    selected = eos.describe_plots(path, slot="pv/a")
    assert selected.selected_plots is not None
    assert selected.record_by_id(selected.selected_record_id).slot_key == "pv/a"
    assert selected.selected_plots.context_by_key("fit_target").values == ("a",)
    assert selected.selected_plots.property_by_key("pressure").unit == "GPa"


def test_eos_cli_lists_public_inventory_descriptions(tmp_path: Path) -> None:
    """The CLI no longer owns a second manual EOS plot-description catalog."""
    path = _quartz_archive(tmp_path)

    result = CliRunner().invoke(main, ["eos", "plot", str(path), "--list-plots"])

    assert result.exit_code == 0, result.output
    assert "fit" in result.output
    assert "Observed data and the fitted P-V or V-T relation." in result.output
    assert "normalized_pressure" in result.output


def test_eos_cli_accepts_canonical_inventory_representation_key(
    tmp_path: Path,
) -> None:
    """A key listed by the public inventory can be passed back to the CLI."""
    path = _quartz_archive(tmp_path)
    output = tmp_path / "canonical_key_plot"

    result = CliRunner().invoke(
        main,
        [
            "eos",
            "plot",
            str(path),
            "--plot",
            "normalized_pressure",
            "--output",
            str(output),
        ],
    )

    assert result.exit_code == 0, result.output
    assert any(output.glob("*.png"))


def test_eos_cli_lists_slots_when_archive_selection_is_ambiguous(
    tmp_path: Path,
) -> None:
    """Archive-level CLI discovery reports slots instead of guessing a record."""
    path = _topaz_archive(tmp_path)

    result = CliRunner().invoke(main, ["eos", "plot", str(path), "--list-plots"])

    assert result.exit_code == 0, result.output
    assert "multiple accepted records" in result.output
    assert "Available EOS result slots and plottable records" in result.output
    assert "pv/a" in result.output
    assert "pv/volume" in result.output
