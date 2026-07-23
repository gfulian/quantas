"""Native EOS HDF5 persistence, immutable history, and accepted state."""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import pytest

from quantas.core.math.fitting import OLSOptions, WLSOptions
from quantas.core.physics.eos import MGDNormalization, PVTEOS, PVTModel
from quantas.modules.eos import (
    EOSArchive,
    EOSDataset,
    EOSFitDomain,
    EOSFitOptions,
    EOSFitRequest,
    EOSFitter,
    EOSResultSlot,
    EOSSlotStatus,
    EOSStateEventType,
    ParameterConstraint,
    infer_result_slots,
    read_eos_input,
)

DATA = Path(__file__).with_name("data")


def _quartz_fit(model: str = "BM3"):
    dataset = read_eos_input(DATA / "PV_quartz.dat")
    request = EOSFitRequest(
        model=model,
        options=EOSFitOptions(solver_options=WLSOptions()),
    )
    return dataset, request, EOSFitter().fit(dataset, request)


def _pvt_fit():
    pressure_parameters = {"K0": 160.0, "KP": 4.2, "KPP": -0.02, "V0": 100.0}
    thermal_parameters = {
        "V0": 100.0,
        "temperature_ref": 300.0,
        "alpha0": 3.0e-5,
        "alpha1": 1.0e-8,
    }
    coupling_parameters = {"dK0_dT": -0.02}
    model = PVTModel("BM3", "linear", "berman:quadratic")
    temperature = np.repeat(np.asarray([300.0, 600.0, 900.0]), 6)
    pressure = np.tile(np.linspace(0.0, 10.0, 6), 3)
    volume = PVTEOS().volume(
        model,
        pressure_parameters,
        thermal_parameters,
        coupling_parameters,
        pressure,
        temperature,
    )
    dataset = EOSDataset(
        jobname="synthetic PVT archive",
        columns={"pressure": pressure, "volume": volume, "temperature": temperature},
        units={"pressure": "GPa", "volume": "angstrom^3", "temperature": "K"},
    )
    constraints = (
        ParameterConstraint.free("K0", 158.0),
        ParameterConstraint.free("KP", 4.0),
        ParameterConstraint.free("V0", 99.8),
        ParameterConstraint.fixed("temperature_ref", 300.0),
        ParameterConstraint.free("alpha0", 2.8e-5),
        ParameterConstraint.free("alpha1", 0.8e-8),
        ParameterConstraint.free("dK0_dT", -0.018),
    )
    request = EOSFitRequest(
        model=model,
        domain=EOSFitDomain.PRESSURE_VOLUME_TEMPERATURE,
        constraints=constraints,
        options=EOSFitOptions(solver_options=OLSOptions(max_iterations=5000)),
    )
    return dataset, request, EOSFitter().fit(dataset, request)


def _mgd_fit():
    """Return one fitted synthetic MGD request for persistence tests."""
    model = PVTModel(
        "BM3",
        "thermal-pressure",
        thermal_pressure_model="mgd",
        mgd_normalization=MGDNormalization.cell(
            formula="NaF", formula_units_per_cell=4
        ),
    )
    temperature = np.repeat(np.asarray([100.0, 295.0, 600.0, 1000.0]), 7)
    pressure = np.tile(np.linspace(0.0, 18.0, 7), 4)
    pressure_parameters = {"K0": 48.0, "KP": 4.5, "KPP": -0.1, "V0": 150.0}
    coupling_parameters = {
        "temperature_ref": 295.0,
        "theta_d0": 459.0,
        "gamma0": 1.547,
        "q": 0.94,
    }
    volume = PVTEOS().volume(
        model,
        pressure_parameters,
        None,
        coupling_parameters,
        pressure,
        temperature,
    )
    dataset = EOSDataset(
        jobname="synthetic MGD archive",
        columns={"pressure": pressure, "volume": volume, "temperature": temperature},
        units={"pressure": "GPa", "volume": "angstrom^3", "temperature": "K"},
    )
    request = EOSFitRequest(
        model=model,
        domain="pvt",
        constraints=(
            ParameterConstraint.free("K0", 47.0),
            ParameterConstraint.free("KP", 4.3),
            ParameterConstraint.free("V0", 149.5),
            ParameterConstraint.fixed("temperature_ref", 295.0),
            ParameterConstraint.free("theta_d0", 440.0, lower_bound=100.0),
            ParameterConstraint.free("gamma0", 1.4),
            ParameterConstraint.free("q", 0.8),
        ),
        options=EOSFitOptions(solver_options=OLSOptions(max_iterations=10000)),
    )
    return dataset, request, EOSFitter().fit(dataset, request)


def test_slot_contract_and_inference() -> None:
    topaz = read_eos_input(DATA / "PV_topaz.dat")
    assert EOSResultSlot.parse("pv/c").key == "pv/c"
    assert [slot.key for slot in infer_result_slots(topaz)] == [
        "pv/a",
        "pv/b",
        "pv/c",
        "pv/volume",
    ]
    rutile = read_eos_input(DATA / "rutile.dat")
    assert [slot.key for slot in infer_result_slots(rutile)] == [
        "vt/a",
        "vt/c",
        "vt/volume",
    ]


def test_archive_can_store_input_without_any_fit(tmp_path: Path) -> None:
    path = tmp_path / "topaz_empty.hdf5"
    dataset = read_eos_input(DATA / "PV_topaz.dat")
    with EOSArchive.create(path, dataset=dataset) as archive:
        states = {state.slot.key: state for state in archive.slots()}
        assert states["pv/volume"].status is EOSSlotStatus.NOT_PROCESSED
        assert states["pv/c"].accepted_record_id is None
        assert archive.record_ids == ()
        archive.record_no_operation("pv/c", note="axis not requested")
        with pytest.raises(KeyError, match="registered result slot"):
            archive.slot_state("vt/c")
        assert archive.slot_state("pv/c").status is EOSSlotStatus.NOT_PROCESSED
    with EOSArchive(path) as archive:
        assert archive.record_ids == ()
        assert archive.events()[-1].event_type is EOSStateEventType.NO_OPERATION


def test_dataset_round_trip_is_self_contained_and_float64(tmp_path: Path) -> None:
    path = tmp_path / "quartz.hdf5"
    dataset = read_eos_input(DATA / "PV_quartz.dat")
    with EOSArchive.create(path, dataset=dataset) as archive:
        restored = archive.dataset(1)
        assert restored.jobname == dataset.jobname
        assert restored.source == dataset.source
        assert restored.metadata == dataset.metadata
        for name in dataset.column_names:
            assert restored.column(name).dtype == np.float64
            np.testing.assert_array_equal(restored.column(name), dataset.column(name))
    with h5py.File(path, "r") as h5:
        node = h5["input/datasets/000001"]
        assert "source_text" in node
        assert len(str(node.attrs["source_sha256"])) == 64
        assert h5["input/datasets/000001/normalized_data/volume"].dtype == np.float64
        assert h5["metadata/numerics"].attrs["storage_dtype"] == "float64"


def test_single_script_fit_creates_record_and_accepted_materialization(
    tmp_path: Path,
) -> None:
    dataset, request, result = _quartz_fit()
    path = tmp_path / "quartz_fit.hdf5"
    with EOSArchive.create(path, dataset=dataset) as archive:
        record = archive.store_fit(1, request, result, provenance={"frontend": "api"})
        state = archive.slot_state("pv/volume")
        assert record.record_id == 1
        assert state.status is EOSSlotStatus.ACCEPTED
        assert state.accepted_record_id == 1
        assert archive.accepted_result("pv/volume") is not None
        assert archive.accepted_result("pv/volume").parameter_values == pytest.approx(
            result.parameter_values
        )
    with h5py.File(path, "r") as h5:
        assert bool(h5["session/records/000001"].attrs["immutable"])
        accepted = h5["results/accepted/pv/volume"]
        assert int(accepted.attrs["record_id"]) == 1
        assert accepted.attrs["source_record_path"] == "/session/records/000001"


def test_rejected_last_attempt_does_not_replace_accepted_result(tmp_path: Path) -> None:
    dataset, request1, result1 = _quartz_fit("BM3")
    _, request2, result2 = _quartz_fit("BM2")
    path = tmp_path / "history.hdf5"
    with EOSArchive.create(path, dataset=dataset) as archive:
        first = archive.store_fit(1, request1, result1)
        second = archive.append_fit(
            1,
            request2,
            result2,
            parent_record_id=first.record_id,
        )
        archive.reject(second.record_id, note="higher residual structure")
        state = archive.slot_state("pv/volume")
        assert state.accepted_record_id == first.record_id
        assert state.last_record_id == second.record_id
        assert state.attempted_record_ids == (1, 2)
        assert state.status is EOSSlotStatus.ACCEPTED
        assert archive.accepted("pv/volume").record_id == first.record_id
        assert tuple(record.record_id for record in archive.records("pv/volume")) == (
            1,
            2,
        )


def test_parent_lineage_candidate_notes_and_initial_guess_events(
    tmp_path: Path,
) -> None:
    dataset, request1, result1 = _quartz_fit("BM2")
    _, request2, result2 = _quartz_fit("BM3")
    path = tmp_path / "lineage.hdf5"
    with EOSArchive.create(path, dataset=dataset) as archive:
        first = archive.store_fit(1, request1, result1)
        second = archive.append_fit(
            1,
            request2,
            result2,
            parent_record_id=first.record_id,
        )
        archive.mark_candidate(second.record_id)
        archive.mark_used_as_initial_guess(
            first.record_id,
            child_record_id=second.record_id,
        )
        archive.add_note(second.record_id, "Kp released from BM2 value")
        assert archive.record(second.record_id).parent_record_id == first.record_id
        assert [event.event_type for event in archive.events()][-3:] == [
            EOSStateEventType.RECORD_CANDIDATE,
            EOSStateEventType.RECORD_USED_AS_INITIAL_GUESS,
            EOSStateEventType.NOTE_ADDED,
        ]


def test_failed_record_is_persistent_but_cannot_be_accepted(tmp_path: Path) -> None:
    dataset, request, result = _quartz_fit()
    result.fit.success = False
    path = tmp_path / "failed.hdf5"
    with EOSArchive.create(path, dataset=dataset) as archive:
        record = archive.append_fit(1, request, result)
        assert not archive.record(record.record_id).successful
        with pytest.raises(ValueError, match="successful"):
            archive.accept(record.record_id)
        assert archive.slot_state("pv/volume").status is EOSSlotStatus.ATTEMPTED


def test_vt_and_pvt_request_result_round_trip(tmp_path: Path) -> None:
    rutile = read_eos_input(DATA / "rutile.dat")
    vt_request = EOSFitRequest(
        model="berman:quadratic",
        domain="vt",
        options=EOSFitOptions(solver_options=WLSOptions(max_iterations=3000)),
    )
    vt_result = EOSFitter().fit(rutile, vt_request)
    pvt_dataset, pvt_request, pvt_result = _pvt_fit()
    path = tmp_path / "multi_domain.hdf5"
    with EOSArchive.create(path, dataset=rutile) as archive:
        pvt_dataset_id = archive.add_dataset(pvt_dataset)
        vt_record = archive.store_fit(1, vt_request, vt_result)
        pvt_record = archive.store_fit(pvt_dataset_id, pvt_request, pvt_result)
    with EOSArchive(path) as archive:
        restored_vt = archive.record(vt_record.record_id)
        restored_pvt = archive.record(pvt_record.record_id)
        assert restored_vt.request.as_dict() == vt_request.as_dict()
        assert restored_vt.result.parameter_values == pytest.approx(
            vt_result.parameter_values
        )
        assert restored_pvt.request.as_dict() == pvt_request.as_dict()
        assert restored_pvt.result.parameter_values == pytest.approx(
            pvt_result.parameter_values
        )
        assert archive.accepted("vt/volume").record_id == vt_record.record_id
        assert archive.accepted("pvt/volume").record_id == pvt_record.record_id


def test_mgd_request_result_round_trip_preserves_normalization(tmp_path: Path) -> None:
    """Persist and restore an executable MGD fit without losing provenance."""
    dataset, request, result = _mgd_fit()
    assert result.fit.success, result.fit.message
    path = tmp_path / "mgd.hdf5"
    with EOSArchive.create(path, dataset=dataset) as archive:
        record = archive.store_fit(1, request, result)
    with EOSArchive(path) as archive:
        restored = archive.record(record.record_id)
        assert restored.request.as_dict() == request.as_dict()
        assert restored.result.parameter_values == pytest.approx(
            result.parameter_values
        )
        model = restored.request.model
        assert isinstance(model, PVTModel)
        assert model.thermal_pressure_spec is not None
        assert model.thermal_pressure_spec.tag == "mie-gruneisen-debye:full"
        assert model.mgd_normalization_spec is not None
        assert model.mgd_normalization_spec.atoms_per_cell == pytest.approx(8.0)
        assert model.mgd_normalization_spec.formula == "NaF"
        assert archive.accepted("pvt/volume").record_id == record.record_id


def test_read_only_archive_rejects_mutation(tmp_path: Path) -> None:
    dataset = read_eos_input(DATA / "PV_quartz.dat")
    path = tmp_path / "readonly.hdf5"
    with EOSArchive.create(path, dataset=dataset):
        pass
    with EOSArchive(path) as archive:
        with pytest.raises(OSError, match="read-only"):
            archive.register_slot("pv/a")
