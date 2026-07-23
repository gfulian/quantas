"""Characterization tests for the elasticity command-line boundary."""

from __future__ import annotations

from pathlib import Path
import importlib

import h5py
import numpy as np
from click.testing import CliRunner

from quantas.core.events import Event, EventLevel
from quantas.core.physics.elasticity import (
    DirectionalExtrema,
    ElasticAverages,
    IsotropicElasticProperties,
    StabilityResult,
)
from quantas.models import ResultData, ResultMetadata
from quantas.cli.elasticity import elasticity
from quantas.cli.elasticity_observer import ElasticityTextObserver
from quantas.modules.elasticity.io.export import ElasticityHDF5Export
from quantas.modules.elasticity.models import (
    ElasticityInput,
    ElasticityOptions,
    ElasticityResult,
)

commands_module = importlib.import_module("quantas.cli.elasticity")

DATA = Path(__file__).parent / "data" / "hydroxylapatite.dat"
CRYSTAL_CALCITE = (
    Path(__file__).parents[2]
    / "interfaces"
    / "data"
    / "calcite_crystal_elastcon_excerpt.out"
)


def _properties(value: float) -> IsotropicElasticProperties:
    """Build deterministic isotropic properties for CLI tests."""
    return IsotropicElasticProperties(value, value + 10.0, value / 2.0, 0.25)


def _averages() -> ElasticAverages:
    """Build deterministic Voigt-Reuss-Hill averages."""
    return ElasticAverages(
        voigt=_properties(100.0),
        reuss=_properties(90.0),
        hill=_properties(95.0),
    )


def _variation() -> DirectionalExtrema:
    """Build deterministic directional extrema."""
    return DirectionalExtrema(
        minimum=100.0,
        maximum=200.0,
        anisotropy=2.0,
        minimum_axis=[1.0, 0.0, 0.0],
        maximum_axis=[0.0, 1.0, 0.0],
        minimum_measurement_axis=[0.0, 0.0, 1.0],
        maximum_measurement_axis=[0.0, 0.0, -1.0],
    )


def _elasticity_result() -> ElasticityResult:
    """Build a complete synthetic elasticity result."""
    angle = np.linspace(0.0, 2.0 * np.pi, 5)
    return ElasticityResult(
        jobname="Synthetic",
        crystal_system="orthorhombic",
        stiffness=np.eye(6) * 100.0,
        compliance=np.eye(6) * 0.01,
        averages=_averages(),
        stability=StabilityResult(True, np.arange(1.0, 7.0)),
        variations={
            "young_modulus": _variation(),
            "linear_compressibility": _variation(),
            "shear_modulus": _variation(),
            "poisson_ratio": _variation(),
        },
        properties_2d={
            "xy": {
                "theta": np.full(5, np.pi / 2.0),
                "phi": angle,
                "young_modulus": np.full(5, 100.0),
            }
        },
    )


def test_text_observer_renders_structured_result_events(tmp_path) -> None:
    """The CLI observer converts calculator events into a persistent report."""
    report_file = tmp_path / "elasticity.txt"
    observer = ElasticityTextObserver(report_file=report_file, silent=True)
    input_data = ElasticityInput(jobname="Synthetic")
    options = ElasticityOptions(calculate_2d=True, ntheta=5)
    result = _elasticity_result()

    events = [
        Event(
            "settings",
            EventLevel.RESULT,
            data={"kind": "settings", "options": options},
        ),
        Event(
            "input",
            EventLevel.RESULT,
            data={"kind": "input", "input": input_data, "result": result},
        ),
        Event(
            "averages",
            EventLevel.RESULT,
            data={"kind": "averages", "result": result},
        ),
        Event(
            "stability",
            EventLevel.RESULT,
            data={"kind": "stability", "result": result},
        ),
        Event(
            "variations",
            EventLevel.RESULT,
            data={"kind": "variations", "result": result},
        ),
        Event("example warning", EventLevel.WARNING),
    ]
    for event in events:
        observer(event)
    observer.save()

    text = observer.text()
    assert "Elasticity options" in text
    assert "Elasticity input" in text
    assert "Voigt-Reuss-Hill average properties" in text
    assert "Mechanical stability" in text
    assert "Isotropic seismic velocities" not in text
    assert "Directional elastic extrema" in text
    assert "WARNING: example warning" in text
    assert report_file.read_text(encoding="utf-8") == text


def test_run_command_maps_cli_flags_to_library_options(tmp_path, monkeypatch) -> None:
    """The Click command remains a thin adapter over the library API."""
    captured: dict[str, object] = {}

    def fake_run_elasticity(filename, options=None, observer=None):
        captured["filename"] = filename
        captured["options"] = options
        captured["observer"] = observer
        return ResultData(
            metadata=ResultMetadata(module="elasticity", method="second_order"),
            results={"elasticity": _elasticity_result()},
        )

    def fake_write_elasticity_hdf5(result, filename, report_text=None):
        output = Path(filename)
        output.write_text("synthetic hdf5 placeholder", encoding="utf-8")
        return output

    monkeypatch.setattr(commands_module, "run_elasticity", fake_run_elasticity)
    monkeypatch.setattr(
        commands_module,
        "write_elasticity_hdf5",
        fake_write_elasticity_hdf5,
    )

    outfile = tmp_path / "result.hdf5"
    report = tmp_path / "result.log"
    runner = CliRunner()
    response = runner.invoke(
        elasticity,
        [
            "run",
            str(DATA),
            "--2d",
            "--ntheta",
            "12",
            "--quiet",
            "--output",
            str(outfile),
            "--report",
            str(report),
        ],
    )

    assert response.exit_code == 0, response.output
    options = captured["options"]
    assert isinstance(options, ElasticityOptions)
    assert options.calculate_2d is True
    assert options.ntheta == 12
    assert not hasattr(options, "calculate_isotropic_velocities")
    assert outfile.exists()
    assert report.exists()


def test_export_command_reads_hdf5_and_writes_requested_table(tmp_path) -> None:
    """The CLI export command exercises the HDF5-to-table path end to end."""
    hdf5_file = tmp_path / "elasticity.hdf5"
    output_file = tmp_path / "polar_output.dat"
    ElasticityHDF5Export().export(
        ResultData(
            metadata=ResultMetadata(module="elasticity", method="second_order"),
            results={"elasticity": _elasticity_result()},
        ),
        hdf5_file,
    )

    response = CliRunner().invoke(
        elasticity,
        ["export", str(hdf5_file), "--output", str(output_file)],
    )

    assert response.exit_code == 0, response.output
    assert output_file.exists()
    assert "Elasticity 2D data exported" in response.output
    assert "theta_deg" in output_file.read_text(encoding="utf-8")


def test_inpgen_preserves_crystal_source_components(tmp_path) -> None:
    """The CLI input generator remains a pure external-output parser."""
    outfile = tmp_path / "calcite.dat"
    response = CliRunner().invoke(
        elasticity,
        [
            "inpgen",
            str(CRYSTAL_CALCITE),
            "--output",
            str(outfile),
            "--interface",
            "crystal",
        ],
        input="Calcite source frame\n",
    )

    assert response.exit_code == 0, response.output
    text = outfile.read_text(encoding="utf-8")
    assert "quantas_tensor_frame" not in text


def test_elasticity_help_uses_current_dimension_flags() -> None:
    """The public CLI uses explicit 2D/3D sampling terminology."""
    response = CliRunner().invoke(elasticity, ["run", "--help"])
    assert response.exit_code == 0, response.output
    assert "--2d" in response.output
    assert "--3d" in response.output
    assert "--ntheta" in response.output
    assert "--nphi" in response.output
    assert "--polar" not in response.output
    assert "--no-isotropic-velocities" not in response.output


def test_run_persists_3d_data_without_rendering(tmp_path) -> None:
    """A 3D run stores numerical data even when figures are not requested."""
    outfile = tmp_path / "elasticity.hdf5"
    response = CliRunner().invoke(
        elasticity,
        [
            "run",
            str(DATA),
            "--3d",
            "--ntheta",
            "5",
            "--nphi",
            "7",
            "--property",
            "young",
            "--no-plot",
            "--quiet",
            "--output",
            str(outfile),
            "--report",
            str(tmp_path / "elasticity-3d.log"),
        ],
    )
    assert response.exit_code == 0, response.output
    assert outfile.exists()
    with h5py.File(outfile, "r") as h5:
        assert "properties_3d" in h5["results"]
        assert "young" in h5["results/properties_3d/surfaces"]


def test_plot_command_generates_transient_3d_surface(tmp_path) -> None:
    """The installed-style plot path rebuilds 3D data from HDF5 stiffness."""
    hdf5_file = tmp_path / "elasticity.hdf5"
    ElasticityHDF5Export().export(
        ResultData(
            metadata=ResultMetadata(module="elasticity", method="second_order"),
            results={"elasticity": _elasticity_result()},
        ),
        hdf5_file,
    )
    base = tmp_path / "surface"
    response = CliRunner().invoke(
        elasticity,
        [
            "plot",
            str(hdf5_file),
            "--3d",
            "--ntheta",
            "5",
            "--nphi",
            "7",
            "--property",
            "young",
            "--output",
            str(base),
        ],
    )
    assert response.exit_code == 0, response.output
    assert (tmp_path / "surface_3d_young.png").exists()


def test_run_command_accepts_xyz_tensor_rotation(tmp_path, monkeypatch) -> None:
    """The CLI converts fixed-axis angles into the shared rotation model."""
    captured: dict[str, object] = {}

    def fake_run_elasticity(filename, options=None, observer=None):
        captured["options"] = options
        return ResultData(
            metadata=ResultMetadata(module="elasticity", method="second_order"),
            results={"elasticity": _elasticity_result()},
        )

    def fake_write_elasticity_hdf5(result, filename, report_text=None):
        output = Path(filename)
        output.write_text("synthetic", encoding="utf-8")
        return output

    monkeypatch.setattr(commands_module, "run_elasticity", fake_run_elasticity)
    monkeypatch.setattr(
        commands_module,
        "write_elasticity_hdf5",
        fake_write_elasticity_hdf5,
    )
    response = CliRunner().invoke(
        elasticity,
        [
            "run",
            str(DATA),
            "--rotate-xyz",
            "0",
            "0",
            "30",
            "--quiet",
            "--output",
            str(tmp_path / "rotated.hdf5"),
            "--report",
            str(tmp_path / "rotated.log"),
        ],
    )

    assert response.exit_code == 0, response.output
    options = captured["options"]
    assert isinstance(options, ElasticityOptions)
    assert options.rotation is not None
    assert options.rotation.angles == (0.0, 0.0, 30.0)


def test_run_command_accepts_matrix_file_and_rejects_two_modes(
    tmp_path,
    monkeypatch,
) -> None:
    """Explicit matrices are supported and mutually exclusive with XYZ angles."""
    matrix_file = tmp_path / "rotation.txt"
    matrix_file.write_text("1 0 0\n0 1 0\n0 0 1\n", encoding="utf-8")
    captured: dict[str, object] = {}

    def fake_run_elasticity(filename, options=None, observer=None):
        captured["options"] = options
        return ResultData(
            metadata=ResultMetadata(module="elasticity", method="second_order"),
            results={"elasticity": _elasticity_result()},
        )

    monkeypatch.setattr(commands_module, "run_elasticity", fake_run_elasticity)
    monkeypatch.setattr(
        commands_module,
        "write_elasticity_hdf5",
        lambda result, filename, report_text=None: Path(filename),
    )
    valid = CliRunner().invoke(
        elasticity,
        [
            "run",
            str(DATA),
            "--rotation-matrix",
            str(matrix_file),
            "--quiet",
            "--output",
            str(tmp_path / "matrix.hdf5"),
            "--report",
            str(tmp_path / "matrix.log"),
        ],
    )
    assert valid.exit_code == 0, valid.output
    options = captured["options"]
    assert isinstance(options, ElasticityOptions)
    np.testing.assert_array_equal(options.rotation.matrix, np.eye(3))

    invalid = CliRunner().invoke(
        elasticity,
        [
            "run",
            str(DATA),
            "--rotate-xyz",
            "0",
            "0",
            "30",
            "--rotation-matrix",
            str(matrix_file),
            "--report",
            str(tmp_path / "invalid-rotation.log"),
        ],
    )
    assert invalid.exit_code != 0
    assert "either --rotate-xyz or --rotation-matrix" in invalid.output


def test_invalid_colormap_fails_before_elasticity_calculation(
    tmp_path, monkeypatch
) -> None:
    """A misspelled colormap is rejected before the scientific workflow."""
    called = False

    def fake_run_elasticity(*args, **kwargs):
        nonlocal called
        called = True
        raise AssertionError("calculation should not start")

    monkeypatch.setattr(commands_module, "run_elasticity", fake_run_elasticity)
    response = CliRunner().invoke(
        elasticity,
        [
            "run",
            str(DATA),
            "--3d",
            "--plot",
            "--cmap",
            "virids",
            "--output",
            str(tmp_path / "unused.hdf5"),
            "--report",
            str(tmp_path / "invalid-colormap.log"),
        ],
    )
    assert response.exit_code != 0
    assert called is False
    assert "Unknown Matplotlib colormap 'virids'" in response.output
    assert "viridis" in response.output


def test_run_maps_persisted_3d_options_to_library(tmp_path, monkeypatch) -> None:
    """The run adapter passes the requested 3D grid to the calculator."""
    captured: dict[str, object] = {}

    def fake_run_elasticity(filename, options=None, observer=None):
        captured["options"] = options
        return ResultData(
            metadata=ResultMetadata(module="elasticity", method="second_order"),
            results={"elasticity": _elasticity_result()},
        )

    def fake_write(result, filename, report_text=None):
        path = Path(filename)
        path.write_text("synthetic", encoding="utf-8")
        return path

    monkeypatch.setattr(commands_module, "run_elasticity", fake_run_elasticity)
    monkeypatch.setattr(commands_module, "write_elasticity_hdf5", fake_write)
    response = CliRunner().invoke(
        elasticity,
        [
            "run",
            str(DATA),
            "--3d",
            "--ntheta",
            "9",
            "--nphi",
            "13",
            "--property",
            "young",
            "--quiet",
            "--output",
            str(tmp_path / "result.hdf5"),
            "--report",
            str(tmp_path / "result-3d.log"),
        ],
    )
    assert response.exit_code == 0, response.output
    options = captured["options"]
    assert isinstance(options, ElasticityOptions)
    assert options.calculate_3d is True
    assert options.surface_options is not None
    assert options.surface_options.ntheta == 9
    assert options.surface_options.nphi == 13
    assert options.surface_options.properties == ("young",)
