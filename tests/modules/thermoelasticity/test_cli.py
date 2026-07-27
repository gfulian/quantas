"""CLI integration tests for the thermoelastic workflow."""

from __future__ import annotations

import csv
from pathlib import Path

from click.testing import CliRunner
import pytest

from quantas.cli.main import main
from quantas.modules.elasticity.io.reader import ElasticityInputFileReader
from quantas.modules.seismic.io.reader import SeismicInputFileReader
from quantas.modules.thermoelasticity import (
    read_thermoelastic_hdf5,
    read_thermoelastic_profile_spec,
)

pytestmark = pytest.mark.cli

_PROJECT_ROOT = Path(__file__).resolve().parents[3]
_THERMOELASTIC_EXAMPLES = _PROJECT_ROOT / "examples" / "thermoelasticity"
_QHA_EXAMPLES = _PROJECT_ROOT / "examples" / "qha" / "crystal-phonons"


def _fit_archive(
    runner: CliRunner,
    tmp_path: Path,
    *,
    plot: bool = False,
) -> Path:
    archive = tmp_path / "dolomite_fit.hdf5"
    command = [
        "thermoelasticity",
        "run",
        str(_THERMOELASTIC_EXAMPLES / "dol_pbe0_thermoelastic.yaml"),
        str(_QHA_EXAMPLES / "dol_pbe0.yaml"),
        "-o",
        str(archive),
        "-r",
        str(tmp_path / "dolomite_fit.log"),
        "--qha-temperature",
        "300",
        "500",
        "100",
        "--qha-pressure",
        "0",
        "5",
        "5",
        "--quiet",
        "--no-progress",
        "--adiabatic",
        "require",
    ]
    if plot:
        command.extend(["--plot", "--plot-output", str(tmp_path / "fit_plots")])
    result = runner.invoke(main, command)
    assert result.exit_code == 0, result.output
    return archive


def test_thermoelasticity_group_exposes_simplified_stages() -> None:
    """The unpublished CLI uses analysis subcommands rather than pt-analysis."""
    result = CliRunner().invoke(main, ["thermoelasticity", "--help"])
    assert result.exit_code == 0, result.output
    for name in ("inpgen", "run", "analysis", "inspect", "table", "export", "plot"):
        assert name in result.output
    assert "pt-analysis" not in result.output

    analysis = CliRunner().invoke(main, ["thermoelasticity", "analysis", "--help"])
    assert analysis.exit_code == 0, analysis.output
    for name in ("point", "grid", "profile"):
        assert name in analysis.output


def test_run_help_prioritizes_presets_and_scientific_choices() -> None:
    """Advanced fit thresholds remain callable but are hidden from ordinary help."""
    result = CliRunner().invoke(main, ["thermoelasticity", "run", "--help"])
    assert result.exit_code == 0, result.output
    for heading in (
        "Scientific model:",
        "Calculation domain:",
        "Validation and diagnostics:",
        "Plotting:",
        "Output and reporting:",
    ):
        assert heading in result.output
    assert "--validation [standard|strict|exploratory]" in result.output
    assert "--plot-preset [screen|publication|monochrome]" in result.output
    assert "--maximum-design-condition" not in result.output


def test_run_can_render_observed_fit_plots(tmp_path: Path) -> None:
    """Model calibration optionally creates its natural fit diagnostics."""
    runner = CliRunner()
    archive = _fit_archive(runner, tmp_path, plot=True)
    assert archive.exists()
    assert (tmp_path / "fit_plots" / "C11_cold_finite_strain_fit.png").exists()


def test_analysis_point_prints_and_writes_shared_state_input(tmp_path: Path) -> None:
    """A point can be passed directly to Elasticity or SEISMIC without HDF5."""
    runner = CliRunner()
    fit = _fit_archive(runner, tmp_path)
    state = tmp_path / "dolomite_P2.5_T400.dat"
    result = runner.invoke(
        main,
        [
            "thermoelasticity",
            "analysis",
            "point",
            str(fit),
            "2.5",
            "400",
            "-o",
            str(state),
        ],
    )
    assert result.exit_code == 0, result.output
    assert "Thermoelastic point (adiabatic)" in result.output
    assert "C11_S" in result.output
    assert "+/-" in result.output
    assert state.exists()

    seismic = SeismicInputFileReader(state)
    elasticity = ElasticityInputFileReader(state)
    assert seismic.completed
    assert elasticity.completed
    assert seismic.density > 0.0
    assert seismic.jobname.startswith("Dolomite PBE0 QSA thermoelasticity")
    assert seismic.stiffness is not None
    assert seismic.stiffness == pytest.approx(elasticity.stiffness)

    numeric_rows = [
        line.split()
        for line in state.read_text(encoding="utf-8").splitlines()
        if len(line.split()) == 6
        and all(_is_float_token(token) for token in line.split())
    ]
    assert len(numeric_rows) == 6


def test_analysis_grid_and_wide_table_round_trip(tmp_path: Path) -> None:
    """Grid analysis stores HDF5 while the wide table avoids duplicated rows."""
    runner = CliRunner()
    fit = _fit_archive(runner, tmp_path)
    archive = tmp_path / "dolomite_grid.hdf5"
    table = tmp_path / "dolomite_grid.csv"
    result = runner.invoke(
        main,
        [
            "thermoelasticity",
            "analysis",
            "grid",
            str(fit),
            "-P",
            "0",
            "5",
            "5",
            "-T",
            "300",
            "500",
            "100",
            "-o",
            str(archive),
            "--table-output",
            str(table),
            "--tensor-condition",
            "both",
            "--force",
        ],
    )
    assert result.exit_code == 0, result.output
    payload = read_thermoelastic_hdf5(archive).results["thermoelasticity"]
    assert payload.stiffness_isothermal is not None
    assert payload.stiffness_isothermal.shape == (3, 2, 6, 6)
    assert payload.stiffness_adiabatic is not None
    with table.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 6
    assert "C11_T_GPa" in rows[0]
    assert "C11_S_GPa" in rows[0]
    assert "C22_T_GPa" not in rows[0]

    exported = tmp_path / "grid_text.dat"
    table_result = runner.invoke(
        main,
        [
            "thermoelasticity",
            "table",
            "grid",
            str(archive),
            "-o",
            str(exported),
            "--format",
            "text",
            "--tensor-condition",
            "isothermal",
            "--force",
        ],
    )
    assert table_result.exit_code == 0, table_result.output
    header = exported.read_text(encoding="utf-8").splitlines()[0]
    assert "C11_T_GPa" in header
    assert "C11_S_GPa" not in header


def test_analysis_profile_is_profile_only_and_writes_wide_table(tmp_path: Path) -> None:
    """Profile analysis does not create a hidden rectangular P-T stiffness grid."""
    runner = CliRunner()
    fit = _fit_archive(runner, tmp_path)
    archive = tmp_path / "crust.hdf5"
    table = tmp_path / "crust.csv"
    result = runner.invoke(
        main,
        [
            "thermoelasticity",
            "analysis",
            "profile",
            str(fit),
            "--linear-profile",
            "crust",
            "--depth",
            "0",
            "20",
            "10",
            "-o",
            str(archive),
            "--table-output",
            str(table),
            "--plot",
            "--force",
        ],
    )
    assert result.exit_code == 0, result.output
    payload = read_thermoelastic_hdf5(archive).results["thermoelasticity"]
    assert payload.stiffness_isothermal is None
    assert payload.metadata["analysis_stage"] == "profiles_only"
    assert set(payload.profiles) == {"crust"}
    with table.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 3
    assert [float(row["depth_km"]) for row in rows] == [0.0, 10.0, 20.0]
    assert "C14_T_GPa" in rows[0]
    assert "C14_S_GPa" in rows[0]
    plot_directory = tmp_path / "crust_plots"
    assert plot_directory.is_dir()
    assert any(plot_directory.glob("*.png"))


def test_inspect_identifies_archive_stage_and_next_commands(tmp_path: Path) -> None:
    """Inspect explains whether a file is a calibration, grid, or profile result."""
    runner = CliRunner()
    fit = _fit_archive(runner, tmp_path)
    fit_inspect = runner.invoke(main, ["thermoelasticity", "inspect", str(fit)])
    assert fit_inspect.exit_code == 0, fit_inspect.output
    assert "model calibration" in fit_inspect.output
    assert "Frame status" in fit_inspect.output
    assert "normalized" in fit_inspect.output
    assert "analysis profile" in fit_inspect.output

    profile = tmp_path / "profile.hdf5"
    evaluated = runner.invoke(
        main,
        [
            "thermoelasticity",
            "analysis",
            "profile",
            str(fit),
            "--linear-profile",
            "crust",
            "--depth",
            "0",
            "10",
            "5",
            "-o",
            str(profile),
            "--force",
        ],
    )
    assert evaluated.exit_code == 0, evaluated.output
    profile_inspect = runner.invoke(main, ["thermoelasticity", "inspect", str(profile)])
    assert profile_inspect.exit_code == 0, profile_inspect.output
    assert "geological profile analysis" in profile_inspect.output
    assert "crust" in profile_inspect.output


def test_plot_cli_exposes_compact_help_and_compare_family(tmp_path: Path) -> None:
    """Plot CLI keeps scientific choices visible and typography out of help."""
    runner = CliRunner()
    listed = runner.invoke(main, ["thermoelasticity", "plot", "--list-plots"])
    assert listed.exit_code == 0, listed.output
    for name in ("fit", "pt", "profile", "compare", "domain"):
        assert name in listed.output

    profile_help = runner.invoke(main, ["thermoelasticity", "plot", "profile", "--help"])
    assert profile_help.exit_code == 0, profile_help.output
    assert "--tensor-condition" in profile_help.output
    assert "--uncertainty" in profile_help.output
    assert "--preset [screen|publication|monochrome]" in profile_help.output
    assert "--marker-edge-width" not in profile_help.output
    assert "--figure-width" not in profile_help.output

    fit = _fit_archive(runner, tmp_path)
    grid = tmp_path / "grid.hdf5"
    analyzed = runner.invoke(
        main,
        [
            "thermoelasticity",
            "analysis",
            "grid",
            str(fit),
            "-P",
            "0",
            "5",
            "5",
            "-T",
            "300",
            "500",
            "100",
            "-o",
            str(grid),
            "--force",
        ],
    )
    assert analyzed.exit_code == 0, analyzed.output
    available = runner.invoke(
        main,
        [
            "thermoelasticity",
            "plot",
            "--list-plots",
            "--archive",
            str(grid),
        ],
    )
    assert available.exit_code == 0, available.output
    assert "Available thermoelastic plot families" in available.output
    for name in ("fit", "pt", "compare", "domain"):
        assert name in available.output
    assert "profile" not in available.output

    invalid_archive_option = runner.invoke(
        main,
        ["thermoelasticity", "plot", "--archive", str(grid)],
    )
    assert invalid_archive_option.exit_code == 2
    assert "valid only with --list-plots" in invalid_archive_option.output

    output = tmp_path / "compare"
    compare = runner.invoke(
        main,
        [
            "thermoelasticity",
            "plot",
            "compare",
            str(grid),
            "--component",
            "C11",
            "--pressure",
            "0",
            "-o",
            str(output),
            "--dpi",
            "72",
        ],
    )
    assert compare.exit_code == 0, compare.output
    assert (output / "C11_T_S_compare_P0GPa.png").exists()


def test_profile_presets_and_composed_specs_use_analysis_profile(
    tmp_path: Path,
) -> None:
    """Scientific presets and independent P(z)/T(z) specs share one command."""
    runner = CliRunner()
    fit = _fit_archive(runner, tmp_path)
    listed = runner.invoke(
        main,
        ["thermoelasticity", "analysis", "profile", str(fit), "--list-presets"],
    )
    assert listed.exit_code == 0, listed.output
    assert "mantle-katsura-2022" in listed.output
    assert "thin-crust" not in listed.output

    temperature = tmp_path / "temperature.dat"
    temperature.write_text(
        "depth_km T_K\n0 300\n35 700\n100 1200\n",
        encoding="utf-8",
    )
    specification = tmp_path / "earth_profile.yaml"
    specification.write_text(
        """schema_version: 1
name: prem-custom-cli
depth:
  min_km: 0
  max_km: 100
  step_km: 10
pressure:
  model: prem
temperature:
  source: table
  file: temperature.dat
  interpolation: pchip
  citation: User-defined test geotherm.
""",
        encoding="utf-8",
    )
    archive = tmp_path / "profile_spec.hdf5"
    result = runner.invoke(
        main,
        [
            "thermoelasticity",
            "analysis",
            "profile",
            str(fit),
            "--profile-spec",
            str(specification),
            "--extrapolation",
            "allow",
            "-o",
            str(archive),
            "--force",
        ],
    )
    assert result.exit_code == 0, result.output
    profile = (
        read_thermoelastic_hdf5(archive)
        .results["thermoelasticity"]
        .profiles["prem-custom-cli"]
    )
    assert profile.depth[-1] == 100.0
    assert profile.temperature[-1] == pytest.approx(1200.0)
    assert profile.metadata["pressure"]["model"] == "prem"


def test_existing_output_prompts_instead_of_only_requiring_force(
    tmp_path: Path,
) -> None:
    """Interactive writes offer the same overwrite behavior as other modules."""
    runner = CliRunner()
    fit = _fit_archive(runner, tmp_path)
    state = tmp_path / "state.dat"
    state.write_text("existing\n", encoding="utf-8")
    declined = runner.invoke(
        main,
        [
            "thermoelasticity",
            "analysis",
            "point",
            str(fit),
            "0",
            "300",
            "-o",
            str(state),
        ],
        input="n\n",
    )
    assert declined.exit_code == 0, declined.output
    assert "already exists. Replace it?" in declined.output
    assert state.read_text(encoding="utf-8") == "existing\n"


def test_full_tensor_export_still_requires_explicit_scope(tmp_path: Path) -> None:
    """The advanced 21-component export remains deliberate and explicit."""
    runner = CliRunner()
    fit = _fit_archive(runner, tmp_path)
    result = runner.invoke(main, ["thermoelasticity", "export", str(fit)])
    assert result.exit_code != 0
    assert "choose one export selection" in result.output


def test_profile_template_command_writes_parseable_spec(tmp_path: Path) -> None:
    """The CLI writes an editable profile specification accepted by the reader."""
    runner = CliRunner()
    template = tmp_path / "profile.yaml"
    result = runner.invoke(
        main,
        ["thermoelasticity", "profile-template", str(template), "--force"],
    )
    assert result.exit_code == 0, result.output
    profile = read_thermoelastic_profile_spec(template)
    assert profile.name == "custom-earth-profile"
    assert profile.depth[0] == pytest.approx(0.0)
    assert profile.depth[-1] == pytest.approx(660.0)


def _is_float_token(value: str) -> bool:
    try:
        float(value)
    except ValueError:
        return False
    return True
