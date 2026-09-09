from __future__ import annotations


from click.testing import CliRunner

from quantas.cli.ha import ha


def test_ha_group_exposes_expected_commands():
    assert set(ha.commands) >= {"run", "export", "plot"}


def test_ha_run_help_is_available():
    runner = CliRunner()
    result = runner.invoke(ha, ["run", "--help"])

    assert result.exit_code == 0
    assert (
        "Harmonic-approximation" in result.output
        or "harmonic-approximation" in result.output
    )
    assert "--eunit" in result.output
    assert "--funit" in result.output
    assert "--kieffer" in result.output


def test_ha_run_activates_embedded_kieffer_cutoffs(tmp_path, monkeypatch) -> None:
    """The HA execution flag reads and forwards the embedded cutoff state."""
    filename = tmp_path / "ha-kieffer.yaml"
    filename.write_text("job: Kieffer CLI\n", encoding="utf-8")
    destination = tmp_path / "result.hdf5"
    cutoffs = object()
    captured = {}

    monkeypatch.setattr(
        "quantas.cli.ha.read_ha_kieffer_input",
        lambda path: cutoffs,
    )

    def fake_run(input_data, *, options, kieffer_cutoffs, observer):
        captured["input"] = input_data
        captured["options"] = options
        captured["cutoffs"] = kieffer_cutoffs
        captured["observer"] = observer
        return object()

    monkeypatch.setattr("quantas.cli.ha.run_ha", fake_run)
    monkeypatch.setattr(
        "quantas.cli.ha.write_ha_hdf5",
        lambda *args, **kwargs: destination,
    )

    result = CliRunner().invoke(
        ha,
        [
            "run",
            str(filename),
            "--kieffer",
            "--output",
            str(destination),
            "--quiet",
            "--no-progress",
        ],
    )

    assert result.exit_code == 0, result.output
    assert captured["input"] == filename
    assert captured["cutoffs"] is cutoffs
    report = filename.with_suffix(".log").read_text(encoding="utf-8")
    assert "Susan Werner Kieffer" in report


def test_ha_export_help_is_available():
    runner = CliRunner()
    result = runner.invoke(ha, ["export", "--help"])

    assert result.exit_code == 0
    assert "--property" in result.output


def test_ha_plot_help_is_available():
    runner = CliRunner()
    result = runner.invoke(ha, ["plot", "--help"])

    assert result.exit_code == 0
    assert "--dpi" in result.output
    assert "--axis" in result.output
    assert "--volume" in result.output
    assert "--temperature" in result.output
    assert "--2d" in result.output


def test_ha_inpgen_help_exposes_formula_units_option():
    runner = CliRunner()
    result = runner.invoke(ha, ["inpgen", "--help"])

    assert result.exit_code == 0
    assert "-Z" in result.output
    assert "--formula-units" in result.output


def test_ha_and_qha_register_the_same_phonon_input_command():
    """HA and QHA must reuse one command object, not duplicate inpgen logic."""
    from quantas.cli.qha import qha

    assert ha.commands["inpgen"] is qha.commands["inpgen"]


def test_ha_and_qha_register_the_same_kieffer_input_command():
    """HA and QHA share one frontend adapter for Kieffer enrichment."""
    from quantas.cli.qha import qha

    assert ha.commands["add-kieffer"] is qha.commands["add-kieffer"]


def test_add_kieffer_help_exposes_pressure_and_quadrature_controls():
    """The enrichment CLI documents its scientific and numerical choices."""
    result = CliRunner().invoke(ha, ["add-kieffer", "--help"])

    assert result.exit_code == 0
    assert "--interface" in result.output
    assert "--elastic-list" in result.output
    assert "--pressure-source" in result.output
    assert "energy-eos" in result.output
    assert "energy-polynomial" in result.output
    assert "--pressure" in result.output
    assert "--eos" in result.output
    assert "--degree" in result.output
    assert "--mu-order" in result.output
    assert "--phi-order" in result.output


def test_qha_add_kieffer_forwards_interface_and_energy_fit_options(
    tmp_path, monkeypatch
) -> None:
    """The shared command preserves established interface/EOS option names."""
    from quantas.cli.qha import qha

    source = tmp_path / "qha.yaml"
    source.write_text("job: test\n", encoding="utf-8")
    elastic = tmp_path / "elastic.out"
    elastic.write_text("test\n", encoding="utf-8")
    destination = tmp_path / "qha-kieffer.yaml"
    captured = {}

    def fake_add(source, destination, outputs, **kwargs):
        captured.update(kwargs)
        return destination

    monkeypatch.setattr(
        "quantas.cli.kieffer_input.qha_api.add_kieffer_input",
        fake_add,
    )
    result = CliRunner().invoke(
        qha,
        [
            "add-kieffer",
            str(source),
            str(elastic),
            "--interface",
            "crystal",
            "--pressure-source",
            "energy-polynomial",
            "--degree",
            "4",
            "--eos",
            "V3",
            "--output",
            str(destination),
        ],
    )

    assert result.exit_code == 0, result.output
    assert captured["interface"] == "crystal"
    assert captured["pressure_policy"] == "energy_polynomial"
    assert captured["polynomial_degree"] == 4
    assert captured["eos"] == "V3"


def test_generated_structure_summary_does_not_require_new_reader_property(
    tmp_path, monkeypatch
):
    """The CLI read-back must work during staged upgrades of Quantas files."""
    from quantas.cli.phonon_input import _read_generated_structure_summary

    filename = tmp_path / "generated.yaml"
    filename.write_text("job: staged upgrade\n", encoding="utf-8")

    class CompatibleReader:
        completed = True
        error = None
        data = {
            "volume": [10.0],
            "energy": [-1.0],
            "phonon": [],
            "structure": {
                "atomic_numbers": [12, 8],
                "volume_series": {"lattice": [[[1.0, 0.0, 0.0]] * 3]},
                "symmetry": {"international_symbol": "Fm-3m"},
            },
        }

        def __init__(self, source):
            self.source = source

    monkeypatch.setattr(
        "quantas.cli.phonon_input.PhononInputFileReader",
        CompatibleReader,
    )

    assert _read_generated_structure_summary(filename) == (2, 1, "Fm-3m")


def test_ha_plot_forwards_exact_grid_sections_to_public_api(
    tmp_path, monkeypatch
) -> None:
    from types import SimpleNamespace

    filename = tmp_path / "ha.h5"
    filename.write_bytes(b"test")
    captured = {}

    monkeypatch.setattr("quantas.cli.ha.read_ha_hdf5", lambda path: object())

    def fake_build(result, *, properties, unit, options):
        captured["properties"] = properties
        captured["unit"] = unit
        captured["options"] = options
        return object()

    monkeypatch.setattr("quantas.cli.ha.build_ha_plots", fake_build)
    monkeypatch.setattr(
        "quantas.cli.ha.render_plot_collection",
        lambda collection, options: SimpleNamespace(artifacts=[]),
    )

    result = CliRunner().invoke(
        ha,
        [
            "plot",
            str(filename),
            "--property",
            "F",
            "--axis",
            "volume",
            "--temperature",
            "100",
            "--temperature",
            "200",
            "--2d",
        ],
    )

    assert result.exit_code == 0
    assert captured["properties"] == "F"
    assert captured["options"].curve_axis == "volume"
    assert captured["options"].selected_temperatures == (100.0, 200.0)
    assert captured["options"].include_contours is True


def test_ha_inpgen_help_exposes_quiet_and_debug_modes():
    result = CliRunner().invoke(ha, ["inpgen", "--help"])

    assert result.exit_code == 0
    assert "--quiet" in result.output
    assert "--debug" in result.output


def test_inpgen_quiet_suppresses_success_output(tmp_path, monkeypatch):
    source = tmp_path / "phonons.out"
    source.write_text("fake", encoding="utf-8")
    destination = tmp_path / "phonons.yaml"

    def fake_create_input(source, destination, **kwargs):
        assert kwargs["observer"] is not None
        return destination

    monkeypatch.setattr(
        "quantas.cli.phonon_input.create_ha_input",
        fake_create_input,
    )
    monkeypatch.setattr(
        "quantas.cli.phonon_input._read_generated_structure_summary",
        lambda filename: None,
    )

    result = CliRunner().invoke(
        ha,
        [
            "inpgen",
            str(source),
            "--output",
            str(destination),
            "--jobname",
            "quiet test",
            "--quiet",
        ],
    )

    assert result.exit_code == 0
    assert result.output == ""


def test_inpgen_debug_renders_mode_tracking_table(tmp_path, monkeypatch):
    from quantas.core.events import Event, EventLevel
    from quantas.models import ReportTable

    source = tmp_path / "phonons.out"
    source.write_text("fake", encoding="utf-8")
    destination = tmp_path / "phonons.yaml"

    summary = ReportTable(
        title="Phonon mode continuity",
        columns=["Volume from", "Volume to", "Unresolved"],
        rows=[[10.0, 11.0, 0]],
        metadata={
            "column_units": ["Å³", "Å³", ""],
            "column_formats": [".4f", ".4f", "integer"],
        },
    )
    detail = ReportTable(
        title="Mode tracking: q #1",
        columns=["Branch", "Frequency from", "Frequency to", "Status"],
        rows=[[1, 100.0, 101.0, "matched"]],
        metadata={
            "column_units": ["", "cm^-1", "cm^-1", ""],
            "column_formats": ["integer", ".4f", ".4f", None],
        },
    )

    def fake_create_input(source, destination, **kwargs):
        observer = kwargs["observer"]
        observer(
            Event(
                "mode tracking complete",
                level=EventLevel.RESULT,
                data={"kind": "mode_tracking_summary", "table": summary},
            )
        )
        observer(
            Event(
                "mode tracking detail",
                level=EventLevel.DEBUG,
                data={"kind": "mode_tracking_detail", "table": detail},
            )
        )
        return destination

    monkeypatch.setattr(
        "quantas.cli.phonon_input.create_ha_input",
        fake_create_input,
    )
    monkeypatch.setattr(
        "quantas.cli.phonon_input._read_generated_structure_summary",
        lambda filename: None,
    )

    result = CliRunner().invoke(
        ha,
        [
            "inpgen",
            str(source),
            "--output",
            str(destination),
            "--jobname",
            "debug test",
            "--debug",
        ],
        terminal_width=240,
    )

    assert result.exit_code == 0
    assert "Phonon mode continuity" in result.output
    assert "Mode tracking: q #1" in result.output
    assert "(Å³)" in result.output
    assert "(cm^-1)" in result.output


def test_inpgen_rejects_quiet_and_debug_together(tmp_path):
    source = tmp_path / "phonons.out"
    source.write_text("fake", encoding="utf-8")

    result = CliRunner().invoke(
        ha,
        [
            "inpgen",
            str(source),
            "--jobname",
            "conflict test",
            "--quiet",
            "--debug",
        ],
    )

    assert result.exit_code != 0
    assert "--quiet and --debug cannot be used together" in result.output
