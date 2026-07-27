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
