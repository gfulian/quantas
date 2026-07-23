"""Contracts for the public CLI, reporting, and citations."""

from __future__ import annotations

from pathlib import Path
import shutil

import click
from click.testing import CliRunner

from quantas.cli.contracts import (
    OUTPUT_GROUP,
    ReportVerbosity,
    default_report_path,
)
from quantas.cli.main import main
from quantas.references import (
    CITATIONS,
    METHOD_CITATION_KEYS,
    MODULE_CITATION_KEYS,
    get_citation,
)


_REMOVED_OPTIONS = {
    "--debug",
    "--edeg",
    "--extended",
    "--fdeg",
    "--outbase",
    "--outfile",
}


def _leaf_commands(
    command: click.Command,
) -> list[tuple[tuple[str, ...], click.Command]]:
    """Return all leaf commands below one Click command tree."""
    leaves: list[tuple[tuple[str, ...], click.Command]] = []

    def visit(current: click.Command, path: tuple[str, ...]) -> None:
        if isinstance(current, click.Group):
            for name, child in current.commands.items():
                visit(child, (*path, name))
            return
        leaves.append((path, current))

    visit(command, ())
    return leaves


def _option_names(command: click.Command) -> set[str]:
    """Return primary and secondary option spellings for one command."""
    names: set[str] = set()
    for parameter in command.params:
        if isinstance(parameter, click.Option):
            names.update(parameter.opts)
            names.update(parameter.secondary_opts)
    return names


def test_public_cli_uses_the_full_thermoelasticity_name() -> None:
    """The pre-release CLI uses the same scientific noun as the API namespace."""
    assert "thermoelasticity" in main.commands
    assert "thermoelastic" not in main.commands


def test_all_public_cli_options_have_help_text() -> None:
    """Generated CLI reference never contains undocumented options."""
    missing: list[str] = []
    for path, command in _leaf_commands(main):
        for parameter in command.params:
            if isinstance(parameter, click.Option) and not parameter.help:
                names = parameter.opts or [parameter.name]
                missing.append(f"{' '.join(path)}: {names[0]}")
    assert not missing, missing


def test_removed_cli_options_are_absent_from_the_public_command_tree() -> None:
    """The pre-release migration leaves no historical option aliases behind."""
    for path, command in _leaf_commands(main):
        overlap = _REMOVED_OPTIONS & _option_names(command)
        assert overlap == set(), f"{' '.join(path)} still exposes {sorted(overlap)}"


def test_significant_run_commands_share_the_output_contract() -> None:
    """Scientific run commands expose one memorable reporting vocabulary."""
    runner = CliRunner()
    for path in (
        ("ha", "run"),
        ("qha", "run"),
        ("elasticity", "run"),
        ("seismic", "run"),
        ("eos", "run"),
        ("thermoelasticity", "run"),
    ):
        result = runner.invoke(main, [*path, "--help"])
        assert result.exit_code == 0, result.output
        assert f"{OUTPUT_GROUP}:" in result.output
        assert "-o, --output" in result.output
        assert "-r, --report" in result.output
        assert "-v, --verbosity" in result.output
        assert "-q, --quiet" in result.output


def test_qha_help_orders_science_before_domain_numerics_and_output() -> None:
    """The most option-rich workflow keeps a stable semantic help order."""
    result = CliRunner().invoke(main, ["qha", "run", "--help"])

    assert result.exit_code == 0
    positions = [
        result.output.index("Scientific model:"),
        result.output.index("Calculation domain:"),
        result.output.index("Numerical method:"),
        result.output.index("Validation and diagnostics:"),
        result.output.index("Units:"),
        result.output.index(f"{OUTPUT_GROUP}:"),
    ]
    assert positions == sorted(positions)

    assert "-D, --energy-degree" in result.output
    assert "-F, --frequency-degree" in result.output


def test_report_contract_uses_log_suffix_and_three_levels(tmp_path: Path) -> None:
    """Automatic reports are deterministic and verbosity is presentation-only."""
    source = tmp_path / "calculation.yaml"
    explicit = tmp_path / "custom-report.txt"

    assert default_report_path(source, None) == tmp_path / "calculation.log"
    assert default_report_path(source, explicit) == explicit
    assert [item.value for item in ReportVerbosity] == [
        "standard",
        "extended",
        "debug",
    ]


def test_ha_and_qha_reference_sets_are_scientifically_complete() -> None:
    """HA and QHA reports expose the approved canonical scientific sources."""
    assert MODULE_CITATION_KEYS["ha"] == (
        "quantas_2022",
        "mcquarrie_simon_1997",
    )
    assert MODULE_CITATION_KEYS["qha"] == (
        "quantas_2022",
        "qha_ulian_valdre_2018",
        "anderson_1995",
        "anderson_masuda_isaak_1995",
        "erba_2014",
        "erba_shahrokhi_moradian_dovesi_2015",
    )
    assert METHOD_CITATION_KEYS["harmonic_statistical_thermodynamics"] == (
        "mcquarrie_simon_1997",
    )
    assert METHOD_CITATION_KEYS["quasi_harmonic_approximation"] == (
        "anderson_1995",
        "anderson_masuda_isaak_1995",
        "erba_2014",
        "erba_shahrokhi_moradian_dovesi_2015",
    )


def test_thermoelasticity_reference_set_includes_general_qha_framework() -> None:
    """Thermoelastic reports cite the general QHA and finite-strain framework."""
    assert "stixrude_lithgow_bertelloni_2005" in MODULE_CITATION_KEYS[
        "thermoelasticity"
    ]
    assert "destefanis_ravoux_cossard_erba_2019" in MODULE_CITATION_KEYS[
        "thermoelasticity"
    ]
    assert METHOD_CITATION_KEYS["cold_finite_strain"] == (
        "stixrude_lithgow_bertelloni_2005",
    )
    assert METHOD_CITATION_KEYS["quasi_static_thermoelasticity"] == (
        "destefanis_ravoux_cossard_erba_2019",
        "stixrude_lithgow_bertelloni_2005",
    )
    citation = get_citation("stixrude_lithgow_bertelloni_2005")
    assert citation.year == 2005
    assert citation.volume == "162"
    assert citation.pages == "610-632"
    assert citation.doi == "10.1111/j.1365-246X.2005.02642.x"


def test_citation_sets_resolve_only_canonical_registry_keys() -> None:
    """Module and method references never embed independent bibliography text."""
    assert all(key == citation.key for key, citation in CITATIONS.items())
    for citation_set in (
        *MODULE_CITATION_KEYS.values(),
        *METHOD_CITATION_KEYS.values(),
    ):
        assert len(citation_set) == len(set(citation_set))
        for key in citation_set:
            assert get_citation(key) is CITATIONS[key]


def test_docs_branding_assets_are_packaged() -> None:
    """The Sphinx logo, banner, and favicon PNG files are packaged."""
    manifest = Path("MANIFEST.in").read_text(encoding="utf-8")
    documentation_rule = next(
        line
        for line in manifest.splitlines()
        if line.startswith("recursive-include docs/source ")
    )
    assert "*.png" in documentation_rule.split()


def test_odrpack_is_a_runtime_dependency() -> None:
    """ODR is part of the supported scientific installation, not a dev extra."""
    text = Path("pyproject.toml").read_text(encoding="utf-8")
    dependencies = text.split("dependencies = [", maxsplit=1)[1].split("]", maxsplit=1)[
        0
    ]
    assert '"odrpack>=0.6,<0.7"' in dependencies


def test_significant_run_writes_automatic_log_when_terminal_is_quiet(
    tmp_path: Path,
) -> None:
    """A scientific run always persists its report beside the primary input."""
    source = (
        Path(__file__).resolve().parents[2]
        / "tests"
        / "modules"
        / "ha"
        / "data"
        / "mgo_b3lyp_qha.yaml"
    )
    input_file = tmp_path / "mgo.yaml"
    shutil.copy2(source, input_file)
    output_file = tmp_path / "mgo_ha.hdf5"

    result = CliRunner().invoke(
        main,
        [
            "ha",
            "run",
            str(input_file),
            "--no-progress",
            "--quiet",
            "-o",
            str(output_file),
        ],
    )

    assert result.exit_code == 0, result.output
    assert result.output == ""
    assert output_file.is_file()
    report = input_file.with_suffix(".log")
    assert report.is_file()
    text = report.read_text(encoding="utf-8")
    assert "Selected options" in text
    assert "Thermodynamic properties" in text


def test_static_plot_commands_share_named_figure_presets() -> None:
    """Every Matplotlib frontend exposes the same high-level style vocabulary."""
    runner = CliRunner()
    commands = (
        ("elasticity", "run"),
        ("elasticity", "plot"),
        ("seismic", "plot"),
        ("ha", "run"),
        ("ha", "plot"),
        ("qha", "plot"),
        ("eos", "plot"),
        ("thermoelasticity", "run"),
        ("thermoelasticity", "plot", "fit"),
        ("thermoelasticity", "plot", "pt"),
        ("thermoelasticity", "plot", "domain"),
        ("thermoelasticity", "plot", "profile"),
        ("thermoelasticity", "plot", "compare"),
    )
    for command in commands:
        result = runner.invoke(main, [*command, "--help"])
        assert result.exit_code == 0, f"{' '.join(command)}\n{result.output}"
        option = "--plot-preset" if command[-1] == "run" else "--preset"
        assert option in result.output
        for preset in ("screen", "publication", "monochrome"):
            assert preset in result.output


def test_elasticity_seismic_references_are_foundational() -> None:
    """Elasticity and SEISMIC reports expose their foundational references."""
    assert MODULE_CITATION_KEYS["elasticity"] == (
        "quantas_2022",
        "elasticity_ulian_moro_valdre_2018",
        "nye_1985",
        "hill_1952",
        "elate_gaillac_pullumbi_coudert_2016",
        "mouhat_coudert_2014",
    )
    assert MODULE_CITATION_KEYS["seismic"] == (
        "quantas_2022",
        "seismic_ulian_valdre_2024",
        "jaeken_cottenier_2016",
        "nye_1985",
    )
    assert METHOD_CITATION_KEYS["voigt_reuss_hill"] == ("hill_1952",)
    assert METHOD_CITATION_KEYS["christoffel_acoustics"] == (
        "jaeken_cottenier_2016",
        "seismic_ulian_valdre_2024",
    )


def test_eos_method_references_follow_the_canonical_review() -> None:
    """EOS method sets cite Angel and the relevant thermodynamic framework."""
    assert METHOD_CITATION_KEYS["isothermal_equations_of_state"] == (
        "eosfit7_angel_gonzalez_platas_alvaro_2014",
        "anderson_1995",
    )
    assert METHOD_CITATION_KEYS["thermal_expansion_equations"] == (
        "eosfit7_angel_gonzalez_platas_alvaro_2014",
    )
    assert METHOD_CITATION_KEYS["pvt_equations_of_state"] == (
        "eosfit7_angel_gonzalez_platas_alvaro_2014",
        "anderson_1995",
    )
    assert METHOD_CITATION_KEYS["mie_gruneisen_debye"] == (
        "anderson_1995",
        "stixrude_lithgow_bertelloni_2005",
    )


def test_public_cli_commands_have_long_form_reference_help() -> None:
    """Every public command explains scope and argument semantics."""
    documented: list[str] = []

    def visit(command: click.Command, path: tuple[str, ...]) -> None:
        help_text = (command.help or "").strip()
        if len(help_text) < 100:
            documented.append(f"{' '.join(path)}: {len(help_text)} characters")
        if isinstance(command, click.Group):
            for name, child in command.commands.items():
                visit(child, (*path, name))

    visit(main, ("quantas",))
    assert not documented, documented


def test_public_cli_options_have_explanatory_reference_help() -> None:
    """Generated option entries contain more than a terse label."""
    insufficient: list[str] = []
    for path, command in _leaf_commands(main):
        for parameter in command.params:
            if not isinstance(parameter, click.Option):
                continue
            help_text = (parameter.help or "").strip()
            if len(help_text) < 80:
                names = parameter.opts or [parameter.name]
                insufficient.append(
                    f"{' '.join(path)}: {names[0]} ({len(help_text)} characters)"
                )
    assert not insufficient, insufficient


def test_long_form_help_is_available_from_direct_group_imports() -> None:
    """sphinx-click imports module groups directly and sees enriched text."""
    from quantas.cli.elasticity import elasticity
    from quantas.cli.eos import eos
    from quantas.cli.ha import ha
    from quantas.cli.qha import qha
    from quantas.cli.seismic import seismic
    from quantas.cli.thermoelastic import thermoelasticity

    groups = (elasticity, eos, ha, qha, seismic, thermoelasticity)
    for group in groups:
        assert len((group.help or "").strip()) >= 100
        for _, command in _leaf_commands(group):
            assert len((command.help or "").strip()) >= 100
            for parameter in command.params:
                if isinstance(parameter, click.Option):
                    assert len((parameter.help or "").strip()) >= 80
