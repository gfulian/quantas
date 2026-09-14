# -*- coding: utf-8 -*-

"""Discover models and capabilities available to Quantas EOS workflows."""

from __future__ import annotations

import click

from quantas.cli.output import CLIOutput
from quantas.core.physics.eos import (
    EOSModel,
    PVTCouplingFamily,
    TemperatureEOSFamily,
    TemperatureEOSModel,
    available_eos_models,
    available_pvt_couplings,
    available_temperature_eos_models,
    thermal_pressure_model_contracts,
)
from quantas.models import ReportTable
from quantas.modules.eos.contracts import EOS_DOMAIN_CAPABILITIES

_DOMAIN_ORDER = ("pv", "ev", "vt", "pvt")
_DOMAIN_NAMES = {
    "pv": "P-V",
    "ev": "E-V",
    "vt": "V-T",
    "pvt": "P-V-T",
}


@click.command(name="show-models")
@click.option(
    "--domain",
    "domains",
    multiple=True,
    type=click.Choice(list(_DOMAIN_ORDER), case_sensitive=False),
    help=(
        "Show only the requested EOS domain. Repeat to show several domains; "
        "when omitted, all EOS domains are shown."
    ),
)
def show_models(domains: tuple[str, ...]) -> None:
    """Show EOS scientific domains and their available model catalogues.

    With no filter, the command reports all EOS domains: pressure-volume,
    energy-volume, volume-temperature, and coupled pressure-volume-temperature.
    Repeated ``--domain`` options reduce the output to the requested sections;
    they do not require one model to span several unrelated domains.
    """
    selected = _selected_domains(domains)
    output = CLIOutput()
    output.table(_domain_table(selected), persist=False)

    if "pv" in selected:
        output.table(_isothermal_table("pv"), persist=False)
    if "ev" in selected:
        output.table(_isothermal_table("ev"), persist=False)
    if "vt" in selected:
        output.table(_temperature_table(), persist=False)
    if "pvt" in selected:
        output.table(_pvt_coupling_table(), persist=False)
        output.table(_thermal_pressure_table(), persist=False)

    output.text_block(
        "Canonical model tags are stored in normalized workflow metadata. "
        "Common isothermal aliases include BM -> BM3, PT/NS -> PT3, "
        "V -> V3, T -> T3, and SJEOS -> SJ; long family names are also "
        "accepted. V-T and P-V-T models use their documented family/variant "
        "or coupling tags.\n"
        "EOS-valued CLI options provide model-aware completion after shell "
        "completion has been registered. On PowerShell run "
        "'quantas completion powershell | Out-String | Invoke-Expression' "
        "for the current session.",
        persist=False,
    )
    output.close()


def _selected_domains(domains: tuple[str, ...]) -> tuple[str, ...]:
    """Return requested domains in stable scientific order."""
    if not domains:
        return _DOMAIN_ORDER
    selected = {value.lower() for value in domains}
    return tuple(value for value in _DOMAIN_ORDER if value in selected)


def _domain_table(domains: tuple[str, ...]) -> ReportTable:
    """Return the public EOS domain-capability summary."""
    capabilities = {
        capability.domain.value: capability for capability in EOS_DOMAIN_CAPABILITIES
    }
    rows = []
    for domain in domains:
        capability = capabilities[domain]
        rows.append(
            [
                domain,
                _DOMAIN_NAMES[domain],
                capability.status.value.replace("_", "-"),
                _yes_no(capability.fitting),
                _yes_no(capability.calculator),
                _yes_no(capability.diagnostics),
                _yes_no(capability.plotting),
            ]
        )
    return ReportTable(
        title="EOS scientific domains",
        columns=["Domain", "Relationship", "Status", "Fit", "Calculate", "Diagnose", "Plot"],
        rows=rows,
    )


def _isothermal_table(domain: str) -> ReportTable:
    """Return one P-V or E-V isothermal-model table."""
    if domain == "pv":
        models = tuple(model for model in _catalog_models() if model.supports_pressure_fit)
        title = "P-V isothermal models"
    else:
        models = tuple(model for model in _catalog_models() if model.supports_energy)
        title = "E-V integrated energy models"
    rows = [
        [
            model.tag,
            model.family_name,
            "-" if model.order is None else model.order,
        ]
        for model in models
    ]
    return ReportTable(
        title=title,
        columns=["Tag", "Formulation", "Order"],
        rows=rows,
    )


def _temperature_table() -> ReportTable:
    """Return all V-T model family/variant combinations."""
    rows = [
        [model.tag, _temperature_family_name(model), model.variant.value]
        for model in available_temperature_eos_models()
        if model.variant is not None
    ]
    return ReportTable(
        title="V-T thermal-expansion models",
        columns=["Tag", "Formulation", "Variant"],
        rows=rows,
    )


def _pvt_coupling_table() -> ReportTable:
    """Return all P-V-T coupling prescriptions."""
    descriptions = {
        PVTCouplingFamily.LINEAR_BULK_MODULUS: (
            "Linear bulk modulus",
            "reference P-V EOS + V-T model",
        ),
        PVTCouplingFamily.ANDERSON_GRUNEISEN: (
            "Anderson-Gruneisen",
            "reference P-V EOS + V-T model",
        ),
        PVTCouplingFamily.THERMAL_PRESSURE: (
            "Thermal pressure",
            "reference P-V EOS + thermal-pressure model",
        ),
    }
    rows = [
        [coupling.value, *descriptions[coupling]]
        for coupling in available_pvt_couplings()
    ]
    return ReportTable(
        title="P-V-T coupling models",
        columns=["Tag", "Formulation", "Components"],
        rows=rows,
    )


def _thermal_pressure_table() -> ReportTable:
    """Return thermal-pressure components available to P-V-T coupling."""
    rows = []
    for model in thermal_pressure_model_contracts():
        family = model.family_name.value.replace("-", " ").title()
        variant = "-" if model.mgd_variant is None else model.mgd_variant.value
        rows.append([model.tag, family, variant])
    return ReportTable(
        title="P-V-T thermal-pressure components",
        columns=["Tag", "Formulation", "Variant"],
        rows=rows,
    )


def _temperature_family_name(model: TemperatureEOSModel) -> str:
    """Return a compact human-readable V-T family name."""
    names = {
        TemperatureEOSFamily.BERMAN: "Berman",
        TemperatureEOSFamily.FEI: "Fei",
        TemperatureEOSFamily.MODIFIED_HOLLAND_POWELL: "Modified Holland-Powell",
        TemperatureEOSFamily.SALJE: "Salje",
        TemperatureEOSFamily.KROLL_HOLLAND_POWELL: "Kroll-Holland-Powell",
    }
    return names[model.family]


def _catalog_models() -> tuple[EOSModel, ...]:
    """Return the stable union of direct P-V and integrated E-V models."""
    models: dict[str, EOSModel] = {}
    for model in available_eos_models():
        models.setdefault(model.tag, model)
    for model in available_eos_models(require_energy=True):
        models.setdefault(model.tag, model)
    return tuple(models.values())


def _yes_no(value: bool) -> str:
    """Return a compact yes/no capability label."""
    return "yes" if value else "no"


__all__ = ["show_models"]
