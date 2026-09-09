# -*- coding: utf-8 -*-

"""Discover equations of state available to Quantas workflows."""

from __future__ import annotations

import click

from quantas.cli.output import CLIOutput
from quantas.core.physics.eos import EOSModel, available_eos_models
from quantas.models import ReportTable


@click.command(name="show-models")
@click.option(
    "--domain",
    "domains",
    multiple=True,
    type=click.Choice(["pv", "ev"], case_sensitive=False),
    help=(
        "Restrict the catalogue to a scientific domain. Repeat to require "
        "support in every selected domain."
    ),
)
def show_models(domains: tuple[str, ...]) -> None:
    """Show the isothermal EOS catalogue and workflow capabilities.

    With no filter, the command lists the union of direct pressure-volume and
    integrated energy-volume models. Repeated ``--domain`` options select the
    intersection of the requested capabilities, which is useful when choosing
    one formulation for several interoperable workflows.
    """
    selected = {value.lower() for value in domains}
    models = tuple(
        model
        for model in _catalog_models()
        if ("pv" not in selected or model.supports_pressure_fit)
        and ("ev" not in selected or model.supports_energy)
    )
    title = _catalog_title(selected)
    rows = [
        [
            model.tag,
            model.family_name,
            "-" if model.order is None else model.order,
            "yes" if model.supports_pressure_fit else "no",
            "yes" if model.supports_energy else "no",
        ]
        for model in models
    ]
    output = CLIOutput()
    output.table(
        ReportTable(
            title=title,
            columns=["Tag", "Formulation", "Order", "P-V fit", "E(V)"],
            rows=rows,
        ),
        persist=False,
    )
    output.text_block(
        "Canonical tags are stored in normalized workflow metadata. "
        "Common aliases include BM -> BM3, PT/NS -> PT3, V -> V3, "
        "T -> T3, and SJEOS -> SJ; long family names such as "
        "birch-murnaghan and natural-strain are also accepted.\n"
        "P-V fit means direct standalone pressure-volume fitting. A model may "
        "still provide analytical P(V) as the derivative of an E(V) form.\n"
        "EOS-valued CLI options support shell completion when Click completion "
        "is enabled for the active shell.",
        persist=False,
    )
    output.close()


def _catalog_models() -> tuple[EOSModel, ...]:
    """Return the stable union of direct P-V and integrated E-V models."""
    models: dict[str, EOSModel] = {}
    for model in available_eos_models():
        models.setdefault(model.tag, model)
    for model in available_eos_models(require_energy=True):
        models.setdefault(model.tag, model)
    return tuple(models.values())


def _catalog_title(domains: set[str]) -> str:
    """Return a readable title for one domain-filtered catalogue."""
    if not domains:
        return "Available equations of state"
    names = {"pv": "P-V", "ev": "E-V"}
    selected = " and ".join(names[item] for item in ("pv", "ev") if item in domains)
    return f"EOS models available for {selected}"


__all__ = ["show_models"]
