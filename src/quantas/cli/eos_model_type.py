# -*- coding: utf-8 -*-

"""Click validation and completion for isothermal EOS model tags.

The command-line interface accepts the compact historical Quantas tags such as
``BM3``, ``NS4``, ``V3``, ``T3``, and ``SJ`` together with the long aliases
understood by the scientific EOS resolver.  Values are normalized to canonical
:class:`~quantas.core.physics.eos.EOSModel` tags before they reach workflow
adapters, while shell completion remains concise and capability-aware.
"""

from __future__ import annotations

import re

import click
from click.shell_completion import CompletionItem

from quantas.core.physics.eos import (
    EOSFamily,
    EOSModel,
    available_eos_models,
    available_eos_tags,
    parse_eos_model,
)


class EOSModelParamType(click.ParamType):
    """Validate and normalize one isothermal EOS model specification.

    Parameters
    ----------
    require_energy : bool, optional
        Require an implemented volume-integrated ``E(V)`` form.
    require_pressure_fit : bool, optional
        Require exposure for direct standalone pressure-volume fitting.

    Notes
    -----
    The converted value is the canonical compact tag, not an :class:`EOSModel`
    instance.  This preserves the existing CLI-to-API contracts and persisted
    representations while centralizing validation in the shared scientific
    resolver.
    """

    name = "model"

    def __init__(
        self,
        *,
        require_energy: bool = False,
        require_pressure_fit: bool = False,
    ) -> None:
        self.require_energy = bool(require_energy)
        self.require_pressure_fit = bool(require_pressure_fit)

    def convert(
        self,
        value: object,
        param: click.Parameter | None,
        ctx: click.Context | None,
    ) -> str | None:
        """Return the canonical EOS tag for one command-line value.

        Parameters
        ----------
        value : object
            Raw Click value.
        param : click.Parameter or None
            Click parameter being converted.
        ctx : click.Context or None
            Active Click context.

        Returns
        -------
        str or None
            Canonical EOS tag such as ``BM3`` or ``SJ``.

        Raises
        ------
        click.BadParameter
            If the model is unknown or lacks a capability required by this
            option.
        """
        if value is None:
            return None
        if isinstance(value, EOSModel):
            model = value
            raw = model.tag
        else:
            raw = str(value).strip()
            try:
                model = parse_eos_model(raw)
            except ValueError as exc:
                self.fail(_model_error(raw), param, ctx)
                raise AssertionError("unreachable") from exc

        if self.require_energy and not model.supports_energy:
            self.fail(
                f"{model.tag} has no implemented E(V) form. "
                "Run 'quantas eos show-models --domain ev' to list compatible models.",
                param,
                ctx,
            )
        if self.require_pressure_fit and not model.supports_pressure_fit:
            self.fail(
                f"{model.name} is not available for direct P-V fitting. "
                "Its analytical P(V) derivative may still be available from an "
                "energy fit. Run 'quantas eos show-models --domain pv' to list "
                "direct P-V models.",
                param,
                ctx,
            )
        return model.tag

    def shell_complete(
        self,
        ctx: click.Context,
        param: click.Parameter,
        incomplete: str,
    ) -> list[CompletionItem]:
        """Return capability-aware shell-completion candidates.

        Parameters
        ----------
        ctx : click.Context
            Active Click context.
        param : click.Parameter
            Parameter requesting completion.
        incomplete : str
            Partially typed value.

        Returns
        -------
        list of click.shell_completion.CompletionItem
            Canonical compact tags and selected historical aliases with model
            descriptions suitable for interactive shell completion.
        """
        del ctx, param
        prefix = incomplete.strip().lower()
        items: list[CompletionItem] = []
        seen: set[str] = set()
        for value, model, alias in _completion_candidates(
            require_energy=self.require_energy,
            require_pressure_fit=self.require_pressure_fit,
        ):
            if not value.lower().startswith(prefix) or value.lower() in seen:
                continue
            help_text = model.name if not alias else f"alias for {model.tag}: {model.name}"
            items.append(CompletionItem(value, help=help_text))
            seen.add(value.lower())
        return items


def _catalog_models() -> tuple[EOSModel, ...]:
    """Return the stable union of direct P-V and integrated E-V models."""
    models: dict[str, EOSModel] = {}
    for model in available_eos_models():
        models.setdefault(model.tag, model)
    for model in available_eos_models(require_energy=True):
        models.setdefault(model.tag, model)
    return tuple(models.values())


def _models_for_capabilities(
    *,
    require_energy: bool,
    require_pressure_fit: bool,
) -> tuple[EOSModel, ...]:
    """Return catalogue models satisfying the requested capabilities."""
    return tuple(
        model
        for model in _catalog_models()
        if (not require_energy or model.supports_energy)
        and (not require_pressure_fit or model.supports_pressure_fit)
    )


def _completion_candidates(
    *,
    require_energy: bool,
    require_pressure_fit: bool,
) -> tuple[tuple[str, EOSModel, bool], ...]:
    """Return canonical and concise historical completion aliases."""
    models = _models_for_capabilities(
        require_energy=require_energy,
        require_pressure_fit=require_pressure_fit,
    )
    by_tag = {model.tag: model for model in models}
    candidates: list[tuple[str, EOSModel, bool]] = []

    for tag in available_eos_tags(
        require_energy=require_energy,
        include_default_aliases=True,
    ):
        try:
            model = parse_eos_model(tag)
        except ValueError:
            continue
        if model.tag in by_tag:
            candidates.append((tag, by_tag[model.tag], tag != model.tag))

    # Historical Natural-Strain and descriptive SJEOS aliases are especially
    # useful at the prompt, but need not inflate normal ``--help`` output.
    for value in ("NS", "NS2", "NS3", "NS4", "SJEOS"):
        try:
            model = parse_eos_model(value)
        except ValueError:
            continue
        if model.tag in by_tag:
            candidates.append((value, by_tag[model.tag], value != model.tag))

    return tuple(candidates)


def _model_error(raw: str) -> str:
    """Return a concise model-resolution error with recovery guidance."""
    text = raw.strip()
    match = re.fullmatch(r"(.+?)(\d+)", text)
    if match is not None:
        base, order_text = match.groups()
        try:
            base_model = parse_eos_model(base)
        except ValueError:
            base_model = None
        if base_model is not None:
            family_models = tuple(
                model for model in _catalog_models() if model.family is base_model.family
            )
            ordered = tuple(model for model in family_models if model.order is not None)
            if not ordered:
                preferred = _orderless_usage(base_model.family)
                return (
                    f"{base_model.family_name} does not define selectable EOS orders. "
                    f"Use {preferred}."
                )
            choices = ", ".join(model.tag for model in ordered)
            return (
                f"Unsupported {base_model.family_name} order {order_text}. "
                f"Available models: {choices}."
            )
    return (
        f"Unknown EOS model {text!r}. "
        "Run 'quantas eos show-models' to list the available models."
    )


def _orderless_usage(family: EOSFamily) -> str:
    """Return concise accepted tags for an orderless family."""
    if family is EOSFamily.STABILIZED_JELLIUM:
        return "SJ or SJEOS"
    if family is EOSFamily.MURNAGHAN:
        return "M or Murnaghan"
    return family.value


ENERGY_EOS_MODEL = EOSModelParamType(require_energy=True)
"""Reusable Click type for options requiring a volume-integrated EOS."""


__all__ = ["ENERGY_EOS_MODEL", "EOSModelParamType"]
