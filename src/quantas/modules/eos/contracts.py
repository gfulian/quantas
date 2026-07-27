# -*- coding: utf-8 -*-

"""Stable public contracts and capability declarations for Quantas EOS.

The EOS workflow differs from single-shot Quantas modules because one native
archive may contain several datasets, immutable fit attempts, and accepted
results.  This module therefore defines an EOS-specific application contract
instead of forcing the workflow into :class:`quantas.models.ModuleContract`.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from enum import Enum
from typing import Any

from .hdf5 import (
    EOS_ARCHIVE_SCHEMA_VERSION,
    EOS_SUPPORTED_ARCHIVE_SCHEMA_VERSIONS,
)
from .models import EOSFitDomain




class EOSCapabilityStatus(str, Enum):
    """Availability level of one EOS scientific domain.

    Attributes
    ----------
    PUBLIC
        Supported by the stable EOS fitting and post-fit workflows.
    CORE_ONLY
        Available as a numerical core used by other Quantas modules, but not
        exposed as an EOS fitting workflow.
    PLANNED
        Reserved for a future public workflow.
    """

    PUBLIC = "public"
    CORE_ONLY = "core_only"
    PLANNED = "planned"


@dataclass(frozen=True, slots=True)
class EOSDomainCapability:
    """Describe the public status of one EOS scientific domain.

    Parameters
    ----------
    domain : EOSFitDomain
        Stable domain identifier.
    status : EOSCapabilityStatus
        Current support level.
    fitting : bool
        Whether :class:`EOSFitter` accepts the domain.
    calculator : bool
        Whether post-fit property evaluation is public.
    diagnostics : bool
        Whether post-fit scientific diagnostics are public.
    plotting : bool
        Whether neutral plot specifications are public.
    note : str
        Short explanation intended for API discovery and documentation.
    """

    domain: EOSFitDomain
    status: EOSCapabilityStatus
    fitting: bool
    calculator: bool
    diagnostics: bool
    plotting: bool
    note: str


EOS_DOMAIN_CAPABILITIES: tuple[EOSDomainCapability, ...] = (
    EOSDomainCapability(
        EOSFitDomain.PRESSURE_VOLUME,
        EOSCapabilityStatus.PUBLIC,
        True,
        True,
        True,
        True,
        "Isothermal volumetric and axial pressure EOS fitting.",
    ),
    EOSDomainCapability(
        EOSFitDomain.VOLUME_TEMPERATURE,
        EOSCapabilityStatus.PUBLIC,
        True,
        True,
        True,
        True,
        "Volumetric and axial thermal-expansion fitting.",
    ),
    EOSDomainCapability(
        EOSFitDomain.PRESSURE_VOLUME_TEMPERATURE,
        EOSCapabilityStatus.PUBLIC,
        True,
        True,
        True,
        True,
        "Coupled P-V-T fitting, including thermal-pressure models.",
    ),
    EOSDomainCapability(
        EOSFitDomain.ENERGY_VOLUME,
        EOSCapabilityStatus.CORE_ONLY,
        False,
        False,
        False,
        False,
        "Integrated energy EOS remain stable numerical core functionality used "
        "by QHA; no public EOS E-V fitting workflow is currently provided.",
    ),
)


def eos_domain_capability(
    domain: EOSFitDomain | str,
) -> EOSDomainCapability:
    """Return the declared capability for one EOS domain.

    Parameters
    ----------
    domain : EOSFitDomain or str
        Stable domain value such as ``"pv"`` or ``"ev"``.

    Returns
    -------
    EOSDomainCapability
        Immutable capability declaration.

    Raises
    ------
    ValueError
        If the domain is unknown.
    """
    resolved = EOSFitDomain(domain)
    for capability in EOS_DOMAIN_CAPABILITIES:
        if capability.domain is resolved:
            return capability
    raise ValueError(f"No EOS capability is declared for {resolved.value!r}")


@dataclass(frozen=True, slots=True)
class EOSModuleContract:
    """Stable application-facing entry points for the EOS subsystem.

    Parameters
    ----------
    name : str
        Stable module identifier.
    archive_schema_version : str
        Schema written by the current native EOS archive implementation.
    supported_archive_schema_versions : tuple of str
        Archive schemas accepted by current readers.
    capabilities : tuple of EOSDomainCapability
        Explicit scientific-domain support matrix.
    read_input, fit, run_batch, open_archive, calculate, diagnose,
    describe_plots, build_plots : callable
        Frontend-neutral public operations.  Their concrete technical
        signatures are documented by the referenced functions.
    """

    name: str
    archive_schema_version: str
    supported_archive_schema_versions: tuple[str, ...]
    capabilities: tuple[EOSDomainCapability, ...]
    read_input: Callable[..., Any]
    fit: Callable[..., Any]
    run_batch: Callable[..., Any]
    open_archive: Callable[..., Any]
    calculate: Callable[..., Any]
    diagnose: Callable[..., Any]
    describe_plots: Callable[..., Any]
    build_plots: Callable[..., Any]

    def capability(self, domain: EOSFitDomain | str) -> EOSDomainCapability:
        """Return the capability declaration for ``domain``."""
        return eos_domain_capability(domain)


__all__ = [
    "EOSCapabilityStatus",
    "EOSDomainCapability",
    "EOSModuleContract",
    "EOS_DOMAIN_CAPABILITIES",
    "eos_domain_capability",
    "EOS_ARCHIVE_SCHEMA_VERSION",
    "EOS_SUPPORTED_ARCHIVE_SCHEMA_VERSIONS",
]
