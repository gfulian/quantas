# -*- coding: utf-8 -*-

"""Dataset-dependent semantic resolution for EOS specification documents."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from quantas.core.physics.eos import (
    PVTCouplingFamily,
    PVTModel,
    ThermalPressureFamily,
    parse_eos_model,
    parse_pvt_coupling,
    parse_temperature_eos_model,
    parse_thermal_pressure_model,
)

from .batch import EOSBatchFailurePolicy, EOSBatchJob, EOSBatchPlan
from .models import (
    EOSDataset,
    EOSFitDomain,
    EOSFitOptions,
    EOSFitRequest,
)
from .report import EOSReportDetail, EOSReportOptions
from .spec import EOSResolvedSpec, EOSSpecDocument, EOSSpecError, _Entry, _Section
from ._spec_fit import (
    _build_constraints,
    _build_mgd_normalization,
    _build_solver_options,
    _reject_mgd_normalization_keys,
)
from ._spec_selection import _resolve_selection, _resolve_targets
from ._spec_values import (
    _entry_error,
    _enum_value,
    _merge_entries,
    _parse_bool_entry,
    _parse_domain,
    _positive_int,
    canonical_parameter_name,
)


def resolve_eos_spec_document(
    document: EOSSpecDocument,
    dataset: EOSDataset,
) -> EOSResolvedSpec:
    """Resolve defaults, targets, models, constraints, and solver contracts.

    Parameters
    ----------
    document : EOSSpecDocument
        Parsed EOS specification.
    dataset : EOSDataset
        Normalized dataset used to expand ``targets = all`` and validate target
        availability.

    Returns
    -------
    EOSResolvedSpec
        Typed batch plan and report options.

    Raises
    ------
    EOSSpecError
        If inherited settings are inconsistent or unavailable for ``dataset``.
    """
    ordinary, job_sections = _partition_sections(document)
    jobs = _resolve_jobs(document, dataset, ordinary, job_sections)
    failure_policy = _resolve_failure_policy(ordinary.get("batch"), document.source)
    try:
        plan = EOSBatchPlan(
            jobs=jobs,
            failure_policy=failure_policy,
            metadata={
                "source": "spec",
                "spec": document.as_dict(),
                **document.metadata,
            },
        )
    except ValueError as exc:
        raise EOSSpecError(str(exc), source=document.source) from exc
    report_options = _parse_report_options(
        ordinary.get("presentation"), document.source
    )
    return EOSResolvedSpec(document=document, plan=plan, report_options=report_options)


def _partition_sections(
    document: EOSSpecDocument,
) -> tuple[dict[str, _Section], tuple[_Section, ...]]:
    ordinary = {
        section.name: section
        for section in document._sections
        if not section.name.startswith("job ")
    }
    jobs = tuple(
        section for section in document._sections if section.name.startswith("job ")
    )
    if not jobs:
        raise EOSSpecError(
            "the specification must contain at least one [job NAME] section",
            source=document.source,
        )
    return ordinary, jobs


def _resolve_failure_policy(
    section: _Section | None,
    source: Path | None,
) -> EOSBatchFailurePolicy:
    if section is None:
        return EOSBatchFailurePolicy.STOP
    entry = section.mapping().get("failure_policy")
    if entry is None:
        return EOSBatchFailurePolicy.STOP
    return _enum_value(
        entry,
        {
            "stop": EOSBatchFailurePolicy.STOP,
            "continue": EOSBatchFailurePolicy.CONTINUE,
        },
        source,
        "failure_policy",
    )


def _resolve_jobs(
    document: EOSSpecDocument,
    dataset: EOSDataset,
    ordinary: dict[str, _Section],
    sections: tuple[_Section, ...],
) -> tuple[EOSBatchJob, ...]:
    jobs: list[EOSBatchJob] = []
    job_ids: set[str] = set()
    for section in sections:
        own = section.mapping()
        domain_entry = own.get("domain")
        if domain_entry is None:
            raise _entry_error(
                "missing required key 'domain'", section, document.source
            )
        domain = _parse_domain(domain_entry, document.source, section.display_name)
        target_entry = own.get("targets")
        if target_entry is None:
            raise _entry_error(
                "missing required key 'targets'", section, document.source
            )
        targets = _resolve_targets(
            target_entry, domain, dataset, document.source, section.display_name
        )
        inherited = _merge_entries(
            ordinary.get("defaults"),
            ordinary.get(f"defaults.{domain.value}"),
            section,
        )
        base_name = section.name[4:].strip()
        for target in targets:
            job_id = base_name if len(targets) == 1 else f"{base_name}-{target}"
            if job_id in job_ids:
                raise EOSSpecError(
                    f"duplicate resolved job identifier {job_id!r}",
                    source=document.source,
                    line=section.line,
                    section=section.display_name,
                )
            job_ids.add(job_id)
            request = _build_request(
                inherited,
                domain=domain,
                target=target,
                request_id=job_id,
                source=document.source,
                section=section.display_name,
                spec_hash=document.source_sha256,
                dataset=dataset,
            )
            jobs.append(
                EOSBatchJob(
                    request=request,
                    job_id=job_id,
                    accept=_parse_bool_entry(
                        inherited.get("accept"),
                        True,
                        document.source,
                        section.display_name,
                    ),
                    replace_accepted=_parse_bool_entry(
                        inherited.get("replace_accepted"),
                        False,
                        document.source,
                        section.display_name,
                    ),
                    note=(
                        None
                        if inherited.get("note") is None
                        else inherited["note"].value
                    ),
                    metadata={
                        "source": "spec",
                        "spec_section": section.display_name,
                        "spec_sha256": document.source_sha256,
                    },
                )
            )
    return tuple(jobs)


def _parse_report_options(
    section: _Section | None, source: Path | None
) -> EOSReportOptions:
    if section is None:
        return EOSReportOptions()
    values = section.mapping()
    detail = EOSReportDetail.SHORT
    if values.get("detail") is not None:
        entry = values["detail"]
        try:
            detail = EOSReportDetail(entry.value.strip().lower())
        except ValueError as exc:
            raise EOSSpecError(
                "detail must be 'short' or 'extended'",
                source=source,
                line=entry.line,
                section=section.display_name,
            ) from exc
    show = _parse_bool_entry(
        values.get("show_uncertainties"), False, source, section.display_name
    )
    maximum: int | None = None
    if values.get("max_data_rows") is not None:
        entry = values["max_data_rows"]
        if entry.value.strip().lower() not in {"all", "none"}:
            maximum = _positive_int(entry, source, section.display_name)
    return EOSReportOptions(
        detail=detail, show_uncertainties=show, max_data_rows=maximum
    )


def _build_request(
    entries: dict[str, _Entry],
    *,
    domain: EOSFitDomain,
    target: str,
    request_id: str,
    source: Path | None,
    section: str,
    spec_hash: str,
    dataset: EOSDataset,
) -> EOSFitRequest:
    model: Any
    try:
        if domain is EOSFitDomain.PRESSURE_VOLUME:
            forbidden = [
                name for name in ("pv_model", "vt_model", "coupling") if name in entries
            ]
            if forbidden:
                raise EOSSpecError(
                    f"key(s) {', '.join(forbidden)} are not valid for a pv job",
                    source=source,
                    line=entries[forbidden[0]].line,
                    section=section,
                )
            model_entry = entries.get("model")
            model = parse_eos_model("BM3" if model_entry is None else model_entry.value)
        elif domain is EOSFitDomain.VOLUME_TEMPERATURE:
            forbidden = [
                name for name in ("pv_model", "vt_model", "coupling") if name in entries
            ]
            if forbidden:
                raise EOSSpecError(
                    f"key(s) {', '.join(forbidden)} are not valid for a vt job",
                    source=source,
                    line=entries[forbidden[0]].line,
                    section=section,
                )
            model_entry = entries.get("model")
            model = parse_temperature_eos_model(
                "berman:quadratic" if model_entry is None else model_entry.value
            )
        else:
            if "model" in entries:
                raise EOSSpecError(
                    "P-V-T jobs use pv_model, vt_model, and coupling instead of model",
                    source=source,
                    line=entries["model"].line,
                    section=section,
                )
            pv_entry = entries.get("pv_model")
            vt_entry = entries.get("vt_model")
            coupling_entry = entries.get("coupling")
            coupling = "linear" if coupling_entry is None else coupling_entry.value
            coupling_family = parse_pvt_coupling(coupling)
            thermal = None
            if coupling_family is not PVTCouplingFamily.THERMAL_PRESSURE:
                thermal = parse_temperature_eos_model(
                    "berman:quadratic" if vt_entry is None else vt_entry.value
                )
            elif vt_entry is not None:
                raise EOSSpecError(
                    "thermal-pressure coupling defines its own thermal model; remove vt_model",
                    source=source,
                    line=vt_entry.line,
                    section=section,
                )
            thermal_pressure_model = None
            normalization = None
            if coupling_family is PVTCouplingFamily.THERMAL_PRESSURE:
                thermal_pressure_entry = entries.get("thermal_pressure_model")
                thermal_pressure_model = parse_thermal_pressure_model(
                    "holland-powell-einstein"
                    if thermal_pressure_entry is None
                    else thermal_pressure_entry.value
                )
                if (
                    thermal_pressure_model.family_name
                    is ThermalPressureFamily.MIE_GRUNEISEN_DEBYE
                ):
                    normalization = _build_mgd_normalization(
                        entries, source=source, section=section
                    )
                else:
                    _reject_mgd_normalization_keys(
                        entries, source=source, section=section
                    )
            else:
                if "thermal_pressure_model" in entries:
                    entry = entries["thermal_pressure_model"]
                    raise EOSSpecError(
                        "thermal_pressure_model is valid only for thermal-pressure coupling",
                        source=source,
                        line=entry.line,
                        section=section,
                    )
                _reject_mgd_normalization_keys(entries, source=source, section=section)
            model = PVTModel(
                pressure_model=parse_eos_model(
                    "BM3" if pv_entry is None else pv_entry.value
                ),
                temperature_model=thermal,
                coupling=coupling_family,
                thermal_pressure_model=thermal_pressure_model,
                mgd_normalization=normalization,
            )
    except EOSSpecError:
        raise
    except (TypeError, ValueError) as exc:
        candidates = [
            entries.get(name) for name in ("model", "pv_model", "vt_model", "coupling")
        ]
        location = next((entry for entry in candidates if entry is not None), None)
        raise EOSSpecError(
            str(exc),
            source=source,
            line=None if location is None else location.line,
            section=section,
        ) from exc

    selection_mask, selection_metadata = _resolve_selection(
        entries, dataset, source=source, section=section
    )
    solver_options = _build_solver_options(entries, source, section)
    constraints = _build_constraints(
        entries,
        source,
        section,
        domain=domain,
        target=target,
        model=model,
    )
    return EOSFitRequest(
        model=model,
        target=target,
        domain=domain,
        constraints=constraints,
        options=EOSFitOptions(solver_options=solver_options),
        mask=selection_mask,
        request_id=request_id,
        metadata={
            "source": "spec",
            "spec_section": section,
            "spec_sha256": spec_hash,
            "selection": selection_metadata,
        },
    )

__all__ = ["canonical_parameter_name", "resolve_eos_spec_document"]
