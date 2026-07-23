# -*- coding: utf-8 -*-

"""Declarative batch execution for EOS fitting workflows.

The classes in this module coordinate existing EOS readers, fitters, archives,
and event observers.  They contain no Click, Rich, terminal-stream, or plotting
logic and are therefore reusable by the Python API, command line, and a future
GUI.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any

from quantas.core.events import Event, EventLevel, EventRecord, NullObserver
from quantas.core.math.fitting import FitResult, FitStatus

from .api import EOSFitter
from .archive import EOSArchive
from .history import EOSResultSlot
from .io import read_eos_input
from .models import EOSDataset, EOSFitRequest, EOSFitResult


class EOSBatchFailurePolicy(str, Enum):
    """Action taken after a failed batch job."""

    STOP = "stop"
    CONTINUE = "continue"


@dataclass(frozen=True, slots=True)
class EOSBatchJob:
    """One ordered EOS fit requested by a batch plan.

    Parameters
    ----------
    request : EOSFitRequest
        Complete scientific and numerical fit request.
    job_id : str or None, optional
        Stable human-readable identifier.  The request identifier is used when
        omitted, followed by a generated positional identifier.
    accept : bool, optional
        Accept a successful record as the current result for its domain/target
        slot.
    replace_accepted : bool, optional
        Allow this job to replace an earlier accepted result in the same slot.
        Replacement is never implicit.
    note : str or None, optional
        Optional scientific note stored with the fit record.
    metadata : dict, optional
        Passive provenance attached to the fit record.
    """

    request: EOSFitRequest
    job_id: str | None = None
    accept: bool = True
    replace_accepted: bool = False
    note: str | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not isinstance(self.request, EOSFitRequest):
            raise TypeError("EOSBatchJob request must be an EOSFitRequest")
        object.__setattr__(
            self, "job_id", None if self.job_id is None else str(self.job_id)
        )
        object.__setattr__(self, "accept", bool(self.accept))
        object.__setattr__(self, "replace_accepted", bool(self.replace_accepted))
        object.__setattr__(self, "note", None if self.note is None else str(self.note))
        object.__setattr__(self, "metadata", dict(self.metadata))

    @property
    def slot(self) -> EOSResultSlot:
        """Return the result slot addressed by the job."""
        return EOSResultSlot.from_request(self.request)

    def as_dict(self) -> dict[str, Any]:
        """Return a serialization-ready job description."""
        return {
            "job_id": self.job_id,
            "slot": self.slot.as_dict(),
            "accept": self.accept,
            "replace_accepted": self.replace_accepted,
            "note": self.note,
            "metadata": dict(self.metadata),
            "request": self.request.as_dict(),
        }


@dataclass(frozen=True, slots=True)
class EOSBatchPlan:
    """Ordered declarative collection of EOS batch jobs."""

    jobs: tuple[EOSBatchJob, ...] = ()
    failure_policy: EOSBatchFailurePolicy = EOSBatchFailurePolicy.STOP
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        jobs = tuple(self.jobs)
        if not all(isinstance(job, EOSBatchJob) for job in jobs):
            raise TypeError("EOSBatchPlan jobs must contain EOSBatchJob objects")
        object.__setattr__(self, "jobs", jobs)
        object.__setattr__(
            self, "failure_policy", EOSBatchFailurePolicy(self.failure_policy)
        )
        object.__setattr__(self, "metadata", dict(self.metadata))
        self._validate_acceptance_order()

    def _validate_acceptance_order(self) -> None:
        accepted: set[str] = set()
        for job in self.jobs:
            if not job.accept:
                continue
            key = job.slot.key
            if key in accepted and not job.replace_accepted:
                raise ValueError(
                    f"batch plan accepts more than one result for {key!r}; "
                    "set replace_accepted=True on the later job"
                )
            accepted.add(key)

    def as_dict(self) -> dict[str, Any]:
        """Return a serialization-ready batch manifest."""
        return {
            "failure_policy": self.failure_policy.value,
            "metadata": dict(self.metadata),
            "jobs": [job.as_dict() for job in self.jobs],
        }


@dataclass(frozen=True, slots=True)
class EOSBatchJobResult:
    """Outcome of one ordered batch job."""

    job_id: str
    request: EOSFitRequest
    result: EOSFitResult
    record_id: int
    accepted: bool

    @property
    def successful(self) -> bool:
        """Return whether the numerical fit succeeded."""
        return bool(self.result.fit.success)


@dataclass(slots=True)
class EOSBatchResult:
    """Frontend-neutral result of one complete batch execution."""

    dataset: EOSDataset
    plan: EOSBatchPlan
    archive_path: Path
    jobs: tuple[EOSBatchJobResult, ...]
    events: tuple[EventRecord, ...] = ()

    @property
    def successful(self) -> bool:
        """Return whether every executed job succeeded."""
        return all(job.successful for job in self.jobs)


class EOSBatchWorkflow:
    """Execute a declarative EOS batch plan and persist every attempt.

    Parameters
    ----------
    fitter : EOSFitter or None, optional
        Scientific fitting service.  A default instance is created when
        omitted.
    observer : callable or None, optional
        Event observer used by frontends.  Progress events are operational and
        are not retained in the returned event history.
    """

    def __init__(
        self,
        fitter: EOSFitter | None = None,
        observer: Callable[[Event], None] | None = None,
    ) -> None:
        self.fitter = EOSFitter() if fitter is None else fitter
        self.observer = NullObserver() if observer is None else observer
        self._events: list[EventRecord] = []

    def run(
        self,
        source: str | Path | EOSDataset,
        plan: EOSBatchPlan,
        archive_path: str | Path,
        *,
        overwrite: bool = False,
        creator: str = "quantas eos run",
        pressure_unit: str | None = None,
        length_unit: str | None = None,
        temperature_unit: str | None = None,
    ) -> EOSBatchResult:
        """Execute the plan and return structured results.

        Parameters
        ----------
        source : str, Path, or EOSDataset
            Keyword-directed data file or preconstructed normalized dataset.
        plan : EOSBatchPlan
            Ordered fit requests and explicit acceptance policy.
        archive_path : str or Path
            Native EOS HDF5 destination.
        overwrite : bool, optional
            Replace an existing destination.
        creator : str, optional
            Archive creator provenance.
        pressure_unit, length_unit, temperature_unit : str or None, optional
            Input-unit overrides applied only while reading a text source.

        Returns
        -------
        EOSBatchResult
            Dataset, records, archive path, and meaningful workflow events.
        """
        self._events = []
        dataset = self._read_source(
            source,
            pressure_unit=pressure_unit,
            length_unit=length_unit,
            temperature_unit=temperature_unit,
        )
        self._emit(
            f"Read EOS dataset '{dataset.jobname}' with {dataset.npoints} observations "
            f"({dataset.selected_npoints} selected by default)",
            data={"kind": "dataset", "dataset": dataset},
        )
        self._validate_plan(dataset, plan)
        destination = Path(archive_path)
        results: list[EOSBatchJobResult] = []
        with EOSArchive.create(
            destination,
            dataset=dataset,
            creator=creator,
            overwrite=overwrite,
        ) as archive:
            archive.write_batch_manifest(plan.as_dict())
            total = len(plan.jobs)
            for index, job in enumerate(plan.jobs, start=1):
                job_id = job.job_id or job.request.request_id or f"job-{index:06d}"
                self._emit(
                    f"Starting EOS batch job {index}/{total}: {job_id}",
                    data={
                        "kind": "job_start",
                        "job_id": job_id,
                        "request": job.request,
                    },
                )
                result = self._execute_job(dataset, job.request)
                record = archive.append_fit(
                    1,
                    result.request,
                    result,
                    note=job.note,
                    provenance={"workflow": "batch", "job_id": job_id, **job.metadata},
                )
                accepted = False
                if result.fit.success and job.accept:
                    current = archive.slot_state(record.slot)
                    if (
                        current.accepted_record_id is not None
                        and not job.replace_accepted
                    ):
                        raise ValueError(
                            f"job {job_id!r} would replace the accepted result for "
                            f"{record.slot.key!r} without explicit permission"
                        )
                    archive.accept(record.record_id)
                    accepted = True
                item = EOSBatchJobResult(
                    job_id=job_id,
                    request=result.request,
                    result=result,
                    record_id=record.record_id,
                    accepted=accepted,
                )
                results.append(item)
                level = EventLevel.RESULT if result.fit.success else EventLevel.WARNING
                self._emit(
                    f"Finished EOS batch job {job_id}: {result.fit.status.value}",
                    level=level,
                    data={"kind": "job_result", "job": item},
                )
                if total:
                    self._emit(
                        f"Completed {index} of {total} EOS batch jobs",
                        level=EventLevel.PROGRESS,
                        progress=index / total,
                        data={
                            "kind": "batch_progress",
                            "current": index,
                            "total": total,
                        },
                    )
                if (
                    not result.fit.success
                    and plan.failure_policy is EOSBatchFailurePolicy.STOP
                ):
                    break
        self._emit(
            f"EOS archive written to {destination}",
            data={"kind": "archive", "path": str(destination)},
        )
        return EOSBatchResult(
            dataset=dataset,
            plan=plan,
            archive_path=destination,
            jobs=tuple(results),
            events=tuple(self._events),
        )

    def _read_source(
        self,
        source: str | Path | EOSDataset,
        **unit_overrides: str | None,
    ) -> EOSDataset:
        if isinstance(source, EOSDataset):
            return source
        return read_eos_input(source, **unit_overrides)

    def _validate_plan(self, dataset: EOSDataset, plan: EOSBatchPlan) -> None:
        available = {state.key for state in _available_slots(dataset)}
        for index, job in enumerate(plan.jobs, start=1):
            if job.slot.key not in available:
                raise ValueError(
                    f"batch job {index} requests unavailable slot {job.slot.key!r}; "
                    f"available slots are {', '.join(sorted(available)) or 'none'}"
                )

    def _execute_job(self, dataset: EOSDataset, request: EOSFitRequest) -> EOSFitResult:
        try:
            return self.fitter.fit(dataset, request)
        except Exception as exc:
            fit = FitResult.failed(
                str(exc),
                status=FitStatus.INVALID_INPUT,
                method=(
                    None
                    if request.options.solver_options is None
                    else request.options.solver_options.method
                ),
            )
            return EOSFitResult(
                request=request,
                fit=fit,
                warnings=[str(exc)],
                metadata={"workflow_exception": type(exc).__name__},
            )

    def _emit(
        self,
        message: str,
        *,
        level: EventLevel = EventLevel.INFO,
        progress: float | None = None,
        data: dict[str, Any] | None = None,
    ) -> None:
        event = Event(message=message, level=level, progress=progress, data=data or {})
        if level is not EventLevel.PROGRESS:
            self._events.append(EventRecord.from_event(event))
        self.observer(event)


def _available_slots(dataset: EOSDataset) -> tuple[EOSResultSlot, ...]:
    from .archive import infer_result_slots

    return infer_result_slots(dataset)


def run_eos_batch(
    source: str | Path | EOSDataset,
    plan: EOSBatchPlan,
    archive_path: str | Path,
    **kwargs: Any,
) -> EOSBatchResult:
    """Execute one EOS batch using the default workflow service."""
    return EOSBatchWorkflow().run(source, plan, archive_path, **kwargs)


__all__ = [
    "EOSBatchFailurePolicy",
    "EOSBatchJob",
    "EOSBatchJobResult",
    "EOSBatchPlan",
    "EOSBatchResult",
    "EOSBatchWorkflow",
    "run_eos_batch",
]
