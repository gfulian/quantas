# -*- coding: utf-8 -*-

"""Live Rich progress rendering for Quantas command-line workflows."""

from __future__ import annotations

import os
from typing import Any

from rich.console import Console, RenderableType
from rich.text import Text
from rich.progress import (
    BarColumn,
    Progress,
    SpinnerColumn,
    ProgressColumn,
    Task,
    TaskID,
    TaskProgressColumn,
    TextColumn,
)

from quantas.core.events import Event


class _PlainTimeRemainingColumn(ProgressColumn):
    """Render an unstyled compact remaining-time estimate."""

    def render(self, task: Task) -> RenderableType:
        """Return the task remaining time without additional color."""
        remaining = task.time_remaining
        if remaining is None:
            return Text("--:--")
        seconds = max(0, int(remaining))
        hours, remainder = divmod(seconds, 3600)
        minutes, seconds = divmod(remainder, 60)
        if hours:
            return Text(f"{hours:d}:{minutes:02d}:{seconds:02d}")
        return Text(f"{minutes:02d}:{seconds:02d}")


class RichProgressDisplay:
    """Render one live workflow progress task on a Rich console.

    The display is intentionally frontend-only. It consumes generic Quantas
    progress events and never writes to report files or scientific results.
    A new task is created automatically when a workflow resets its numerical
    counter, as happens between distinct directional grids.

    Parameters
    ----------
    console : rich.console.Console
        Console used by all other CLI rendering.
    enabled : bool, optional
        Whether live progress is requested. The display is still disabled for
        non-interactive or explicitly non-ANSI terminals.
    transient : bool, optional
        Whether completed bars are removed from the terminal.
    """

    def __init__(
        self,
        console: Console,
        *,
        enabled: bool = True,
        transient: bool = True,
    ) -> None:
        self.console = console
        self.enabled = bool(enabled) and _supports_live_progress(console)
        self.transient = bool(transient)
        self._progress: Progress | None = None
        self._task_id: TaskID | None = None
        self._last_current = 0.0
        self._last_total = 0.0

    @property
    def active(self) -> bool:
        """Return whether a live task is currently displayed."""
        return self._progress is not None and self._task_id is not None

    def update(self, event: Event) -> None:
        """Update the live display from one Quantas progress event.

        Parameters
        ----------
        event : Event
            Progress event containing a fractional progress value and,
            preferably, ``current`` and ``total`` payload fields.
        """
        if not self.enabled:
            return
        label, current, total = _progress_values(event)
        if total <= 0.0:
            total = 100.0
        current = min(max(current, 0.0), total)

        reset = (
            self._task_id is None
            or current < self._last_current
            or (not _same_total(total, self._last_total) and current <= 1.0)
        )
        if reset:
            self._replace_task(label, total, current)
        else:
            assert self._progress is not None
            assert self._task_id is not None
            self._progress.update(
                self._task_id,
                description=label,
                total=total,
                completed=int(round(current)),
                refresh=True,
            )

        self._last_current = current
        self._last_total = total
        if current >= total or (event.progress is not None and event.progress >= 1.0):
            self.finish()

    def finish(self) -> None:
        """Complete and stop the current live task, if any."""
        if self._progress is None:
            return
        if self._task_id is not None and self._last_total > 0.0:
            self._progress.update(
                self._task_id,
                completed=self._last_total,
                refresh=True,
            )
        self.stop()

    def stop(self) -> None:
        """Stop the current live task without claiming completion."""
        if self._progress is None:
            return
        self._progress.stop()
        self._progress = None
        self._task_id = None
        self._last_current = 0.0
        self._last_total = 0.0

    close = stop

    def _replace_task(self, label: str, total: float, current: float) -> None:
        """Replace the current task with a new progress sequence."""
        if self._progress is not None:
            self.stop()
        self._progress = _build_progress(
            self.console,
            transient=self.transient,
        )
        self._progress.start()
        self._task_id = self._progress.add_task(
            label,
            total=total,
            completed=int(round(current)),
        )
        self._progress.refresh()


def _build_progress(console: Console, *, transient: bool) -> Progress:
    """Return the standard Quantas live progress layout."""
    return Progress(
        SpinnerColumn(style="green", finished_text=" "),
        TextColumn("{task.description}", style="none"),
        BarColumn(
            bar_width=32,
            style="green",
            complete_style="green",
            finished_style="green",
            pulse_style="green",
        ),
        TaskProgressColumn(
            text_format="{task.percentage:>3.0f}%",
            style="none",
            markup=False,
        ),
        TextColumn(
            "{task.completed:.0f}/{task.total:.0f}",
            style="none",
            markup=False,
        ),
        _PlainTimeRemainingColumn(),
        console=console,
        auto_refresh=False,
        transient=transient,
        redirect_stdout=False,
        redirect_stderr=False,
        expand=False,
    )


def _progress_values(event: Event) -> tuple[str, float, float]:
    """Return normalized label, current, and total values for an event."""
    data: dict[str, Any] = event.data
    label = str(data.get("label") or event.message or "Working")
    raw_total = data.get("total")
    raw_current = data.get("current")

    try:
        total = 100.0 if raw_total is None else float(raw_total)
    except (TypeError, ValueError):
        total = 100.0
    try:
        current = 0.0 if raw_current is None else float(raw_current)
    except (TypeError, ValueError):
        fraction = 0.0 if event.progress is None else float(event.progress)
        current = fraction * total
    return label, current, total


def _supports_live_progress(console: Console) -> bool:
    """Return whether a console can safely host a live progress display."""
    if os.environ.get("QUANTAS_NO_ANSI") or os.environ.get("NO_COLOR"):
        return False
    if os.environ.get("TERM", "").lower() == "dumb":
        return False
    return bool(console.is_terminal)


def _same_total(left: float, right: float) -> bool:
    """Return whether two progress totals represent the same task extent."""
    tolerance = max(abs(left), abs(right), 1.0) * 1.0e-12
    return abs(left - right) <= tolerance


__all__ = ["RichProgressDisplay"]
