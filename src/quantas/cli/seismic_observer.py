# -*- coding: utf-8 -*-

"""Terminal observer for seismic CLI workflow events."""

from __future__ import annotations

from quantas.cli.output import CLIOutput
from quantas.core.events import Event, EventLevel


class SeismicProgressObserver:
    """Render seismic messages and live progress through shared CLI output.

    Parameters
    ----------
    silent : bool, optional
        Suppress terminal output while retaining report collection.
    show_progress : bool, optional
        Enable the live Rich progress display.
    output : CLIOutput or None, optional
        Shared output service.
    """

    def __init__(
        self,
        *,
        silent: bool = False,
        show_progress: bool = True,
        output: CLIOutput | None = None,
    ) -> None:
        self.output = (
            CLIOutput(silent=silent, show_progress=show_progress)
            if output is None
            else output
        )
        self.show_progress = bool(show_progress)

    def __call__(self, event: Event) -> None:
        """Receive and render one Quantas event."""
        if event.level == EventLevel.PROGRESS and self.show_progress:
            self.output.progress(event)
        elif event.level in {
            EventLevel.INFO,
            EventLevel.WARNING,
            EventLevel.ERROR,
        }:
            self.output.message(event.message, level=event.level)

    def close(self) -> None:
        """Release terminal progress resources."""
        self.output.close()
