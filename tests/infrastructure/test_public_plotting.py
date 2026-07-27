"""Contracts for :mod:`quantas.api.plotting`."""

from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

from quantas.api import common, plotting


def test_public_plotting_aliases_authoritative_contracts() -> None:
    """The public namespace reuses existing contracts without wrappers."""
    from quantas.models import (
        ContourPlotSpec as InternalContourPlotSpec,
        LinePlotSpec as InternalLinePlotSpec,
        PlotCollection as InternalPlotCollection,
        PolarPlotSpec as InternalPolarPlotSpec,
        SurfacePlotSpec as InternalSurfacePlotSpec,
    )

    assert plotting.LinePlotSpec is InternalLinePlotSpec
    assert plotting.ContourPlotSpec is InternalContourPlotSpec
    assert plotting.PolarPlotSpec is InternalPolarPlotSpec
    assert plotting.SurfacePlotSpec is InternalSurfacePlotSpec
    assert plotting.PlotCollection is InternalPlotCollection
    assert common.PlotCollection is plotting.PlotCollection


def test_public_plotting_supports_typed_renderer_dispatch() -> None:
    """A frontend can inspect concrete specifications through public types."""
    specification = plotting.LinePlotSpec(
        key="test-curve",
        title="Test curve",
        filename_stem="test_curve",
        x_axis=plotting.PlotAxis(key="temperature", label="Temperature", unit="K"),
        y_axis=plotting.PlotAxis(key="volume", label="Volume", unit="angstrom^3"),
        series=[
            plotting.PlotSeries(
                key="volume",
                label="Volume",
                x=np.asarray([0.0, 300.0], dtype=np.float64),
                y=np.asarray([10.0, 10.2], dtype=np.float64),
            )
        ],
    )
    collection = plotting.PlotCollection(plots=[specification])

    assert isinstance(collection.plots[0], plotting.LinePlotSpec)
    assert collection.plots[0].x_axis.unit == "K"
    np.testing.assert_allclose(collection.plots[0].series[0].y, [10.0, 10.2])


def test_public_plotting_import_is_renderer_independent(tmp_path: Path) -> None:
    """Importing neutral contracts must not load concrete frontend stacks."""
    script = tmp_path / "check_plotting_import.py"
    script.write_text(
        """
import json
import sys
from quantas.api import plotting
print(json.dumps({
    name: name in sys.modules
    for name in ("matplotlib", "plotly", "dash", "rich")
}))
""",
        encoding="utf-8",
    )
    env = os.environ.copy()
    source_root = str(Path(__file__).resolve().parents[2] / "src")
    env["PYTHONPATH"] = os.pathsep.join(
        item for item in (source_root, env.get("PYTHONPATH", "")) if item
    )
    completed = subprocess.run(
        [sys.executable, str(script)],
        env=env,
        check=True,
        capture_output=True,
        text=True,
    )

    assert json.loads(completed.stdout) == {
        "matplotlib": False,
        "plotly": False,
        "dash": False,
        "rich": False,
    }
