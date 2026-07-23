"""Tests for the supported rendering bridge in :mod:`quantas.api`."""

from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys

from quantas.api import common, elasticity, rendering


def test_public_plain_text_renderer_accepts_neutral_tables() -> None:
    """Report tables can be rendered without importing an internal backend."""
    table = common.ReportTable(
        title="Public table",
        columns=["Property", "Value"],
        rows=[["Bulk modulus", 72.5]],
        metadata={"column_units": ["", "GPa"]},
    )

    assert rendering.render_table(table) == rendering.render_tables([table])
    assert "Bulk modulus" in rendering.render_table(table)


def test_public_rendering_import_does_not_load_matplotlib(tmp_path: Path) -> None:
    """Plain-text rendering must remain usable without the optional plot stack."""
    script = tmp_path / "check_rendering_import.py"
    script.write_text(
        """
import json
import sys
from quantas.api import rendering
print(json.dumps({"matplotlib": "matplotlib" in sys.modules}))
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

    assert json.loads(completed.stdout) == {"matplotlib": False}


def test_public_plot_renderer_returns_figures_and_paths(tmp_path: Path) -> None:
    """Neutral plot collections are rendered through the supported API bridge."""
    source = Path(__file__).resolve().parents[2] / "examples/elasticity/calcite.dat"
    result = elasticity.run(
        source,
        options=elasticity.Options(calculate_2d=True, ntheta=37),
    )
    collection = elasticity.build_2d_plots(result, properties=("young",))

    rendered = rendering.render_plots(
        collection,
        output_dir=tmp_path,
        preset="publication",
        image_format="png",
        close=True,
    )

    assert len(rendered.artifacts) == 1
    assert len(rendered.figures) == 1
    assert len(rendered.paths) == 1
    assert all(path.is_file() for path in rendered.paths)
    assert all(artifact.kind == "polar" for artifact in rendered.artifacts)
