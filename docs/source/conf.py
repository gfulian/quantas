"""Sphinx configuration for the Quantas documentation."""

from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from quantas import __version__  # noqa: E402

project = "Quantas"
author = "Gianfranco Ulian and Giovanni Valdrè"
copyright = "2020–2026, Gianfranco Ulian and Giovanni Valdrè"
version = __version__
release = __version__

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.mathjax",
    "numpydoc",
    "sphinx_click.ext",
]

autosummary_generate = True
autodoc_typehints = "description"
autodoc_member_order = "bysource"
autoclass_content = "class"
numpydoc_show_class_members = False
numpydoc_class_members_toctree = False

source_suffix = {".rst": "restructuredtext"}
master_doc = "index"
language = "en"
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
templates_path = ["_templates"]

html_theme = "sphinx_rtd_theme"
html_title = f"Quantas {release} documentation"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_logo = "_static/branding/Quantas-logo-symbol-128.png"
html_favicon = "_static/branding/Quantas-favicon-32.png"
html_theme_options = {
    "collapse_navigation": False,
    "navigation_depth": 4,
    "sticky_navigation": False,
    "logo_only": False,
}

htmlhelp_basename = "Quantasdoc"


def _prepare_generated_assets(app: object) -> None:
    """
    Generate reproducible tutorial assets before Sphinx reads sources.
    This function should be called only when explicitly requested.
    """
    del app

    if os.environ.get("QUANTAS_REGENERATE_DOC_ASSETS") != "1":
        return

    for script in (
        "generate_elasticity_seismic_assets.py",
        "generate_thermoelasticity_assets.py",
    ):
        subprocess.run(
            [sys.executable, str(ROOT / "docs" / "tools" / script)],
            cwd=ROOT,
            check=True,
        )


def setup(app: object) -> dict[str, bool]:
    """Register optional documentation asset preparation/regeneration."""
    app.connect("builder-inited", _prepare_generated_assets)
    return {"parallel_read_safe": True, "parallel_write_safe": True}

