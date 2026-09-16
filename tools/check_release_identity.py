#!/usr/bin/env python3
"""Validate the Quantas release identity across source metadata and an optional tag."""

from __future__ import annotations

import argparse
import runpy
from pathlib import Path
import re
from typing import Sequence


PROJECT_ROOT = Path(__file__).resolve().parents[1]
VERSION_MODULE = PROJECT_ROOT / "src" / "quantas" / "_version.py"


def _parser() -> argparse.ArgumentParser:
    """Return the command-line parser."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--tag",
        help="Optional release tag; when supplied it must equal v<package-version>.",
    )
    return parser


def _require(pattern: str, text: str, *, source: str) -> None:
    """Raise ``RuntimeError`` when *pattern* is absent from *text*."""
    if re.search(pattern, text, flags=re.MULTILINE) is None:
        raise RuntimeError(f"Release identity mismatch in {source}: {pattern!r}")


def validate_release_identity(*, tag: str | None = None) -> str:
    """Validate release-facing version metadata and return the package version.

    Parameters
    ----------
    tag : str, optional
        Git release tag to validate. When provided, it must equal ``v`` followed
        by the authoritative package version.

    Returns
    -------
    str
        Authoritative Quantas package version.

    Raises
    ------
    RuntimeError
        If release-facing source metadata or *tag* disagrees with the
        authoritative package version.
    """
    namespace = runpy.run_path(str(VERSION_MODULE))
    version = str(namespace["__version__"])

    citation = (PROJECT_ROOT / "CITATION.cff").read_text(encoding="utf-8")
    project_state = (PROJECT_ROOT / "PROJECT_STATE.md").read_text(encoding="utf-8")
    roadmap = (PROJECT_ROOT / "ROADMAP.md").read_text(encoding="utf-8")
    changelog = (PROJECT_ROOT / "CHANGELOG.md").read_text(encoding="utf-8")

    _require(
        rf"^version:\s*{re.escape(version)}\s*$",
        citation,
        source="CITATION.cff",
    )
    _require(
        rf"^\| Current development version \| `{re.escape(version)}` \|$",
        project_state,
        source="PROJECT_STATE.md",
    )
    _require(
        rf"Development is now on ``{re.escape(version)}``",
        roadmap,
        source="ROADMAP.md",
    )
    _require(
        rf"^## \[{re.escape(version)}\] - Unreleased$",
        changelog,
        source="CHANGELOG.md",
    )

    if tag is not None:
        expected_tag = f"v{version}"
        if tag != expected_tag:
            raise RuntimeError(
                f"Release tag {tag!r} does not match package version {expected_tag!r}."
            )

    return version


def main(argv: Sequence[str] | None = None) -> int:
    """Validate the source/release identity and print the accepted version."""
    args = _parser().parse_args(argv)
    version = validate_release_identity(tag=args.tag)
    print(f"Release identity verified for Quantas {version}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
