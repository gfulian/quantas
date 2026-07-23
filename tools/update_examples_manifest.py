#!/usr/bin/env python3
"""Generate or verify the deterministic manifest of distributed examples."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[1]
EXAMPLES_ROOT = PROJECT_ROOT / "examples"
MANIFEST_PATH = EXAMPLES_ROOT / "manifest.json"
CHECKSUM_PATH = EXAMPLES_ROOT / "MANIFEST.sha256"
_IGNORED_FILES = {MANIFEST_PATH.name, CHECKSUM_PATH.name}
_IGNORED_DIRECTORIES = {"__pycache__", ".pytest_cache", ".mypy_cache", ".ruff_cache"}
_IGNORED_SUFFIXES = {".pyc", ".pyo"}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def build_manifest() -> dict[str, object]:
    """Return deterministic metadata for every curated example file."""
    files: list[dict[str, object]] = []
    for path in sorted(EXAMPLES_ROOT.rglob("*")):
        if not path.is_file():
            continue
        relative_path = path.relative_to(EXAMPLES_ROOT)
        if path.name in _IGNORED_FILES:
            continue
        if any(part in _IGNORED_DIRECTORIES for part in relative_path.parts):
            continue
        if path.suffix in _IGNORED_SUFFIXES:
            continue
        relative = relative_path.as_posix()
        files.append(
            {
                "path": relative,
                "size_bytes": path.stat().st_size,
                "sha256": _sha256(path),
            }
        )
    return {
        "schema": "quantas-examples-manifest",
        "version": 1,
        "files": files,
    }


def _serialized(manifest: dict[str, object]) -> bytes:
    return (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode("utf-8")


def write_manifest() -> None:
    """Write the JSON manifest and its detached SHA-256 record."""
    payload = _serialized(build_manifest())
    MANIFEST_PATH.write_bytes(payload)
    digest = hashlib.sha256(payload).hexdigest()
    CHECKSUM_PATH.write_text(f"{digest}  {MANIFEST_PATH.name}\n", encoding="utf-8")


def check_manifest() -> bool:
    """Return whether checked-in manifest files match the current examples."""
    if not MANIFEST_PATH.is_file() or not CHECKSUM_PATH.is_file():
        return False
    payload = _serialized(build_manifest())
    expected = MANIFEST_PATH.read_bytes()
    if payload != expected:
        return False
    digest = hashlib.sha256(expected).hexdigest()
    return CHECKSUM_PATH.read_text(encoding="utf-8") == (
        f"{digest}  {MANIFEST_PATH.name}\n"
    )


def main() -> int:
    """Run the manifest generator command-line interface."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="verify the checked-in manifest instead of rewriting it",
    )
    args = parser.parse_args()
    if args.check:
        if check_manifest():
            print("Examples manifest is current.")
            return 0
        print("Examples manifest is missing or stale.", file=sys.stderr)
        return 1
    write_manifest()
    print(f"Wrote {MANIFEST_PATH.relative_to(PROJECT_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
