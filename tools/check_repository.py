#!/usr/bin/env python3
"""Validate Git-only release conditions when a worktree is available."""

from __future__ import annotations

from pathlib import Path
import subprocess
from typing import Sequence


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def _git(*arguments: str, capture: bool = False) -> subprocess.CompletedProcess[str]:
    """Run one Git command from the project root."""
    return subprocess.run(
        ["git", *arguments],
        cwd=PROJECT_ROOT,
        check=False,
        text=True,
        capture_output=capture,
    )


def main(argv: Sequence[str] | None = None) -> int:
    """Check whitespace and cleanliness, or skip outside a Git worktree."""
    del argv
    probe = _git("rev-parse", "--show-toplevel", capture=True)
    if probe.returncode != 0:
        print("No Git worktree detected; repository-only checks were skipped.")
        return 0
    worktree_root = Path(probe.stdout.strip()).resolve()
    if worktree_root != PROJECT_ROOT.resolve():
        print("The source tree is not the root of a Git worktree; repository-only checks were skipped.")
        return 0

    whitespace = _git("diff", "--check")
    if whitespace.returncode != 0:
        return whitespace.returncode

    status = _git("status", "--porcelain", capture=True)
    if status.returncode != 0:
        return status.returncode
    if status.stdout.strip():
        print("Git worktree is not clean:")
        print(status.stdout, end="")
        return 1

    print("Git worktree is clean.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
