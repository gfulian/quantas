"""Repository-aware validation tool behavior."""

from __future__ import annotations

from pathlib import Path
import subprocess
import sys


def test_repository_checker_skips_plain_source_directories(tmp_path: Path) -> None:
    """Git-only validation is optional for extracted public source archives."""
    tool = Path(__file__).resolve().parents[2] / "tools" / "check_repository.py"
    copied = tmp_path / "check_repository.py"
    copied.write_text(tool.read_text(encoding="utf-8"), encoding="utf-8")
    completed = subprocess.run(
        [sys.executable, str(copied)],
        cwd=tmp_path,
        check=False,
        text=True,
        capture_output=True,
    )
    assert completed.returncode == 0
    assert "skipped" in completed.stdout.lower()
