#!/usr/bin/env python3
"""Run deterministic Quantas pytest suites in isolated subprocesses."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from typing import Sequence


PROJECT_ROOT = Path(__file__).resolve().parents[1]
LIBRARY_STAGES: tuple[str, ...] = (
    "core",
    "elasticity",
    "seismic",
    "ha",
    "qha",
    "eos",
    "thermoelasticity",
)

PLOTTING_STAGES: tuple[str, ...] = (
    "elasticity-plotting",
    "seismic-plotting",
    "ha-plotting",
    "qha-plotting",
    "eos-plotting",
    "thermoelasticity-plotting",
)

STAGED_ALL_SUITES: tuple[str, ...] = (
    *LIBRARY_STAGES,
    "cli",
    *PLOTTING_STAGES,
    "examples",
)

SUITE_PATHS: dict[str, tuple[str, ...]] = {
    "core": (
        "tests/infrastructure",
        "tests/physics/units",
        "tests/physics/eos",
        "tests/physics/thermodynamics",
        "tests/physics/earth_profiles",
        "tests/math",
        "tests/chemistry",
        "tests/interfaces",
    ),
    "elasticity": (
        "tests/physics/elasticity",
        "tests/modules/elasticity",
    ),
    "seismic": (
        "tests/physics/seismic",
        "tests/modules/seismic",
    ),
    "ha": ("tests/modules/ha",),
    "qha": ("tests/modules/qha",),
    "eos": ("tests/modules/eos",),
    "thermoelasticity": ("tests/modules/thermoelasticity",),
    "examples": ("tests/examples",),
    "cli": ("tests/infrastructure", "tests/modules"),
    "elasticity-plotting": ("tests/modules/elasticity",),
    "seismic-plotting": ("tests/modules/seismic",),
    "ha-plotting": ("tests/modules/ha",),
    "qha-plotting": ("tests/modules/qha",),
    "eos-plotting": ("tests/modules/eos",),
    "thermoelasticity-plotting": ("tests/modules/thermoelasticity",),
}


def _parser() -> argparse.ArgumentParser:
    """Return the command-line parser for the test-suite runner."""
    parser = argparse.ArgumentParser(
        description=(
            "Run one named Quantas pytest suite. The special 'all' target "
            "runs domain-scoped library, CLI, and plotting suites in isolated subprocesses."
        )
    )
    parser.add_argument(
        "suite",
        nargs="?",
        default="all",
        help="Suite name accepted by pytest --suite (default: all).",
    )
    parser.add_argument(
        "-k",
        dest="keyword",
        help="Forward a pytest -k expression to every selected stage.",
    )
    parser.add_argument("-v", "--verbose", action="count", default=0)
    parser.add_argument(
        "--allow-plugin-autoload",
        action="store_true",
        help=(
            "Allow third-party pytest plugin autoloading. It is disabled by "
            "default to keep collection and shutdown deterministic."
        ),
    )
    parser.add_argument(
        "--fail-fast",
        action="store_true",
        help="Stop the staged complete run after the first failed suite.",
    )
    parser.add_argument(
        "--stage-timeout",
        type=float,
        default=600.0,
        metavar="SECONDS",
        help=(
            "Maximum runtime for each isolated pytest stage. Use 0 to "
            "disable the timeout (default: 600 seconds)."
        ),
    )
    return parser


def _pytest_command(
    suite: str,
    *,
    verbosity: int,
    keyword: str | None,
    extra: Sequence[str],
) -> list[str]:
    """Build one pytest subprocess command.

    Parameters
    ----------
    suite
        Quantas suite passed to ``pytest --suite``.
    verbosity
        Number of requested ``-v`` flags.
    keyword
        Optional pytest ``-k`` expression.
    extra
        Additional pytest arguments.

    Returns
    -------
    list of str
        Complete subprocess command.
    """
    command = [sys.executable, "-m", "pytest", "--suite", suite]
    command.extend(SUITE_PATHS.get(suite, ()))
    if verbosity:
        command.append("-" + "v" * verbosity)
    if keyword:
        command.extend(["-k", keyword])
    command.extend(arg for arg in extra if arg != "--")
    return command


def _selected_suites(suite: str) -> tuple[str, ...]:
    """Return the subprocess stages associated with one public target."""
    if suite == "all":
        return STAGED_ALL_SUITES
    if suite == "library":
        return LIBRARY_STAGES
    return (suite,)


def _test_environment(*, allow_plugin_autoload: bool) -> dict[str, str]:
    """Return the deterministic environment used by pytest subprocesses."""
    environment = os.environ.copy()
    environment.setdefault("MPLBACKEND", "Agg")
    for variable in (
        "OPENBLAS_NUM_THREADS",
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        environment[variable] = "1"
    if not allow_plugin_autoload:
        environment["PYTEST_DISABLE_PLUGIN_AUTOLOAD"] = "1"
    return environment


def _spawn_stage(
    command: Sequence[str],
    *,
    environment: dict[str, str],
    timeout: float | None,
) -> int:
    """Run one pytest stage with the platform's most reliable primitive."""
    arguments = [str(part) for part in command]
    timeout_executable = shutil.which("timeout")
    if sys.platform.startswith("linux") and timeout_executable is not None:
        if timeout is not None:
            arguments = [
                timeout_executable,
                "--foreground",
                "--signal=TERM",
                "--kill-after=2",
                f"{timeout:.17g}",
                *arguments,
            ]
        return int(
            os.spawnve(
                os.P_WAIT,
                arguments[0],
                arguments,
                environment,
            )
        )

    try:
        completed = subprocess.run(
            arguments,
            cwd=PROJECT_ROOT,
            env=environment,
            check=False,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        return 124
    return int(completed.returncode)


def _run_command(
    command: Sequence[str],
    *,
    environment: dict[str, str],
    timeout: float | None,
) -> int:
    """Run one pytest command and return its exit status.

    Returns status 124 when the isolated stage exceeds its configured timeout.
    """
    print("\n$ " + " ".join(command), flush=True)
    return _spawn_stage(
        command,
        environment=environment,
        timeout=timeout,
    )


def main(argv: Sequence[str] | None = None) -> int:
    """Run the requested suite or staged complete validation.

    Parameters
    ----------
    argv
        Optional command-line arguments used mainly by tests.

    Returns
    -------
    int
        Zero when every selected stage passes, otherwise the first non-zero
        pytest exit status.
    """
    args, pytest_args = _parser().parse_known_args(argv)
    environment = _test_environment(
        allow_plugin_autoload=args.allow_plugin_autoload,
    )
    statuses: list[tuple[str, int]] = []

    try:
        with tempfile.TemporaryDirectory(
            prefix="quantas-matplotlib-",
        ) as matplotlib_root:
            for stage_index, suite in enumerate(_selected_suites(args.suite), start=1):
                stage_environment = environment.copy()
                matplotlib_config = Path(matplotlib_root) / f"{stage_index:02d}-{suite}"
                matplotlib_config.mkdir(parents=True, exist_ok=True)
                stage_environment["MPLCONFIGDIR"] = str(matplotlib_config)
                command = _pytest_command(
                    suite,
                    verbosity=args.verbose,
                    keyword=args.keyword,
                    extra=tuple(pytest_args),
                )
                timeout = None if args.stage_timeout <= 0.0 else args.stage_timeout
                status = _run_command(
                    command,
                    environment=stage_environment,
                    timeout=timeout,
                )
                statuses.append((suite, status))
                if status != 0 and args.fail_fast:
                    break
    except KeyboardInterrupt:
        print("\nQuantas test run interrupted.", file=sys.stderr)
        return 130

    if len(statuses) > 1:
        print("\nQuantas staged test summary:")
        for suite, status in statuses:
            label = "passed" if status == 0 else f"failed ({status})"
            print(f"  {suite:<9} {label}")

    return next((status for _, status in statuses if status != 0), 0)


if __name__ == "__main__":
    raise SystemExit(main())
