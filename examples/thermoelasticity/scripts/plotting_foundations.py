#!/usr/bin/env python
"""Build thermoelastic fit, P-T, profile, and domain plots through the API.

Examples
--------
Render selected plot families from a thermoelastic HDF5 archive::

    python plotting_foundations.py dolomite_profile.hdf5 --output plots
"""

from __future__ import annotations

import argparse
from pathlib import Path

from quantas.api import rendering, thermoelasticity


def parser() -> argparse.ArgumentParser:
    """Return the command-line parser for this standalone API example."""
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("archive", type=Path)
    result.add_argument("--output", type=Path, default=Path("thermoelastic_plots"))
    result.add_argument("--profile-name")
    result.add_argument(
        "--components",
        nargs="+",
        default=["C11", "C33", "C44"],
    )
    return result


def main() -> None:
    """Read one archive and render the foundational plot families."""
    args = parser().parse_args()
    result = thermoelasticity.read_result(args.archive)
    args.output.mkdir(parents=True, exist_ok=True)

    collections = [
        thermoelasticity.build_fit_plots(result, components=args.components),
        thermoelasticity.build_pt_plots(result, components=args.components[:2]),
        thermoelasticity.build_domain_plot(result),
    ]
    payload = thermoelasticity.get_result(result)
    if payload.profiles:
        profile_name = args.profile_name or next(iter(payload.profiles))
        collections.append(
            thermoelasticity.build_profile_plots(
                result,
                profile_name=profile_name,
                components=args.components,
                options=thermoelasticity.ProfilePlotOptions(
                    mode="relative",
                    uncertainty="auto",
                    color_by="temperature",
                ),
            )
        )

    for collection in collections:
        rendered = rendering.render_plots(
            collection,
            output_dir=args.output,
            image_format="png",
            close=True,
        )
        for warning in rendered.warnings:
            print(f"WARNING: {warning}")
        for artifact in rendered.artifacts:
            if artifact.path is not None:
                print(artifact.path)


if __name__ == "__main__":
    main()
