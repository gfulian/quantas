"""Regenerate thermoelasticity tutorial figures and reference tables."""

from __future__ import annotations

from pathlib import Path
import os
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[2]
STATIC_OUTPUT = ROOT / "docs" / "source" / "_static" / "tutorials" / "thermoelasticity"
DOWNLOAD_OUTPUT = ROOT / "docs" / "source" / "_downloads" / "tutorials" / "thermoelasticity"


def main() -> None:
    """Execute the public API tutorial and copy stable documentation assets."""
    STATIC_OUTPUT.mkdir(parents=True, exist_ok=True)
    DOWNLOAD_OUTPUT.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="quantas-thermoelasticity-docs-") as tmp:
        output = Path(tmp) / "tutorial"
        environment = os.environ.copy()
        source_path = str(ROOT / "src")
        environment["PYTHONPATH"] = (
            source_path
            if not environment.get("PYTHONPATH")
            else source_path + os.pathsep + environment["PYTHONPATH"]
        )
        subprocess.run(
            [
                sys.executable,
                str(ROOT / "examples" / "thermoelasticity" / "tutorial_api.py"),
                "--output-dir",
                str(output),
            ],
            cwd=ROOT,
            check=True,
            stdout=subprocess.DEVNULL,
            env=environment,
        )
        copies = {
            output / "figures" / "fit_C11_cold_finite_strain_fit.png":
                STATIC_OUTPUT / "dolomite_c11_fit.png",
            output / "figures" / "pt_thermoelastic_pt_value.png":
                STATIC_OUTPUT / "dolomite_pt_c11_c33.png",
            output / "figures" / "compare_thermoelastic_T_S_compare_P0GPa.png":
                STATIC_OUTPUT / "dolomite_isothermal_adiabatic.png",
            output / "figures" / "profile_profile_dolomite_continental_demo_C11_C33_C44_relative.png":
                STATIC_OUTPUT / "dolomite_profile_relative.png",
            output / "figures" / "domain_thermoelastic_PT_domain.png":
                STATIC_OUTPUT / "dolomite_pt_domain.png",
            output / "dolomite_grid.csv":
                DOWNLOAD_OUTPUT / "dolomite_grid.csv",
            output / "dolomite_profile.csv":
                DOWNLOAD_OUTPUT / "dolomite_profile.csv",
            output / "dolomite_5GPa_800K.dat":
                DOWNLOAD_OUTPUT / "dolomite_5GPa_800K.dat",
        }
        for source, destination in copies.items():
            if not source.is_file():
                raise RuntimeError(f"Expected tutorial asset was not generated: {source}")
            shutil.copyfile(source, destination)


if __name__ == "__main__":
    main()
