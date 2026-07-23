# -*- coding: utf-8 -*-

"""Text readers for geological thermoelastic pressure-temperature profiles."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from quantas.core.physics.earth_profiles import read_earth_profile_spec
from quantas.modules.thermoelasticity.models import ThermoelasticDepthProfile


def read_thermoelastic_depth_profile(
    filename: str | Path,
    *,
    name: str | None = None,
) -> ThermoelasticDepthProfile:
    """Read a depth-pressure-temperature profile from CSV or whitespace text.

    Parameters
    ----------
    filename : str or Path
        Input file. A header containing ``depth_km``, ``P_GPa`` and ``T_K`` is
        required. Comma-separated and whitespace-separated files are accepted.
    name : str or None, optional
        Profile identifier. The filename stem is used by default.

    Returns
    -------
    ThermoelasticDepthProfile
        Sorted and validated profile.

    Raises
    ------
    ValueError
        If required columns are absent, values are invalid, or depths repeat.
    """
    path = Path(filename)
    text = path.read_text(encoding="utf-8")
    lines = [
        line
        for line in text.splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    if not lines:
        raise ValueError("depth-profile file contains no data")
    if "," in lines[0]:
        rows = list(csv.DictReader(lines))
    else:
        headers = lines[0].split()
        rows = [dict(zip(headers, line.split(), strict=True)) for line in lines[1:]]
    aliases = {
        "depth": ("depth_km", "depth", "z_km"),
        "pressure": ("P_GPa", "pressure_GPa", "pressure", "P"),
        "temperature": ("T_K", "temperature_K", "temperature", "T"),
    }
    if not rows:
        raise ValueError("depth-profile file contains no numerical rows")

    def column(kind: str) -> np.ndarray:
        key = next(
            (candidate for candidate in aliases[kind] if candidate in rows[0]), None
        )
        if key is None:
            raise ValueError(f"missing required {kind} column")
        try:
            return np.asarray([float(row[key]) for row in rows], dtype=np.float64)
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"invalid values in {kind} column") from exc

    depth = column("depth")
    order = np.argsort(depth, kind="stable")
    return ThermoelasticDepthProfile(
        name=path.stem if name is None else name,
        depth=depth[order],
        pressure=column("pressure")[order],
        temperature=column("temperature")[order],
        metadata={"kind": "tabulated", "source": str(path)},
    )


def read_thermoelastic_profile_spec(
    filename: str | Path,
) -> ThermoelasticDepthProfile:
    """Read a composed geothermobarometric profile from YAML.

    Parameters
    ----------
    filename : str or Path
        Earth-profile YAML specification.  Pressure and temperature models are
        built independently by :mod:`quantas.core.physics.earth_profiles`;
        relative tabulated paths are resolved from the YAML location.

    Returns
    -------
    ThermoelasticDepthProfile
        Evaluated profile adapted to the thermoelastic workflow contract.

    Raises
    ------
    ValueError
        If the specification or any referenced model is invalid.
    """
    profile = read_earth_profile_spec(filename)
    return ThermoelasticDepthProfile(
        name=profile.name,
        depth=profile.depth,
        pressure=profile.pressure,
        temperature=profile.temperature,
        metadata=profile.metadata,
    )


def write_thermoelastic_profile_template(
    filename: str | Path,
    *,
    overwrite: bool = False,
) -> Path:
    """Write an editable Earth-profile YAML specification template.

    Parameters
    ----------
    filename : str or Path
        Destination YAML file.
    overwrite : bool, optional
        Replace an existing file.

    Returns
    -------
    Path
        Written template path.

    Notes
    -----
    The template combines PREM pressure with a two-part temperature model: a
    user-editable continental conductive lithosphere and the Katsura (2022)
    mantle adiabat.  All values are examples and should be reviewed for the
    intended tectonic setting.
    """
    path = Path(filename).with_suffix(".yaml")
    if path.exists() and not overwrite:
        raise FileExistsError(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    text = """# Quantas composed geothermobarometric profile specification
# Review every lithospheric parameter before quantitative use.
schema_version: 1
name: custom-earth-profile

depth:
  min_km: 0.0
  max_km: 660.0
  step_km: 1.0
  include_critical_depths: true

pressure:
  model: prem
  max_depth_km: 2891.0
  integration_step_km: 0.25

temperature:
  model: piecewise
  name: custom-lithosphere-mantle-temperature
  segments:
    - depth_min_km: 0.0
      depth_max_km: 100.0
      model: continental-conductive
      name: editable-continental-lithosphere
      surface_temperature_K: 288.15
      surface_heat_flow_mW_m2: 40.0
      layers:
        - name: upper-crust
          thickness_km: 15.0
          conductivity_W_mK: 2.5
          heat_production_uW_m3: 1.0
        - name: lower-crust
          thickness_km: 20.0
          conductivity_W_mK: 2.5
          heat_production_uW_m3: 0.3
        - name: lithospheric-mantle
          thickness_km: 65.0
          conductivity_W_mK: 3.3
          heat_production_uW_m3: 0.02
    - depth_min_km: 100.0
      depth_max_km: 660.0
      model: katsura-2022
      join:
        mode: continuous-offset
"""
    path.write_text(text, encoding="utf-8")
    return path


__all__ = [
    "read_thermoelastic_depth_profile",
    "read_thermoelastic_profile_spec",
    "write_thermoelastic_profile_template",
]
