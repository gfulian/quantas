# Quantas

<p align="center"><img src="docs/source/_static/branding/Quantas_banner-1600.png" alt="Quantas — quantitative analysis of solid-state properties" width="900" ></p>

[![Quantas CI](https://github.com/gfulian/quantas/actions/workflows/ci.yml/badge.svg?branch=dev%2Frefactor)](https://github.com/gfulian/quantas/actions/workflows/ci.yml)

Quantas is a typed Python library for quantitative analysis of solid-state
thermodynamics, elasticity, acoustic-wave propagation, and equations of state.
The same scientific workflows are available through a public Python API and a
Click-based command-line interface. Frontend-neutral results, reports, plots,
and events provide the contract for future graphical applications.

The numerical core is independent of Click, Rich, Matplotlib and GUI
frameworks. Native workflow results are stored in HDF5 using `float64` and
`complex128` numerical values.

## Scientific modules

| Module | Main purpose |
| --- | --- |
| Elasticity | Tensor validation, stability, VRH averages, anisotropy, directional properties, and rotations |
| SEISMIC | Christoffel analysis, phase and group velocities, polarization, degeneracies, and caustics |
| HA | Harmonic vibrational thermodynamics |
| QHA | Quasi-harmonic equilibrium properties over pressure and temperature |
| EOS | PV, VT, and PVT fitting, uncertainty treatment, diagnostics, and derived properties |
| Thermoelasticity | Quasi-static elastic tensors over pressure, temperature, grids, and geological profiles |

## Requirements and installation

Quantas requires Python 3.10 or newer. Install the library from the repository
with:

```bash
python -m pip install .
```

Install optional plotting support or the complete development environment with:

```bash
python -m pip install ".[plot]"
python -m pip install ".[dev]"
```

ODRPACK is a runtime dependency because orthogonal-distance regression is part
of the supported EOS workflow. The future Quantas GUI is developed separately;
GUI frameworks are not dependencies of the scientific library.

## Python API

The supported public interface is organized by scientific domain under
`quantas.api`:

```python
from quantas.api import qha
from quantas.renderers.tables import render_tables

input_data = qha.read_input("material.yaml")
options = qha.Options(
    temperature_min=0.0,
    temperature_max=1000.0,
    temperature_step=10.0,
    pressure_min=0.0,
    pressure_max=10.0,
    pressure_step=1.0,
)

result = qha.run(input_data, options=options)
qha.write_result(result, "material_qha.hdf5")
print(render_tables(qha.build_report(result)))
```

The public namespaces are:

```python
from quantas.api import elasticity, eos, ha, qha, seismic, thermoelasticity
from quantas.api import interop, profiles, registry
```

`registry` exposes frontend-neutral module capabilities. `interop` contains
validated transformations between workflows, such as QHA to thermoelasticity
and thermoelasticity to SEISMIC.

## Command-line interface

List the top-level commands and module-specific help with:

```bash
quantas --help
quantas elasticity --help
quantas seismic --help
quantas ha --help
quantas qha --help
quantas eos --help
quantas thermoelasticity --help
```

Typical workflows separate calculation, inspection, plotting, and export:

```bash
quantas elasticity run calcite.dat -o calcite_elasticity.hdf5
quantas elasticity plot calcite_elasticity.hdf5 --preset publication

quantas qha run material.yaml -o material_qha.hdf5
quantas qha plot material_qha.hdf5 --preset publication

quantas thermoelasticity run material_soec.yaml material_qha.hdf5 \
  -o material_thermoelastic.hdf5
quantas thermoelasticity analysis profile material_thermoelastic.hdf5 \
  --preset mantle-katsura-2022
```

Significant calculations create a plain-text `.log` report automatically.
`--quiet` controls terminal output without suppressing the report. Reporting
uses the common levels `standard`, `extended`, and `debug`.

Plotting commands share the presets `screen`, `publication`, and `monochrome`.
Explicit figure dimensions and DPI values override preset defaults.

## Results and reproducibility

Native HDF5 results contain:

- program, version, module, method, and schema metadata;
- normalized input data and calculation options;
- numerical results and uncertainties;
- warnings and persistent workflow events;
- numerical precision and scientific provenance.

Readers and exporters inspect metadata rather than infer result types from file
names. Operational progress events are observer-only and are not persisted.

Curated real-data examples are distributed in `examples/` with provenance and
SHA-256 manifests. They support end-to-end scientific regression tests without
replacing small unit tests that isolate formulas and error conditions.

## Development and validation

Run the complete isolated test matrix with:

```bash
python tools/run_tests.py all -- -q
```

Useful developer checks include:

```bash
python tools/check_architecture.py
python tools/check_distribution.py
python tools/project_stats.py
```

The project uses Ruff, mypy, pytest, Sphinx, wheel and source-distribution
installation checks. Contributions should preserve API/CLI equivalence,
frontend independence, native HDF5 contracts, and numerical behavior.

See [CONTRIBUTING.md](CONTRIBUTING.md), [RELEASE.md](RELEASE.md), and
[ROADMAP.md](ROADMAP.md) for development and release policies.

## Documentation

The full manual is built with Sphinx:

```bash
python -m pip install ".[docs]"
sphinx-build -E -a -b html -W --keep-going docs/source docs/_build/html
```

The hosted documentation is available at
<https://quantas.readthedocs.io/>.

## Citation and license

Quantas is academic, open-source software. Please cite Quantas and the relevant
scientific methods when publishing results derived from the package. Canonical
software and method references are provided through `CITATION.cff` and the
internal citation registry.

Quantas is distributed under the BSD 3-Clause license. See [LICENSE](LICENSE).
