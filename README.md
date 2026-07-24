# Quantas

<p align="center"><img src="docs/source/_static/branding/Quantas_banner-1600.png" alt="Quantas — quantitative analysis of solid-state properties" width="900" ></p>

[![Quantas CI](https://github.com/gfulian/quantas/actions/workflows/ci.yml/badge.svg?branch=dev%2Frefactor)](https://github.com/gfulian/quantas/actions/workflows/ci.yml)
[![Development Status: Beta](https://img.shields.io/badge/status-beta-orange)](https://quantas.readthedocs.io/en/2.0.0-beta/)
[![Quantas 2.0 beta documentation](https://readthedocs.org/projects/quantas/badge/?version=2.0.0-beta)](https://quantas.readthedocs.io/en/2.0.0-beta/)
[![Python 3.10+](https://img.shields.io/badge/Python-%3E%3D3.10-blue)](https://www.python.org/downloads/)
[![License: BSD 3-Clause](https://img.shields.io/badge/License-BSD--3--Clause-blue)](LICENSE)

Quantas is an open-source scientific toolkit for the quantitative analysis of
solid-state properties. It was created to make demanding thermodynamic,
elastic, seismic, and equation-of-state calculations easier to perform,
inspect, reproduce, and reuse.

This is a complete refactoring of the original program (version 0.9.1) into a
modern Python library with a command-line interface and frontend-neutral
scientific workflows. The same numerical calculation can be used from Python,
from the terminal, or by future graphical applications.

> **Quantas 2 is currently in beta.**
>
> The numerical workflows are extensively tested, but the public API, file
> schemas, and user interface may still receive changes before the first stable
> 2.0 release.

A graphical user interface is in active development, and it will be provided as
a separate package.

## Scientific capabilities

Quantas currently provides workflows for:

- harmonic and quasi-harmonic thermodynamics;
- second-order elasticity, stability, VRH averages, rotations, and directional
  properties;
- Christoffel analysis, phase and group velocities, polarization, anisotropy,
  and seismic surfaces;
- pressure-volume, volume-temperature, and pressure-volume-temperature
  equations of state;
- thermoelastic properties over pressure, temperature, and geological profiles,
  within the quasi-static approximation (QSA).

Results are stored in self-describing HDF5 archives and can be inspected,
plotted, exported, and passed between compatible workflows.

## Requirements and installation

Quantas requires Python 3.10 or newer.

Clone the repository and install the current beta with:

```bash
git clone --branch dev/refactor https://github.com/gfulian/quantas.git
cd quantas
python -m pip install .
```

Install optional plotting support or the complete development environment with:

```bash
python -m pip install ".[plot]"
python -m pip install ".[dev]"
```

## Documentation

The complete Quantas 2 beta manual includes installation instructions,
tutorials, command-line reference, Python API reference, interoperability
examples, scientific methods, and validation material:

[Read the Quantas 2.0 beta documentation](https://quantas.readthedocs.io/en/2.0.0-beta/)

Documentation for the historical Quantas 0.9.1 release remains available
through the Read the Docs version selector.

## Citation

Quantas is academic, open-source software. Please cite Quantas and the relevant
scientific methods when publishing results derived from the package. Canonical
software and method references are provided through `CITATION.cff` and the
internal citation registry.



## Contributing and license

Contributions, bug reports, validation datasets, and scientific comparisons
are welcome. See [How to contribute](CONTRIBUTING.md) and the
[expected development roadmap](ROADMAP.md).

Quantas is distributed under the BSD 3-Clause license. See [LICENSE](LICENSE).
