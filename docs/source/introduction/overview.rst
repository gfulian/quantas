Overview
========

Quantas is a scientific software project for the quantitative analysis of
solid-state properties from theoretical or experimental data.  The current
Quantas 2.0 series turns the historical program into a reusable Python library
with equivalent access from the command line and, in the future, from a
separate graphical interface.

Scientific domains currently covered by Quantas include:

- **Harmonic thermodynamics (HA)** from phonon frequencies and q-point
  weights;
- **Quasi-harmonic thermodynamics (QHA)**, including equilibrium volume,
  expansion, and derived thermoelastic properties;
- **Elasticity**, including tensor analysis, stability, averaging, directional
  properties, and 2D/3D visualization;
- **SEISMIC**, based on the Christoffel equation and dedicated to phase/group
  velocities, polarization, and derived seismic observables;
- **Equations of state (EOS)** in the P--V, V--T, and P--V--T domains;
- **Thermoelasticity**, with quasi-static pressure-temperature evolution of the
  elastic tensor and interoperability with SEISMIC.

Design principles
-----------------

Quantas 2.0 follows a small number of explicit architectural principles.

1. **Scientific reliability before cleanup or optimization.**
   Numerical conventions, units, and historical behavior are preserved unless a
   correction is scientifically justified and validated.
2. **One scientific workflow, multiple frontends.**
   API, CLI, and Quantas GUI must agree numerically and consume the same
   passive contracts whenever possible.
3. **Frontend-neutral core.**
   The numerical core must not depend on Click, Rich, Dash, Matplotlib,
   Plotly, terminal streams, or progress bars.
4. **Native HDF5 persistence.**
   Quantas stores normalized inputs, options, results, warnings, and
   methodological provenance in a consistent result envelope.
5. **Stable public Python API.**
   The only supported Python surface is :mod:`quantas.api`.

How to read this manual
-----------------------

A recommended reading order is:

#. :doc:`capabilities`
#. :doc:`../getting_started/installation`
#. :doc:`../getting_started/index`
#. :doc:`../workflows/concepts/index`
#. :doc:`../theory/index`
#. :doc:`../workflows/index`
#. :doc:`../tutorials/index`

This structure intentionally separates:

- **scientific background**;
- **workflow logic**;
- **tutorial execution**;
- **normative API and CLI reference**.

That separation is essential to keep the manual readable for researchers,
command-line users, Python users, and future contributors.
