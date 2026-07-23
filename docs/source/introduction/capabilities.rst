Capabilities
============

This page summarizes what Quantas can currently do, which approximations are
implemented, and which outputs are produced.

Scientific modules
------------------

HA
^^

The harmonic workflow computes vibrational thermodynamic properties from phonon
frequencies and q-point weights.  Quantas reports quantities such as Helmholtz
free energy, internal energy, entropy, zero-point energy, and isochoric heat
capacity.

QHA
^^^

The quasi-harmonic workflow extends the harmonic approximation by considering
how phonon frequencies or thermodynamic quantities evolve with volume.  Quantas
supports frequency-based and thermodynamic interpolation schemes, free-energy
minimization, and derived properties such as :math:`\alpha_V`, :math:`C_P`,
:math:`K_T`, and :math:`K_S`.

Elasticity
^^^^^^^^^^

The elasticity workflow analyzes second-order elastic tensors.  It covers
stability, Voigt--Reuss--Hill averages, bulk and shear moduli, compliance,
directional mechanical properties, exact directional extrema, tensor rotation,
and 2D/3D visualizations.

SEISMIC
^^^^^^^

The seismic workflow uses the elastic tensor and density to solve the
Christoffel equation and derive phase velocity, group velocity, polarization,
enhancement, and other directional observables relevant to acoustic and
geophysical analysis.

EOS
^^^

The equation-of-state workflow fits and evaluates models in the P--V, V--T,
and P--V--T domains.  Quantas supports OLS, WLS, effective variance, and ODR,
with native HDF5 archives, diagnostics, and plotting.

Thermoelasticity
^^^^^^^^^^^^^^^^

The thermoelastic workflow couples quasi-harmonic and elastic information to
construct :math:`C_{ijkl}(P,T)` in the quasi-static approximation.  Point,
Grid, and Profile analyses are supported, together with isothermal-to-
adiabatic conversion and export to SEISMIC.

Shared outputs
--------------

Across workflows, Quantas can provide:

- frontend-neutral report tables;
- frontend-neutral plot specifications;
- native HDF5 result archives with normalized inputs and provenance;
- deterministic plain-text logs;
- workflow-specific export files such as CSV tables or SEISMIC-ready inputs.

What is new in Quantas 2.0
--------------------------

Compared with the historical Quantas 0.9 series, Quantas 2.0 introduces:

- a stable public Python API under :mod:`quantas.api`;
- complete separation between numerical core and frontend adapters;
- Click/Rich CLI with grouped options and neutral reporting contracts;
- native HDF5 persistence across workflows;
- explicit interoperability between modules;
- a stronger documentation architecture oriented around theory, workflows,
  tutorials, and reference sections;
- a codebase prepared for a future GUI without changing the scientific core.

For a concise historical comparison, see :doc:`changes_from_0_9`.
