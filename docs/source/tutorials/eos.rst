EOS fitting tutorials
=====================

Equation-of-state analysis is not one single fit.  Quantas distinguishes three
scientific domains that answer different questions:

- **P--V**: how a structure compresses along an isotherm;
- **V--T**: how a structure expands along an isobar;
- **P--V--T**: how pressure and temperature act together.

The tutorials below use curated experimental datasets and follow the same
sequence in every domain:

#. inspect the input and its uncertainties;
#. choose a physically defensible EOS formulation;
#. choose a regression method consistent with the uncertainty model;
#. validate the complete batch before fitting;
#. preserve every candidate in a native HDF5 archive;
#. inspect parameters, residuals, covariance, and extrapolation;
#. calculate derived properties and generate figures.

The examples deliberately compare several formulations and solvers.  The
smallest residual is not automatically the best scientific model: parameter
stability, uncertainty, correlations, physical limits, and the sampled domain
must be considered together.

Recommended order
-----------------

.. toctree::
   :maxdepth: 1

   eos_batch
   eos_pv
   eos_vt
   eos_pvt
   eos_diagnostics_calculator
   eos_plotting
   eos_mgd_fit_api

The first page explains the common batch and solver concepts.  The following
three pages apply them in P--V, V--T, and P--V--T order.
