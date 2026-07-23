Quantas command line
====================

The command reference is generated directly from the public Click application.
Consequently, command names, option types, accepted choices, defaults, option
groups, and long-form help are version-matched to the installed Quantas source.
The surrounding pages add workflow context without copying the parser
specification by hand.

Start by selecting a scientific group:

.. code-block:: console

   quantas ha --help
   quantas qha --help
   quantas elasticity --help
   quantas seismic --help
   quantas eos --help
   quantas thermoelasticity --help

Most calculations follow ``run -> inspect/plot/export``.  EOS uses a persistent
record archive because one dataset commonly produces several candidate fits;
Thermoelasticity separates calibration from point, grid, and profile analysis.
Those differences are explained in the corresponding module pages.

Before using the generated option lists, read :doc:`conventions` for range
syntax, repeatable options, units, reports, overwrite behavior, plotting
presets, and exit statuses.

Command groups
--------------

.. list-table::
   :header-rows: 1
   :widths: 22 46 32

   * - Group
     - Primary purpose
     - Typical persisted result
   * - ``ha``
     - Harmonic thermodynamics at fixed sampled volumes
     - HA HDF5 archive
   * - ``qha``
     - Equilibrium properties over pressure and temperature
     - QHA HDF5 archive
   * - ``elasticity``
     - Stability, aggregate moduli, and directional elasticity
     - Elasticity HDF5 archive
   * - ``seismic``
     - Christoffel phase/group fields, polarization, and enhancement
     - SEISMIC HDF5 archive
   * - ``eos``
     - P--V, V--T, and P--V--T candidate fitting and record management
     - Persistent EOS archive
   * - ``thermoelasticity``
     - Cold finite-strain calibration and QSA evaluation
     - Fit, grid, or profile HDF5 archive

Top-level command
-----------------

.. click:: quantas.cli.main:main
   :prog: quantas
