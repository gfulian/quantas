Energy--volume EOS from static electronic-structure data
========================================================

This tutorial exercises the public ``ev/energy`` workflow with a seven-point
MgO static energy--volume dataset generated from CRYSTAL calculations.  The
input is already normalized to the ordinary Quantas EOS text format so the
example focuses on fitting and post-fit use.  In a production calculation the
same file can be created directly from backend outputs with
:command:`quantas eos inpgen`.

The example data are available as
:download:`EV_mgo_pbe.dat <../_downloads/EV_mgo_pbe.dat>`.

Inspecting the data
-------------------

The file declares volume, the six cell parameters, and the authoritative total
energy for each state::

   UNITS V=angstrom^3 A=angstrom B=angstrom C=angstrom \
         ALPHA=degree BETA=degree GAMMA=degree E=Ha
   FORMAT V A B C ALPHA BETA GAMMA E

The primary Energy EOS fit uses ``V`` and ``E``.  When the complete cell metrics
are present, Quantas also treats them as a volume-aligned structural path and
derives the crystallographic response from the same E--V fit.  During input
generation, pressure is not taken from an arbitrary stress tensor;
pressure is reconstructed from the fitted energy relation as

.. math::

   P(V) = -\frac{dE}{dV}.

Run a BM3 fit
-------------

Fit a third-order Birch--Murnaghan integrated EOS with ordinary least squares:

.. code-block:: console

   quantas eos run EV_mgo_pbe.dat --domain ev --ev-eos BM3 \
      --output mgo_ev.hdf5 --report mgo_ev.log --force

The public E--V parameter convention is:

- ``E0`` in Ha;
- ``V0`` in angstrom cubed;
- ``K0`` in GPa;
- ``KP`` dimensionless;
- ``KPP`` in GPa :math:`^{-1}`.

For the bundled MgO dataset the BM3 fit gives approximately:

.. list-table:: MgO BM3 reference fit
   :header-rows: 1
   :widths: 28 28 44

   * - Parameter
     - Value
     - Unit
   * - ``E0``
     - -275.173937178
     - Ha
   * - ``V0``
     - 18.817428245
     - angstrom cubed
   * - ``K0``
     - 178.761458
     - GPa
   * - ``KP``
     - 3.815504
     - 1
   * - ``KPP``
     - -0.02091296
     - GPa :math:`^{-1}`

The energy RMSE is about :math:`3.28\times10^{-6}` Ha.  These values are a
regression target for the public workflow, not an experimental benchmark for
MgO or a recommendation that BM3 is universally preferable to other models.

The standard E--V report also lists every sampled volume together with the
ab initio energy, the EOS-calculated energy, the energy residual, and the
pressure reconstructed from :math:`-dE/dV`.  This makes the quality of the
energy fit and the associated static pressure scale visible without requiring
a separate diagnostics export.

Model discovery and alternatives
--------------------------------

List the E--V models exposed by the active installation with:

.. code-block:: console

   quantas eos show-models --domain ev

The integrated Murnaghan, Birch--Murnaghan, natural-strain,
Vinet, modified-Tait, and stabilized-jellium families are available according
to their registered model orders.  SJEOS is an E--V model; its analytical
pressure derivative is available to the calculator even though direct P--V
fitting is not exposed for SJEOS.

A useful sensitivity check is therefore to repeat the fit with, for example,
``T3`` or ``SJ`` and compare residual structure and derived pressure rather than
selecting a model from energy RMSE alone.


Structural response
-------------------

The curated primitive-cell MgO dataset declares ``SYSTEM cubic``, space group
225, ``CRYSTAL_REFERENCE primitive``, and ``CELL_MULTIPLICITY 4``.  The shared
structural-path backend therefore identifies the exact cubic branch and derives

.. math::

   \eta_a=\frac{\partial\ln a}{\partial\ln V}=\frac13,
   \qquad M_a=\frac{K}{\eta_a}=3K.

At the BM3 equilibrium point the primitive-cell fit gives approximately
:math:`M_a=536.2844` GPa.  The companion
:download:`EV_mgo_pbe_crystallographic.dat <../_downloads/EV_mgo_pbe_crystallographic.dat>`
contains the same physical path normalized to the conventional FCC cell.  Its
E0 and V0 are four times the primitive-cell values, while K0 and KP are
unchanged; the equilibrium conventional cubic axis is about 4.2222125 angstrom.

The primary structural response is available for any integrated Energy EOS
with the required derivatives, including ``SJ``.  It is not obtained by fitting
energy directly against :math:`a^3`; doing so would only be thermodynamically
equivalent to E(V) for special geometries such as a cubic conventional cell.

Optional pressure-form axial EOS parameterization
-------------------------------------------------

For direct comparison with pressure-based axial EOS work, request a secondary
pressure-form model independently of the Energy EOS:

.. code-block:: console

   quantas eos run EV_mgo_pbe_crystallographic.dat --domain ev --ev-eos SJ \
      --axial-eos BM3 --output mgo_ev_axial.hdf5 --report mgo_ev_axial.log --force

Quantas first fits ``SJ`` to E(V), derives the pressure and its propagated
covariance, and only then fits BM3 to :math:`P(a^3)`.  The report labels this as
a secondary axial fit.  The full covariance among the derived pressures is
stored in HDF5, while the current WLS solver uses the marginal pressure
uncertainties; that diagonal approximation is recorded explicitly.

Diagnostics
-----------

Inspect the accepted E--V record with:

.. code-block:: console

   quantas eos diagnose mgo_ev.hdf5 --slot ev/energy \
      --output mgo_ev_diagnostics.csv

The diagnostic table contains the observed and calculated energy, the energy
residual, and the EOS-derived pressure at each sampled volume.  If an external
dataset also supplies a scientifically meaningful pressure column, Quantas
retains it separately as ``source_pressure`` and reports
``source_pressure - eos_pressure``.  The CRYSTAL input generator deliberately
does not manufacture such pressures for constrained-volume optimizations.

Post-fit calculation
--------------------

Evaluate the EOS on a pressure grid:

.. code-block:: console

   quantas eos calculate mgo_ev.hdf5 --slot ev/energy \
      --pressure-range 0:20:2 --output mgo_ev_properties.csv

For an E--V record the calculator solves the inverse relation when pressure is
provided and reports volume, energy, bulk modulus, ``K'``, and ``K''``.  A
volume can instead be supplied directly with ``--coordinate`` or
``--coordinate-range``.

At the fitted equilibrium volume, the reconstructed pressure is zero within
numerical fit tolerance.  Tests should therefore compare zero pressure with the
*fitted* ``V0`` rather than with an antecedent theoretical value used to
construct a synthetic dataset.

Plots
-----

Inspect the available representations:

.. code-block:: console

   quantas eos plot mgo_ev.hdf5 --slot ev/energy --list-plots

The ordinary E--V inventory contains the fitted energy curve, the derived
``P(V)`` curve, and energy residuals.  Render all available plots with:

.. code-block:: console

   quantas eos plot mgo_ev.hdf5 --slot ev/energy

Python API
----------

The same fit is available through :mod:`quantas.api.eos` without invoking the
CLI:

:download:`Download ev_fit_api.py <../_downloads/ev_fit_api.py>`

.. literalinclude:: ../_downloads/ev_fit_api.py
   :language: python
   :linenos:

A second API example shows the shared structural response and the optional
secondary axial fit:

:download:`Download ev_structural_response_api.py <../_downloads/ev_structural_response_api.py>`

.. literalinclude:: ../_downloads/ev_structural_response_api.py
   :language: python
   :linenos:

Generating the dataset from CRYSTAL
-----------------------------------

For separate CRYSTAL outputs, create a text file containing one output path per
line and run:

.. code-block:: console

   quantas eos inpgen crystal-files.txt --interface crystal --list \
      --jobname "MgO static energy-volume EOS" -o mgo_ev.dat

The default keeps the primitive-cell normalization.  To write the conventional
crystallographic cell instead, use:

.. code-block:: console

   quantas eos inpgen crystal-files.txt --interface crystal --list \
      --crystal-reference crystallographic -o mgo_ev_crystallographic.dat

Energy and volume are scaled together by the cell multiplicity, and the cell
parameters are transformed with the same fixed reference transformation.

One CRYSTAL source may itself be a native multi-volume ``EOS`` calculation.  A
list may therefore combine single-state and multi-state sources.  Quantas
requires compatible atom count, chemical composition, total-energy correction
semantics, and non-duplicate volumes before writing one normalized Energy EOS
dataset.

Generating the dataset from VASP
--------------------------------

For VASP, put one completed calculation directory per line in the list file::

   eos/01
   eos/02
   eos/03
   eos/00
   eos/04
   eos/05
   eos/06

Then generate the same Quantas Energy EOS format with:

.. code-block:: console

   quantas eos inpgen vasp-runs.txt --interface vasp --list \
      --jobname "MgO static energy-volume EOS" -o mgo_vasp_ev.dat

Each directory must contain ``vasprun.xml`` and may also contain ``OUTCAR``.
The current adapter accepts exactly one ionic state from each source: use a
dedicated static E--V series rather than inserting an optimization history.
For the zero-electronic-temperature E--V workflow Quantas selects VASP
``energy(sigma->0)`` and rejects a list whose smearing, Brillouin-zone sampling,
pseudopotentials, or other relevant electronic settings define incompatible
energy surfaces.  This means that an optimized reference geometry may be reused
for a static EOS calculation, but its optimization-run energy is not mixed into
the static EOS dataset.

VASP cells are reduced to their primitive normalization where translational
symmetry permits it.  If the run already uses a primitive cell, its original
basis is preserved.  ``--crystal-reference crystallographic`` can then be used
exactly as for CRYSTAL to express the generated table in one fixed
crystallographic-cell normalization.
