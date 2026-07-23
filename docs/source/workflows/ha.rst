Harmonic Approximation: implementation and workflow
===================================================

Purpose and scope
-----------------

The harmonic-approximation workflow evaluates vibrational thermodynamic
properties from phonon frequencies supplied at one or more fixed volumes.  It
implements the independent harmonic-oscillator expressions described in
:doc:`../theory/ha`, but it does **not** minimize the free energy with respect
to volume.  Every sampled volume remains an independent thermodynamic state.

This distinction is operationally important:

- HA answers *what are the harmonic properties at each supplied volume?*;
- QHA additionally asks *which volume minimizes* :math:`F(V,T)+PV` *at each
  pressure and temperature?*

The HA workflow is therefore both a complete calculation in its own right and
the first numerical stage of a QHA calculation.

Computational pipeline
----------------------

The workflow follows the sequence

.. code-block:: text

   phonon YAML
       │
       ├─ validate volumes, static energies, frequencies, and q-point weights
       ├─ normalize q-point weights
       ├─ convert temperature to K and frequencies to Hz
       ├─ evaluate harmonic oscillator sums on the complete T × V grid
       ├─ convert properties to the requested output units
       └─ build ResultData → HDF5 / report / plots / table export

The numerical thermodynamic backend is vectorized with NumPy.  Each property
is evaluated over all temperatures, q-points, modes, and volumes without a
Python loop over individual oscillators.

Scientific input contract
-------------------------

A normalized HA input contains at least:

``volume``
   One value for every sampled structure, with shape ``(nvol,)``.

``energy``
   Static electronic energy at every sampled volume, also with shape
   ``(nvol,)``.

``frequencies``
   Mode-resolved phonon frequencies with shape
   ``(nq, 3 * natoms, nvol)``.

``weights``
   One q-point weight for every q-point, with shape ``(nq,)``.

``natoms`` and ``formula_units``
   Cell normalization information used by reports and molar conversions.

Quantas checks that the number of modes is exactly :math:`3N`, that every
volume has one static energy and one complete phonon spectrum, and that the
q-point and weight axes agree.

Multi-volume HA
~~~~~~~~~~~~~~~

HA is natively multi-volume.  For every sampled volume :math:`V_i`, Quantas
calculates

.. math::

   U_0(V_i),\quad U_{\mathrm{zp}}(V_i),

and the temperature-dependent arrays

.. math::

   U_{\mathrm{th}}(T,V_i),\quad U(T,V_i),\quad
   S(T,V_i),\quad F_{\mathrm{vib}}(T,V_i),\quad
   F(T,V_i),\quad C_V(T,V_i).

The result therefore preserves the full sampled ``temperature × volume``
surface.  No volume is labelled as an equilibrium volume by HA.

q-point weights
~~~~~~~~~~~~~~~

The stored q-point weights are normalized internally by their sum before any
oscillator sum is evaluated.  Multiplying all weights by the same positive
constant therefore does not change the result.

Quantas does not currently perform an independent symmetry analysis to infer
q-point equivalence or multiplicity.  Consequently:

- a full q-point set may use equal weights when scientifically appropriate;
- an irreducible set must carry the correct multiplicities from its source;
- weights imported from a code such as Phonopy are consumed as supplied;
- users should not invent multiplicities solely from visual similarity of
  q-point coordinates.

The format and provenance requirements are described in
:doc:`../formats/phonon_yaml`.

Treatment of zero and imaginary frequencies
-------------------------------------------

The harmonic numerical functions exclude frequencies less than or equal to
zero from the oscillator sums.  This prevents undefined logarithms and Bose
factors from contaminating all properties, but it must not be mistaken for a
physical correction of an unstable structure.

A non-positive frequency may represent:

- a translational acoustic mode affected by numerical noise;
- an insufficiently converged force-constant calculation;
- an inconsistent supercell or sampling setup;
- a genuine dynamical instability.

Quantas does not decide which interpretation is correct.  Before accepting a
thermodynamic result, inspect how many modes were omitted, their magnitude,
their q-point location, and whether the behavior persists across volumes.
Removing a genuinely unstable mode can make a table finite without making the
harmonic model scientifically valid.

Temperature grid and the zero-temperature limit
------------------------------------------------

The public temperature range is defined by ``MIN MAX STEP``.  Quantas uses the
same inclusive grid convention in HA and QHA:

.. code-block:: python

   np.arange(MIN, MAX + STEP, STEP)

Choose ``MAX`` consistently with ``STEP``.  When the interval is not an exact
multiple of the step, this convention may include the next grid point beyond
the nominal maximum.

At :math:`T=0` the implementation uses explicit limits:

- :math:`U_{\mathrm{th}}=0`;
- :math:`S=0`;
- :math:`C_V=0`;
- :math:`F_{\mathrm{vib}}=U_{\mathrm{zp}}`.

The zero-point energy is evaluated once as a volume-dependent quantity and is
stored with shape ``(1, nvol)``.  Reports and plots broadcast it in temperature
only for presentation; the native result is not duplicated.

Units and normalization
-----------------------

Before the oscillator sums are evaluated:

- temperatures are converted to kelvin;
- frequencies are converted to hertz;
- q-point weights are normalized.

The low-level harmonic functions return energy-like quantities in
``kJ mol^-1`` and entropy/heat capacity in ``J mol^-1 K^-1``.  The workflow
then converts them to the requested Quantas energy convention while preserving
``float64`` values.  Display precision is applied only by renderers.

Default behavior
----------------

.. list-table:: HA defaults and rationale
   :header-rows: 1
   :widths: 24 20 56

   * - Control
     - Default
     - Rationale
   * - Temperature range
     - ``298.15 298.15 1.0`` K
     - Produces one ambient-temperature state unless the user requests a
       thermodynamic curve.
   * - Energy unit
     - ``Ha``
     - Preserves the native electronic-structure energy scale used by the
       normalized input.
   * - Volume length unit
     - ``A``
     - Cell volumes are expressed as the cube of this length unit.
   * - Frequency unit
     - ``cm^-1``
     - Conventional output of vibrational calculations and spectroscopy.
   * - Thermodynamic calculation
     - enabled
     - The public CLI always evaluates the requested HA properties.

The defaults define a safe minimal calculation, not a recommended production
temperature range.  The scientific question determines the appropriate lower
and upper temperatures and resolution.

Performance and memory use
--------------------------

The dominant work and temporary-array size scale approximately as

.. math::

   N_T\,N_q\,N_{\mathrm{mode}}\,N_V.

The calculation evaluates several thermodynamic functions sequentially, so a
large ``T × q × mode × volume`` product affects both execution time and peak
memory.

Practical acceleration strategies are:

1. use a coarse temperature grid during input validation;
2. calculate only the temperature interval relevant to the study;
3. avoid creating many nearly redundant temperature points for a smooth plot;
4. verify the final grid only after the phonon dataset is scientifically sound;
5. use ``--benchmark`` to report backend timings when performance is being
   investigated.

Reducing the number of sampled volumes or q-points is not a neutral numerical
optimization: it changes the physical input and must be justified at the
phonon-calculation level.

.. note::

   The current HA workflow has no numerical parameter with a default value of
   ``512``.  Its resolution is controlled by the supplied phonon dataset and
   the requested temperature grid.  A value named ``512`` encountered in a
   different Quantas module or in historical material must not be transferred
   to HA.

Reports, plots, and persistence
-------------------------------

The calculation returns a generic :class:`quantas.api.common.ResultData`
envelope containing a typed HA payload.  From the same result, the public API
can build:

- neutral report tables;
- neutral plot specifications;
- native HDF5 output;
- selected property tables.

The CLI follows the same sequence and automatically writes a deterministic log
for a significant calculation.  Plotting occurs after the numerical run and
does not alter or recompute the stored thermodynamic arrays.

Warnings and failure modes
--------------------------

The workflow stops before calculation when:

- required arrays are missing;
- volumes and static energies have inconsistent lengths;
- the phonon array does not have shape ``(nq, 3 * natoms, nvol)``;
- weights do not match the q-point axis;
- the sum of weights is not positive;
- the temperature range or step is invalid;
- a requested unit conversion is unsupported.

Non-positive modes are not a hard numerical failure because they are excluded
from the sums.  They remain a scientific validation responsibility, as
explained above.

Decision guide
--------------

.. list-table:: Practical HA decisions
   :header-rows: 1
   :widths: 34 32 34

   * - Situation
     - Recommended action
     - Reason
   * - First check of a new input
     - Run one or a few temperatures and inspect all volumes.
     - Detect shape, unit, and unstable-mode problems quickly.
   * - Low-temperature study
     - Include :math:`T=0` and use a finer grid where heat capacity changes
       rapidly.
     - The thermodynamic functions are most curved at low temperature.
   * - High-temperature trend only
     - Use a wider but coarser preliminary grid, then refine only if needed.
     - Harmonic curves are usually smoother at high temperature.
   * - Irreducible q-point input
     - Verify source-provided multiplicities before running.
     - Quantas normalizes weights but does not infer them.
   * - Several non-positive optical modes
     - Stop and re-evaluate the phonon calculation.
     - Excluding them may conceal a real instability.
   * - Need equilibrium volume or :math:`C_P`
     - Continue with QHA.
     - HA does not minimize volume and directly provides only :math:`C_V`.

Related documentation
---------------------

- Scientific theory: :doc:`../theory/ha`
- Complete worked example: :doc:`../tutorials/ha`
- Input specification: :doc:`../formats/phonon_yaml`
- CLI syntax: :doc:`../cli/ha`
- Public Python surface: :doc:`../api/ha`
