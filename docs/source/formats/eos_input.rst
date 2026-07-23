EOS text input
==============

EOS reads keyword-directed tabular files. The ``FORMAT`` declaration is the
single authority that defines the meaning and order of the numeric columns.
The table continues until end of file, so a row count is not required. This
keeps experimental data independent from the fit plan: EOS family, order,
parameter constraints, and numerical solver belong to an
:class:`~quantas.api.eos.FitRequest` or a batch specification. The
measurement table may carry only non-destructive data organization through
``GROUP`` and ``USE`` columns or a trailing ``*`` exclusion marker.

The reader accepts both historical Quantas files and common EosFit-style files.
Keyword matching is case-insensitive, commas and whitespace can both separate
columns, and a trailing colon is accepted. Thus all of the following are
valid::

   FORMAT P V SIGP SIGV
   Format 1, P, V, SigP, SigV
   format: T, sigT, P, sigP, V, sigV

Minimal examples
----------------

EosFit-style P-V table:

.. code-block:: text

   TITLE Quartz compression data
   FORMAT 1 P V SIGP SIGV
   0.001 112.981 0.000001 0.002
   0.430 111.725 0.008000 0.014

Historical Quantas layout:

.. code-block:: text

   JOB
   Quartz compression data
   FORMAT
   P V sigmaP sigmaV
   DATA
   0.001 112.981 0.000001 0.002

Simple isobaric V-T table without a pressure column:

.. code-block:: text

   TITLE Triclinic phase thermal expansion
   FORMAT 1 T V
   91.8534 667.631
   95.1401 667.641

A constant pressure may also be included as an experimental condition:

.. code-block:: text

   TITLE Rutile relative expansion
   TSCALE K
   VSCALE V/V0
   LSCALE L/L0
   FORMAT: T, SIGT, P, SIGP, V, SIGV, A, SIGA, C, SIGC
   303 3 0.0001 0.000001 1.0000 0.0005 1.0000 0.0003 1.0001 0.0003

Here the pressure column records an isobaric reference condition. Because it is
constant, it cannot constrain a P-V EOS; it remains valid and useful metadata
for a V-T analysis.

Comments
--------

Text following ``#`` is excluded from numeric parsing. Full-line comments can
also use the case-insensitive ``COMMENT`` keyword:

.. code-block:: text

   COMMENT Data selected from the literature and rescaled to V0
   # Individual references can also be grouped in the table

   FORMAT T V
   # Smith et al. (2020)
   300 1.0000
   400 1.0031  # point checked manually

Full-line ``#`` comments, inline ``#`` text, and ``COMMENT`` values are
retained in ``dataset.metadata["comments"]`` while remaining excluded from
numeric parsing. This preserves source notes for future HDF5 persistence
without treating terminal prose as scientific workflow history.

Keywords
--------

``TITLE`` or ``JOB``
   Optional free-text dataset description. The value may follow on the same
   line or on the next non-comment line.

``COMMENT``
   Optional free-text note retained in dataset metadata. Repeated declarations
   are allowed.

``FORMAT``
   Mandatory ordered column declaration. It may appear on the same line or the
   next non-comment line. ``FORMAT 1 ...`` is accepted for compatibility with
   one-line EosFit records. Multi-line records are not yet supported.

``DATA``
   Optional marker for the start of the numeric table. When omitted, the first
   non-keyword line after ``FORMAT`` starts the table.

``SYSTEM``
   Optional crystal-system label normalized to one of ``triclinic``,
   ``monoclinic``, ``orthorhombic``, ``tetragonal``, ``trigonal``,
   ``hexagonal``, or ``cubic``. ``rhombohedral`` is accepted as an alias of
   ``trigonal``. The current fitting workflows retain the system as scientific
   metadata and record the independent conventional axes; they do not invent
   missing cell parameters or impose directional constraints.

``TSCALE``
   Optional temperature scale: ``K``, ``C``, or ``F``. Celsius and Fahrenheit
   values are converted to Kelvin. Original values and their input scale remain
   available in the raw columns and metadata.

``VSCALE``
   Declares how volume-like values are represented. ``VSCALE V/V0`` marks
   ``V`` and ``SIGV`` as dimensionless ratios relative to a reference volume.
   ``VSCALE absolute`` restores the default dimensional interpretation.

``LSCALE``
   Declares how linear quantities are represented. ``LSCALE L/L0`` marks
   ``A``, ``B``, ``C`` and their sigma columns as dimensionless ratios relative
   to reference lengths. ``LSCALE absolute`` is the default.

``UNITS``
   Optional ``COLUMN=UNIT`` declarations, for example
   ``UNITS P=GPa V=angstrom^3``. No silent conversion is currently applied
   except for temperature scales. A normalized ``VSCALE`` or ``LSCALE`` is
   incompatible with a dimensional unit declaration and is rejected rather
   than silently reinterpreted.

``PROVENANCE``
   Optional free-text source description attached to every provided column.

Columns
-------

Column names are case-insensitive and accept short aliases.

==================== ================================
Canonical quantity   Common input aliases
==================== ================================
pressure             ``P``, ``PRESSURE``
sigma_pressure       ``SIGP``, ``SIGMAP``
temperature          ``T``, ``TEMP``, ``TEMPERATURE``
sigma_temperature    ``SIGT``, ``SIGMAT``
volume               ``V``, ``VOL``, ``VOLUME``
sigma_volume         ``SIGV``, ``SIGMAV``
energy               ``E``, ``ENERGY``
sigma_energy         ``SIGE``, ``SIGMAE``
a, b, c              ``A``, ``B``, ``C``
sigma_a/b/c          ``SIGA``, ``SIGB``, ``SIGC``
alpha, beta, gamma   full angle names
sigma angles         ``SIGALPHA``, ``SIGBETA``, ``SIGGAMMA``
group                ``GROUP``, ``GRP``
use                  ``USE``, ``INCLUDE``
==================== ================================

All normalized numerical arrays use ``float64``. Every uncertainty column must
have its corresponding measured quantity, and negative uncertainties are
rejected. Raw and normalized columns are both retained by ``EOSDataset`` for
future self-contained HDF5 persistence.


Non-destructive selection and groups
------------------------------------

A dataset may identify experimental groups and its default fitting selection
without deleting observations. ``GROUP`` values are positive integers. ``USE``
accepts ``1``/``0``, ``yes``/``no``, or equivalent boolean words. A trailing
``*`` excludes a row regardless of the ``USE`` value and is convenient when
editing an existing EosFit-style table:

.. code-block:: text

   TITLE Two experimental campaigns
   SYSTEM tetragonal
   FORMAT P V SIGP SIGV GROUP USE
   0.0 100.00 0.01 0.02 1 1
   1.0  99.30 0.01 0.02 1 1 *
   0.0 100.05 0.02 0.03 2 1
   2.0  98.61 0.02 0.03 2 0

The table contains four observations, but only the first and third are selected
by default. Quantas reports total, selected, and excluded counts overall and for
each group. Every row remains available through :class:`EOSDataset` and is
stored in the native HDF5 archive.

The Python API uses the dataset default selection whenever an
:class:`EOSFitRequest` has no explicit mask. A strict EOS batch specification
can start from either the default selection or all rows and can then select
groups or one-based row numbers. See :doc:`eos_spec`.

.. code-block:: python

   from quantas.api import eos

   dataset = eos.read_input("grouped.dat")
   print(dataset.selected_npoints, dataset.excluded_npoints)
   print(dataset.group_summary())

An explicit API mask is a final selection and may intentionally override the
input default. This permits a GUI or analysis script to re-enable a previously
excluded observation without modifying the source table.

Absolute and normalized structural data
---------------------------------------

Without scale keywords, the default units are ``angstrom^3`` for volume and
``angstrom`` for cell lengths. With ``VSCALE V/V0`` or ``LSCALE L/L0``, the
stored normalized values and uncertainties are dimensionless. No multiplication
by an unknown reference value is attempted.

Relative data remain sufficient for thermal-expansion coefficients because a
constant scale factor cancels in logarithmic derivatives:

.. math::

   \alpha_V = \frac{1}{V}\left(\frac{\partial V}{\partial T}\right)_P,
   \qquad
   \alpha_x = \frac{1}{x}\left(\frac{\partial x}{\partial T}\right)_P.

A V-T fit result reports whether its reference
quantity is an absolute volume/length or a dimensionless ratio. The reader
therefore stores a per-column scale map in
``dataset.metadata["column_scales"]``.

Coordinate classification
-------------------------

The presence of a column does not prove that it can act as an independent
fitting coordinate. EOS classifies each measured quantity using its selected
numeric range:

.. code-block:: python

   from quantas.api import eos

   dataset = eos.read_input("rutile.dat")
   classification = dataset.classify()

   classification.is_isobaric
   classification.reference_pressure
   classification.variable_coordinates
   classification.constant_coordinates

For each coordinate, classification records the minimum, maximum, span, unit,
constancy tolerance, and—only for a constant coordinate—the reference value. Classification can be recomputed for a mask,
which is important when a subset of a larger P-T table represents one isobar or
one isotherm.

A constant pressure is interpreted as an isobaric reference condition. It is
valid for V-T or L-T data but cannot constrain a pressure-volume equation. The
P-V API rejects such a request before solver dispatch with a scientific error
message. A constant temperature is handled analogously as an isothermal
reference condition.

Selecting V-T data
------------------

The reader and public dataset contract support frontend-neutral V--T series
independently of the thermal model selected for fitting:

.. code-block:: python

   from quantas.api import eos

   dataset = eos.read_input("rutile.dat")
   vt = dataset.series("volume", independent="temperature")
   at = dataset.series("a", independent="temperature")

The returned series contains available sigma columns, units, selection mask,
and absolute or normalized scale metadata.  Selecting a series does not choose
a thermal-expansion equation or numerical solver.

Data and fit plans
------------------

The tabular file describes measurements, units, scales, reference conditions,
and provenance, plus an optional default non-destructive selection. EOS
family, order, fixed parameters, bounds, job-specific masks, and statistical
method belong to an :class:`~quantas.api.eos.FitRequest` or a batch
specification. This separation allows one experimental dataset to be reused in
several fits without deleting or rewriting observations.

Compatibility data
------------------

The test suite contains unchanged reference files in all currently supported
layouts:

``PV_quartz.dat``
   EosFit-style ``FORMAT 1 P V SigP SigV`` P-V data.

``PV_topaz.dat``
   Historical Quantas ``JOB``/``FORMAT``/``DATA`` P-V data with volume and
   ``a``, ``b``, and ``c``.

``T_triclinic.dat``
   Simple ``FORMAT 1 T V`` thermal-expansion table without pressure.

``rutile.dat``
   Mixed-case keywords, ``COMMENT`` records, ``FORMAT:``, constant pressure,
   and normalized ``V/V0`` and ``L/L0`` quantities.


P--V--T tables
--------------

A global P--V--T fit uses the same declared-column syntax. The minimum table is

.. code-block:: text

   TITLE Example P-V-T dataset
   FORMAT P SIGP V SIGV T SIGT
   0.0  0.02  100.00  0.002  300.0  0.2
   5.0  0.02   97.10  0.002  300.0  0.2
   0.0  0.02  101.00  0.002  600.0  0.2
   5.0  0.02   98.00  0.002  600.0  0.2

Column order is arbitrary because ``FORMAT`` remains authoritative. Pressure,
volume, and temperature must all vary within the selected observations for a
global P--V--T fit. A constant-pressure table is classified as isobaric and is
valid for V--T; a constant-temperature table is classified as isothermal and
is valid for P--V. Storing the constant coordinate is useful provenance but
does not add fitting information.

``sigma_pressure`` is used by WLS. Effective variance can use all available
``sigma_pressure``, ``sigma_volume``, and ``sigma_temperature`` terms. P--V--T
ODR requires all three uncertainty columns to be finite and strictly positive
for every selected observation.

The data file does not select the P--V EOS, V--T EOS, coupling prescription,
FREE/FIXED parameters, or numerical solver. Those choices belong to
``EOSFitRequest`` and, later, to the batch plan, so the same measurements can be
reused in global and controlled analyses.


CLI unit overrides
------------------

:command:`quantas eos run` accepts ``--pressure-unit``, ``--length-unit``, and
``--temperature-unit``. Explicit command-line units override declarations in
the input file. Without either source, EOS assumes GPa, Angstrom/Angstrom
cubed, and kelvin. Absolute values and their standard uncertainties are
normalized together; the raw arrays and original unit labels remain available
in :class:`~quantas.api.eos.Dataset` and in the native HDF5 archive.

The length unit defines both cell-length columns and the cube used by absolute
volume columns. Supported metric length scales include metre, centimetre,
millimetre, micrometre, nanometre, and picometre, in addition to Angstrom and
Bohr. Relative declarations ``VSCALE V/V0`` and ``LSCALE L/L0`` remain
dimensionless and are not converted by a physical-unit override.
