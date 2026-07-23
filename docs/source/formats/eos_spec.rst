EOS batch specification
=======================

An EOS specification is a strict, hand-written text file that describes an
ordered :class:`~quantas.api.eos.BatchPlan`.  It is separate from the
experimental data file: the data file contains measured values, units, and
experimental provenance, while the specification contains analysis choices.

The file extension has no meaning.  ``analysis.spec``, ``analysis.txt``,
``analysis.dat``, and a file without a suffix are treated identically.  Quantas
recognizes the format from the required first non-empty line::

   # QUANTAS EOS SPEC 1

Version 1 is deliberately small and strict.  Unknown sections, unknown keys,
duplicate sections, duplicate keys, invalid values, incompatible solver
options, unavailable targets, and invalid model parameters are errors.  No
scientific option is ignored silently.

Running a specification
-----------------------

Pass the data file as the positional argument and the independent specification
with ``--spec``:

.. code-block:: console

   quantas eos run PV_topaz.dat --spec topaz.spec

The specification becomes the authority for scientific settings.  Therefore
it cannot be mixed with scientific command-line options such as ``--domain``,
``--fit``, ``--solver``, ``--pv-eos``, ``--fix``, or input-unit overrides.
Operational and presentation options remain available:

.. code-block:: console

   quantas eos run PV_topaz.dat --spec topaz.spec \
      --output topaz.hdf5 --report topaz.log --quiet

``--verbosity extended``, ``--show-uncertainties``, and ``--max-data-rows`` may increase
or override presentation detail without changing the scientific plan.

Generating the commented template
---------------------------------

Create the complete version-1 template with:

.. code-block:: console

   quantas eos spec-template

The default destination is ``eos.spec``.  A path with any suffix, or no suffix,
may be supplied as the optional positional argument:

.. code-block:: console

   quantas eos spec-template topaz-analysis.txt

The command never overwrites an existing file silently.  Use ``--force`` only
when replacement is intentional:

.. code-block:: console

   quantas eos spec-template topaz-analysis.txt --force

A synchronized copy is also available as :download:`eos.spec <../_downloads/eos.spec>`.

The generated file is deliberately verbose and self-documenting.  It contains:

* the mandatory signature and syntax rules;
* every recognized section and ordinary key;
* pressure, length, and temperature unit overrides;
* all common solver and covariance controls;
* P--V, V--T, and P--V--T model lists;
* fixed, initial, and bounded parameter examples;
* acceptance and replacement semantics;
* short and extended presentation settings;
* commented examples for homogeneous, heterogeneous, thermal, and coupled jobs.

Only one minimal P--V volume job is active.  All larger examples are commented,
so the template is a valid specification immediately after generation and can
be checked safely with ``--dry-run``.  The distributed copy under
``examples/eos/eos.spec`` is generated from the same public template function
used by the command.

Dry-run validation
------------------

Use ``--dry-run`` before a real batch:

.. code-block:: console

   quantas eos run PV_topaz.dat --spec topaz.spec --dry-run

A dry run:

#. validates the signature and complete syntax;
#. reads and normalizes the data using ``[input]`` declarations;
#. resolves general and domain-specific defaults;
#. expands ``targets = all`` against the quantities actually present;
#. constructs typed EOS models, constraints, and solver options;
#. validates every request against the selected data and uncertainties;
#. prints the resolved plan and normalized requests;
#. performs no fit and creates no HDF5 archive.

A plain-text report may still be requested during a dry run.  This is useful for
recording the exact resolved plan before launching a calculation.

Syntax rules
------------

The syntax is INI-like, but it is parsed by Quantas rather than by a generic INI
configuration service.

* The first non-empty line must be exactly ``# QUANTAS EOS SPEC 1``.
* Section names use square brackets.
* Entries use ``key = value``.
* Blank lines are ignored.
* A line whose first non-whitespace character is ``#`` is a comment.
* Section and ordinary-key names are case-insensitive.
* Job names are retained as stable identifiers after normalization.
* Parameter names follow the public EOS fitting contracts.
* Duplicate keys and duplicate sections are errors.
* Inline comments are not supported in version 1; place comments on their own
  lines.

Minimal example
---------------

A BM2 volume fit using effective variance requires only:

.. code-block:: ini

   # QUANTAS EOS SPEC 1

   [defaults]
   solver = effective-variance

   [job volume]
   domain = pv
   targets = volume
   model = BM2

Unspecified values use documented Quantas defaults.  The default solver is OLS,
the default P--V model is BM3, the default V--T model is
``berman:quadratic``, and the default P--V--T coupling is ``linear``.

Sections
--------

``[metadata]``
~~~~~~~~~~~~~~

This optional section contains free-form textual provenance.  Values are stored
in the resolved batch manifest and native HDF5 archive.

.. code-block:: ini

   [metadata]
   title = Topaz pressure-volume analysis
   description = Independent volume and axial fits
   author = Research group name

Metadata does not modify the numerical calculation.

``[input]``
~~~~~~~~~~~

This optional section overrides input units before the data file is normalized.
Recognized keys are:

``pressure_unit``
   Pressure unit accepted by the Quantas unit converter.  The normalized unit
   is GPa.

``length_unit``
   Cell-length unit.  The corresponding cubic unit is applied to absolute
   volumes.  The normalized units are Angstrom and Angstrom cubed.

``temperature_unit``
   ``K``, ``C``, or ``F``.  Temperatures and their uncertainties are normalized
   consistently to kelvin.

.. code-block:: ini

   [input]
   pressure_unit = kbar
   length_unit = bohr
   temperature_unit = C

The original values and original unit labels remain in the HDF5 input record.
Relative ``V/V0`` and ``L/L0`` quantities remain dimensionless.

``[batch]``
~~~~~~~~~~~

The batch section currently accepts one key:

``failure_policy``
   ``stop`` or ``continue``.  A failed numerical attempt is persisted before
   the policy is applied.

.. code-block:: ini

   [batch]
   failure_policy = continue

``[defaults]``
~~~~~~~~~~~~~~

General defaults apply to every job unless a domain-specific default or job
entry overrides the same key.

Common keys are:

.. list-table:: General fit keys
   :header-rows: 1
   :widths: 28 72

   * - Key
     - Meaning
   * - ``solver``
     - ``ols``, ``wls``, ``effective-variance``, or ``odr``.
   * - ``covariance_scaling``
     - ``absolute``, ``reduced-chi-square``, or ``inflate-only``.
   * - ``max_iterations``
     - Positive solver iteration or evaluation limit.
   * - ``ftol``, ``xtol``, ``gtol``
     - Positive convergence tolerances where supported.
   * - ``inner_max_iterations``
     - Inner WLS evaluation limit for effective variance only.
   * - ``odr_difference``
     - ``central`` or ``forward`` for ODR only.
   * - ``odr_ndigit``
     - Positive number of reliable model digits for ODR only.
   * - ``accept``
     - ``yes`` or ``no``; successful jobs are accepted by default.
   * - ``replace_accepted``
     - ``yes`` or ``no``; replacement is never implicit.
   * - ``selection_base``
     - ``default`` respects input ``USE``/``*`` exclusions; ``all`` starts from every row.
   * - ``groups`` or ``include_groups``
     - Comma-separated positive group identifiers, or ``all`` for ``groups``.
   * - ``exclude_groups``
     - Comma-separated positive group identifiers removed from the selection.
   * - ``include_rows``
     - Comma-separated one-based data-row identifiers retained after the base/group selection.
   * - ``exclude_rows``
     - Comma-separated one-based data-row identifiers removed last.

Example:

.. code-block:: ini

   [defaults]
   solver = effective-variance
   covariance_scaling = inflate-only
   max_iterations = 200
   accept = yes
   replace_accepted = no

``[defaults.pv]``, ``[defaults.vt]``, and ``[defaults.pvt]``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Domain defaults are applied after ``[defaults]`` and before the individual job.
They may contain the common fit keys and parameter declarations.

P--V uses ``model``:

.. code-block:: ini

   [defaults.pv]
   model = BM3
   fix.KP = 4.0

V--T uses ``model``:

.. code-block:: ini

   [defaults.vt]
   model = berman:quadratic

P--V--T uses an explicit composition:

.. code-block:: ini

   [defaults.pvt]
   pv_model = BM3
   vt_model = berman:quadratic
   coupling = linear

For ``coupling = thermal-pressure``, ``vt_model`` must be omitted because the
coupling supplies its own thermal-pressure formulation.

The subordinate ``thermal_pressure_model`` key selects
``holland-powell-einstein`` (the default), ``mgd``, or
``mgd:q-compromise``.  ``mgd`` is the full model and exposes refinable
``theta_d0``, ``gamma0``, and ``q``.  ``mgd:q-compromise`` is the distinct
EosFit-style compromise approximation: it omits ``q`` rather than assigning a
hidden numerical value, and therefore should be selected explicitly when
reproducing a q-compromise reference fit.  A cell-volume MGD fit requires
either ``atoms_per_cell`` or the pair ``formula`` and
``formula_units_per_cell``.

.. code-block:: ini

   [job naf-mgd]
   domain = pvt
   targets = volume
   pv_model = BM4
   coupling = thermal-pressure
   thermal_pressure_model = mgd
   volume_basis = cell
   formula = NaF
   formula_units_per_cell = 4
   solver = effective-variance
   fix.temperature_ref = 295
   initial.theta_d0 = 459
   bound.theta_d0 = 1 : 5000
   initial.gamma0 = 1.5
   bound.gamma0 = 0.01 : 10
   initial.q = 1.0
   bound.q = -10 : 10

This fit uses only pressure, volume, and temperature observations. Adiabatic
bulk-modulus data and multi-observable objectives are outside the current MGD
workflow.


Data selection and groups
-------------------------

Selection settings are valid in ``[defaults]``, domain defaults, and individual
jobs. Resolution is deterministic:

#. start from ``selection_base = default`` (the input ``USE`` and ``*`` mask)
   or ``selection_base = all``;
#. retain ``groups``/``include_groups`` when declared;
#. remove ``exclude_groups``;
#. retain ``include_rows`` when declared;
#. remove ``exclude_rows``.

Rows are numbered from one in data-table order. ``groups`` and
``include_groups`` are aliases and cannot both appear in the same resolved job.
Unknown groups, out-of-range rows, and an empty final selection are errors.

.. code-block:: ini

   [job campaign-2]
   domain = pv
   targets = volume
   model = BM3
   groups = 2
   exclude_rows = 17

By default, an input point marked ``USE=0`` or with a trailing ``*`` cannot be
re-enabled by a group selection because the base is ``default``. Use
``selection_base = all`` explicitly when a job must reconsider every archived
observation:

.. code-block:: ini

   [job sensitivity-check]
   domain = pv
   targets = volume
   model = BM3
   selection_base = all
   include_groups = 1, 2
   exclude_rows = 4, 19
   accept = no

The resolved boolean mask, selection counts, and per-group counts are stored in
the request and HDF5 fit record. No data row is removed from the input dataset.

``[presentation]``
~~~~~~~~~~~~~~~~~~

Presentation settings do not change fitted values or HDF5 scientific content.
Detailed solver tracing is intentionally not a persistent specification setting:
request it operationally with ``quantas eos run ... --verbosity debug`` so the same
scientific specification can be used for normal and diagnostic executions.

``detail``
   ``short`` or ``extended``.

``show_uncertainties``
   ``yes`` or ``no``.

``max_data_rows``
   A positive integer, ``all``, or ``none``.  This limits display only.

.. code-block:: ini

   [presentation]
   detail = extended
   show_uncertainties = yes
   max_data_rows = 30

``[job NAME]``
~~~~~~~~~~~~~~

Every specification must contain at least one named job section.  Jobs are
executed in file order.  ``domain`` and ``targets`` are mandatory.

.. code-block:: ini

   [job volume]
   domain = pv
   targets = volume
   model = BM2

A job section may override inherited model, solver, constraint, acceptance, and
presentation-independent numerical settings.  ``note`` is persisted with the
fit record.

When one section selects several targets, Quantas expands it into one job per
target while retaining a stable identifier:

.. code-block:: ini

   [job axes]
   domain = pv
   targets = a, b, c
   model = PT3
   solver = ols

The resolved identifiers are ``axes-a``, ``axes-b``, and ``axes-c``.

Default precedence
------------------

For each job, settings are applied in the following order:

#. internal Quantas defaults;
#. ``[defaults]``;
#. the matching ``[defaults.pv]``, ``[defaults.vt]``, or ``[defaults.pvt]``;
#. ``[job NAME]``.

A later value replaces an earlier value with the same key.  Version 1 has no
inheritance between job sections and no conditional logic.

Targets
-------

Recognized targets are ``volume``, ``a``, ``b``, and ``c``.  Separate targets
with commas:

.. code-block:: ini

   targets = volume, a, b

``targets = all`` expands to every supported quantity present in the data file.
It cannot be mixed with explicit targets.  P--V--T currently supports only
``targets = volume``.

P--V models
-----------

Compact tags are recommended:

``M``
   Murnaghan.

``BM2``, ``BM3``, ``BM4``
   Birch--Murnaghan.

``PT2``, ``PT3``, ``PT4`` or ``NS2``, ``NS3``, ``NS4``
   Natural strain/Poirier--Tarantola.

``V2``, ``V3``
   Vinet.

``T2``, ``T3``, ``T4``
   Tait.

The documented full family names and their normalized aliases are also valid.
The compact tags above are preferred because they make the selected order
immediately visible in a batch plan.

.. code-block:: ini

   [job volume]
   domain = pv
   targets = volume
   model = birch-murnaghan3

V--T models
-----------

Supported tags are:

* ``berman:linear`` and ``berman:quadratic``;
* ``fei:linear`` and ``fei:inverse-square``;
* ``mhp:simplified`` and ``mhp:general``;
* ``salje``;
* ``khp``.

.. code-block:: ini

   [job thermal-volume]
   domain = vt
   targets = volume
   model = berman:quadratic
   solver = wls

P--V--T models
--------------

A P--V--T job declares each component explicitly:

.. code-block:: ini

   [job pvt-volume]
   domain = pvt
   targets = volume
   pv_model = BM3
   vt_model = berman:quadratic
   coupling = anderson-gruneisen
   solver = effective-variance

Available couplings are ``linear``, ``anderson-gruneisen``, and
``thermal-pressure``.

Parameter constraints
---------------------

Constraint keys use the same names as the public fitting functions.

``fix.NAME``
   Hold a parameter at the supplied finite value.

``initial.NAME``
   Override the initial value of a free parameter.

``bound.NAME``
   Set ``LOW : HIGH`` bounds for a free parameter.  ``none`` denotes an open
   side.  A bound declaration requires a matching ``initial.NAME``.

.. code-block:: ini

   initial.K0 = 160.0
   bound.K0 = 1.0 : 500.0
   fix.KP = 4.0

A fixed parameter cannot also declare an initial value or bounds.  Parameters
must exist in the selected model and target.  Parameters implied by an EOS
order cannot be made free or fixed; the dry run detects this when the complete
parameter map is prepared.

For volumetric P--V fits, common names are ``K0``, ``KP``, ``KPP``, and ``V0``.
For axial P--V fits, use ``M0``, ``MP``, ``MPP``, and ``L0``.  For V--T
fits, the reference quantity is always ``V0`` for volume and ``L0`` for a
cell axis.  Other thermal and P--V--T names follow the API, such as
``temperature_ref``, ``alpha0``, ``theta_e``, ``dK0_dT``, and ``delta``.
The historical ``X0`` and generic ``value_ref`` names are rejected so that
specification files have one stable public nomenclature.

Solver-specific validation
--------------------------

Method-specific keys are accepted only by the corresponding solver.  For
example, this is invalid:

.. code-block:: ini

   solver = ols
   odr_ndigit = 12

Likewise, WLS requires positive dependent-variable uncertainties, effective
variance requires at least one usable uncertainty component, and ODR requires
positive uncertainties for every fitted coordinate.  These requirements are
checked by ``--dry-run`` without invoking the ODR backend.

Acceptance and replacement
--------------------------

A successful job is accepted by default.  Set ``accept = no`` to persist the
attempt without making it the current result.

Two jobs may target the same scientific slot, but replacement must be explicit:

.. code-block:: ini

   [job volume-bm2]
   domain = pv
   targets = volume
   model = BM2
   accept = yes

   [job volume-bm3]
   domain = pv
   targets = volume
   model = BM3
   accept = yes
   replace_accepted = yes

Both immutable records remain in the archive; the BM3 record becomes the
current accepted result.

Complete example
----------------

.. code-block:: ini

   # QUANTAS EOS SPEC 1

   [metadata]
   title = Topaz EOS analysis
   description = Volume and axial fits with independent models

   [input]
   pressure_unit = GPa
   length_unit = angstrom
   temperature_unit = K

   [batch]
   failure_policy = continue

   [defaults]
   solver = effective-variance
   covariance_scaling = inflate-only
   max_iterations = 200
   accept = yes
   replace_accepted = no

   [defaults.pv]
   model = BM3

   [presentation]
   detail = extended
   show_uncertainties = yes
   max_data_rows = 30

   [job volume]
   domain = pv
   targets = volume
   model = BM2

   [job axis-a]
   domain = pv
   targets = a
   model = PT3
   solver = odr
   odr_difference = central
   odr_ndigit = 12
   initial.M0 = 450.0
   initial.MP = 10.0
   bound.M0 = 1.0 : 2000.0
   bound.MP = 0.0 : 30.0

   [job axes-bc]
   domain = pv
   targets = b, c
   model = BM3
   solver = ols

Error diagnostics
-----------------

Diagnostics identify the source path, source line, and section whenever that
information is available.  Examples include:

* ``unknown section [solver]``;
* ``unknown key 'solvert'``;
* ``unknown solver 'odl'``;
* ``odr_ndigit is valid only for solver = odr``;
* ``parameter 'banana' is not available for pv/volume``;
* ``dataset does not contain target(s): c``;
* ``thermal-pressure coupling defines its own thermal model; remove vt_model``.

Version 1 intentionally reports unknown names directly.  Automatic
``did-you-mean`` suggestions are not part of the initial parser.

Python API
----------

Parsing and dataset-dependent resolution are separate:

.. code-block:: python

   from quantas.api import eos

   document = eos.read_spec("topaz.spec")
   dataset = eos.read_input(
       "PV_topaz.dat",
       **document.input_options.as_dict(),
   )
   resolved = eos.resolve_spec(document, dataset)

   plan = resolved.plan
   report_options = resolved.report_options

The parser does not execute fits and does not depend on Click or Rich.  The
resolved plan is the same :class:`~quantas.api.eos.BatchPlan` consumed
by :func:`quantas.api.eos.run_batch`.

Template generation API
-----------------------

The same commented template used by :command:`quantas eos spec-template` can
be written through :func:`quantas.api.eos.write_spec_template`.  Template
generation remains separate from parsing and resolution, so producing a
starter file never executes or modifies a fitting plan.
