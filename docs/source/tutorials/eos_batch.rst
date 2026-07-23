Controlling EOS fitting with a batch specification
==================================================

Why use a specification file?
-----------------------------

A direct :command:`quantas eos run` command is convenient for one fit.  A
scientific comparison usually requires several related fits that differ in
model, solver, constraints, selected rows, or fitted target.  Repeating long
commands by hand makes the comparison difficult to audit.

A ``QUANTAS EOS SPEC 1`` file records the complete fitting plan.  It is a plain
text, strict, declarative document.  It contains no Python code and does not
change the input data.  The same specification can be resolved by the CLI or by
:mod:`quantas.api.eos`.

Generate an annotated template with

.. code-block:: console

   quantas eos spec-template eos.spec

The three specifications used by this tutorial can be downloaded here:

- :download:`Quartz P--V specification <../_downloads/tutorials/eos/quartz_pv_tutorial.spec>`
- :download:`Rutile V--T specification <../_downloads/tutorials/eos/rutile_vt_tutorial.spec>`
- :download:`NaF P--V--T specification <../_downloads/tutorials/eos/naf_pvt_tutorial.spec>`

Resolution order
----------------

Settings are resolved in a predictable order:

.. code-block:: text

   Quantas defaults
       → [defaults]
       → [defaults.pv], [defaults.vt], or [defaults.pvt]
       → [job NAME]

A job-specific value therefore overrides a domain default, which overrides a
common default.  This makes it possible to define one standard solver policy
and change only the few jobs needed for a sensitivity test.

Dry-run validation
------------------

Always validate a non-trivial batch before running it:

.. code-block:: console

   quantas eos run examples/eos/PV_quartz.dat \
       --spec examples/eos/specs/quartz_pv_tutorial.spec \
       --dry-run

The dry run:

- reads and normalizes the dataset;
- resolves defaults and job overrides;
- expands ``targets = all`` or target lists;
- validates model/domain compatibility;
- checks required uncertainty columns;
- validates fixed values, starting values, and bounds;
- reports the selected and excluded observations;
- does **not** call a solver;
- does **not** create an HDF5 archive.

This is especially important for P--V--T fits, where model composition and
normalization choices are more complex.

Choosing a regression method
----------------------------

Quantas exposes four fitting strategies.  They share the same model evaluator
and result contract, but they represent different assumptions about the
observations.

Ordinary least squares: ``ols``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

OLS minimizes unweighted residuals in the dependent quantity.  For a P--V fit,
pressure is the fitted response and volume is the coordinate.  Input standard
uncertainties are ignored.

Use OLS when:

- uncertainties are unavailable;
- the data are already homogeneous and similarly precise;
- a deliberately unweighted reference fit is required.

Do not interpret OLS parameter errors as a complete experimental uncertainty
model when the observations have strongly unequal precision.

Weighted least squares: ``wls``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

WLS weights the response residual by its standard uncertainty.  It requires a
positive response uncertainty for every selected observation.  Coordinate
uncertainty is not projected into the response.

For P--V data this normally means using :math:`\sigma_P` while treating volume
as exact.  For V--T data it normally means using :math:`\sigma_V` while treating
temperature as exact.

Effective variance: ``effective-variance``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Effective variance iteratively projects coordinate uncertainties into the
response using the local derivative of the EOS.  In a P--V fit, both
:math:`\sigma_P` and :math:`\sigma_V` contribute to the effective pressure
variance.  In a V--T fit, both structural and temperature uncertainties can
contribute.  In P--V--T, both volume and temperature coordinate terms are
projected into pressure.

This method is often a practical EosFit-like choice when complete standard
uncertainties are available and the explicit EOS form is retained.

Orthogonal distance regression: ``odr``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

ODR minimizes a weighted distance in response and coordinate space.  It can
adjust the coordinates as well as the response and therefore represents a
different objective from effective variance.  Quantas uses the ODRPACK95
backend and requires complete positive standard uncertainties for every fitted
coordinate and response.

A representative job is:

.. code-block:: ini

   [job bm3-odr]
   domain = pv
   targets = volume
   model = BM3
   solver = odr
   odr_difference = central
   odr_ndigit = 12
   accept = no

ODR should not be selected simply because it is the most elaborate solver.
Compare it with WLS and effective variance and inspect whether the fitted
parameters change by more than their uncertainties.

Covariance scaling
------------------

Weighted fits in these tutorials use

.. code-block:: ini

   covariance_scaling = inflate-only

When the reduced chi-square is greater than one, parameter covariance is
inflated.  When it is less than one, the reported covariance is not reduced.
This EosFit-like policy avoids making parameter errors artificially smaller
because the supplied experimental uncertainties were conservative.

Candidate records and accepted results
--------------------------------------

Every attempted fit is stored as an immutable record.  A successful record may
also become the accepted result for its ``domain/target`` slot.

In the tutorial specifications:

- sensitivity jobs use ``accept = no``;
- one reference job uses ``accept = yes``;
- no previous result is silently replaced.

The archive therefore retains the complete comparison while post-fit commands
can use one unambiguous accepted slot such as ``pv/volume``.

Running the same batch from Python
----------------------------------

:download:`Download run_spec_api.py <../_downloads/tutorials/eos/run_spec_api.py>`

.. literalinclude:: ../_downloads/tutorials/eos/run_spec_api.py
   :language: python
   :linenos:

For example:

.. code-block:: console

   python examples/eos/scripts/run_spec_api.py \
       examples/eos/PV_quartz.dat \
       examples/eos/specs/quartz_pv_tutorial.spec \
       quartz_pv.hdf5 \
       --report quartz_pv.log \
       --force

The script uses only :mod:`quantas.api.eos` and the public rendering entry
points.  It validates every expanded request before executing the batch.

Exercise: add ODR without losing the reference fit
---------------------------------------------------

Copy ``quartz_pv_tutorial.spec`` and add the ``bm3-odr`` job shown above.
Keep ``accept = no`` so that the effective-variance BM3 record remains the
accepted result.  Run the new batch and complete the table:

.. list-table::
   :header-rows: 1

   * - Solver
     - :math:`V_0`
     - :math:`K_0`
     - :math:`K'_0`
     - RMSE
     - Reduced :math:`\chi^2`
   * - OLS
     - 112.967525
     - 37.285421
     - 5.933527
     - 0.010898
     - not defined
   * - WLS
     - 112.981000
     - 37.145494
     - 5.976661
     - 0.011033
     - 1.571859
   * - Effective variance
     - 112.980884
     - 37.125926
     - 5.988260
     - 0.011133
     - 0.952062
   * - ODR
     - ``...``
     - ``...``
     - ``...``
     - ``...``
     - ``...``

Then answer two questions:

#. Are the solver-dependent shifts larger than the parameter E.S.D.s?
#. Which uncertainty assumptions best represent the experiment?
