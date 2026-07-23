EOS MGD volume-only fitting through the Python API
==================================================

This example fits a Mie--Grüneisen--Debye P--V--T model using pressure, volume,
and temperature observations only. Adiabatic bulk-modulus observations are not
required and no multi-observable objective is constructed.

:download:`Download mgd_fit_api.py <../_downloads/mgd_fit_api.py>`

Input requirements
------------------

The input table must contain ``P``, ``V``, and ``T`` columns. Effective
variance additionally requires their standard uncertainties.  The distributed
NaF dataset declares ``VSCALE cm^3/mol``; its values are therefore molar
volumes per formula unit.  The MGD normalization must match that convention:

.. code-block:: python

   from quantas.api import eos

   normalization = eos.MGDNormalization.molar_formula_unit(
       formula="NaF",
   )

Do not convert this dataset to a crystallographic cell basis by supplying an
arbitrary ``Z`` value.  Cell-volume normalization is appropriate only when the
input volume is actually reported per crystallographic cell.

Build the compositional model
-----------------------------

.. code-block:: python

   model = eos.PVTModel(
       pressure_model="BM4",
       coupling="thermal-pressure",
       thermal_pressure_model="mgd",
       mgd_normalization=normalization,
   )

Use ``thermal_pressure_model="mgd:q-compromise"`` to select the explicit
q-compromise variant. That variant has no refinable ``q`` parameter.

Choose refinable MGD parameters
--------------------------------

The reference temperature is normally fixed to the temperature of the
reference isotherm. ``theta_d0``, ``gamma0``, and ``q`` may be fixed or refined
through the ordinary :class:`quantas.api.eos.ParameterConstraint`
contract. The example releases all three full-MGD parameters:

.. code-block:: python

   constraints = (
       eos.ParameterConstraint.fixed("temperature_ref", 295.0),
       eos.ParameterConstraint.free(
           "theta_d0",
           459.0,
           lower_bound=1.0,
           upper_bound=5000.0,
       ),
       eos.ParameterConstraint.free(
           "gamma0",
           1.5,
           lower_bound=0.01,
           upper_bound=10.0,
       ),
       eos.ParameterConstraint.free(
           "q",
           1.0,
           lower_bound=-10.0,
           upper_bound=10.0,
       ),
   )

For limited datasets, fixing ``theta_d0`` or ``q`` is often more defensible
than releasing every parameter. Correlation coefficients should always be
inspected before interpreting the fitted values.

Run the fit
-----------

.. code-block:: python

   dataset = eos.read_input("NaF_PVT.dat")
   request = eos.FitRequest(
       model=model,
       domain="pvt",
       target="volume",
       constraints=constraints,
       options=eos.FitOptions(
           solver_options=eos.EffectiveVarianceOptions(
               max_iterations=50,
               inner_max_iterations=10000,
           )
       ),
   )
   result = eos.fit(dataset, request)

The complete downloadable script prints the parameter state, fitted value,
E.S.D., residual statistics, reduced chi-square, and strongest parameter
correlations. Change ``SOLVER`` to ``"ols"`` when the input does not contain
complete uncertainty columns.

.. code-block:: console

   python mgd_fit_api.py NaF_PVT.dat

The same scientific request can be expressed through ``quantas eos run`` or a
``QUANTAS EOS SPEC 1`` file. All routes use the same evaluator and fitting
contracts.
