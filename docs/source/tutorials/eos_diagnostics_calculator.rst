EOS diagnostics and calculator
==============================

This tutorial starts from an existing native EOS archive produced by
:command:`quantas eos run` or :class:`quantas.api.eos.Session`. The
post-fit tools never refit the data. They resolve one immutable successful
record, reconstruct its model and parameter covariance, and derive diagnostics
or properties from that record.

Selecting a result
------------------

When an archive contains exactly one accepted result, no explicit selector is
required. With several accepted targets, use a stable result slot such as
``pv/volume``, ``pv/a``, ``vt/volume``, or ``pvt/volume``. An immutable
``record_id`` may be selected instead when a historical candidate or superseded
fit must be inspected.

Residual and finite-strain diagnostics
--------------------------------------

The command

.. code-block:: console

   quantas eos diagnose quartz_EOS.hdf5 --slot pv/volume \
      --output quartz_diagnostics.csv

reports observed and calculated values, physical residuals, standardized
residuals when the solver supplied them, input selection/group metadata, and
finite-strain diagnostics where a scientifically established transformation is
available.

Normalized-pressure representations are provided for Birch--Murnaghan,
natural-strain/Poirier--Tarantola, and Vinet pressure EOS families. Murnaghan
and Tait still receive residual diagnostics, but Quantas does not invent a
generic normalized-pressure transformation for them. At the exact reference
state the normalized pressure is mathematically singular and is therefore
reported as unavailable rather than replaced by a fitted modulus.

For an axial EOS, Quantas applies the same cubed-length convention used during
the fit. The finite strain is evaluated from :math:`L^3/L_0^3`; consequently a
normalized-pressure intercept corresponds to the auxiliary modulus
:math:`M/3`, not directly to the reported linear modulus :math:`M`.

Property calculation
--------------------

Evaluate a pressure grid from an accepted P--V record with

.. code-block:: console

   quantas eos calculate quartz_EOS.hdf5 --slot pv/volume \
      --pressure-range 0:10:1 --output quartz_properties.csv

For a V--T result, provide temperatures:

.. code-block:: console

   quantas eos calculate rutile_EOS.hdf5 --slot vt/volume \
      --temperature-range 20:1200:20

For P--V--T, pressure/volume and temperature ranges form a Cartesian grid by
default:

.. code-block:: console

   quantas eos calculate spessartine_EOS.hdf5 --slot pvt/volume \
      --pressure-range 0:15:1 --temperature-range 300:1200:100 \
      --output spessartine_grid.csv

Use ``--pairwise`` to pair equal-length coordinate and temperature vectors
instead. ``--coordinate`` means volume for a volumetric record and the fitted
physical length for an axial P--V record.

The calculator returns the properties available for the selected domain:

* P--V: volume or axis, modulus, and its first and second pressure derivatives;
* V--T: structural value, thermal-expansion coefficient, and temperature
  derivative;
* P--V--T: pressure, temperature, volume, :math:`K_T`, :math:`K'_T`,
  :math:`K''_T`, :math:`\alpha(P,T)`, :math:`(\partial K/\partial T)_P`,
  zero-pressure reference properties, and thermal pressure when the selected
  coupling defines one.

Uncertainties and extrapolation
-------------------------------

By default, one-sigma property uncertainties are propagated from the complete
fitted parameter covariance through a first-order numerical delta method.
Values supplied as independent state coordinates are not assigned parameter
uncertainties. Use ``--no-uncertainty`` to skip propagation or
``--relative-step`` to control the parameter perturbation used by the delta
method.

These uncertainties represent fitted-parameter covariance only. They do not
include uncertainty in the requested pressure, volume, or temperature, model
selection, systematic experimental errors, or extrapolation. States outside
the sampled ranges of the observations used by the fit are flagged explicitly.

Python API
----------

The downloadable example performs the same workflow through the public API and
writes both CSV files:

:download:`Download eos_diagnostics_calculator_api.py <../_downloads/eos_diagnostics_calculator_api.py>`

.. literalinclude:: ../_downloads/eos_diagnostics_calculator_api.py
   :language: python
   :linenos:
