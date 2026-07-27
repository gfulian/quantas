Thermoelasticity plotting API
=============================

These functions build frontend-neutral plot specifications.  They do not
import or configure Matplotlib directly.  Pass the returned
:class:`quantas.api.plotting.PlotCollection` to
:func:`quantas.api.rendering.render_plots`.


Result-aware discovery
----------------------

.. autofunction:: quantas.api.thermoelasticity.describe_plots

The returned :class:`quantas.api.plotting.PlotInventory` describes cumulative
workflow capabilities rather than assigning one exclusive result type.  A
thermoelastic archive may expose any compatible subset of:

``fit``
   Independent elastic-volume calibration fits and residual diagnostics.
``pt``
   Isothermal or adiabatic stiffness values and uncertainty quantities on a
   two-dimensional pressure-temperature grid.
``profile``
   Absolute or reference-relative stiffness along archived depth paths.
``compare``
   Isothermal--adiabatic sections at fixed pressure or temperature when both
   tensor fields are valid.
``domain``
   Equilibrium-volume coverage, extrapolation masks, and optional profile paths.

Point and one-dimensional analysis archives do not advertise contour families.
The public default builder falls back to the archived calibration fits instead
of attempting to construct an invalid P--T contour.  Comparison coordinates
are evaluated through the calibrated public analysis engine; they are not
nearest-grid selections.

Component resolution
--------------------

.. autofunction:: quantas.api.thermoelasticity.resolve_components

P--T maps
---------

.. autofunction:: quantas.api.thermoelasticity.build_pt_plots

Profile plots
-------------

.. autofunction:: quantas.api.thermoelasticity.build_profile_plots

Isothermal--adiabatic comparison
--------------------------------

.. autofunction:: quantas.api.thermoelasticity.build_compare_plots

Domain diagnostics
------------------

.. autofunction:: quantas.api.thermoelasticity.build_domain_plot
