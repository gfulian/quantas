Thermoelasticity plotting API
=============================

These functions build frontend-neutral plot specifications.  They do not
import or configure Matplotlib directly.  Pass the returned
:class:`quantas.api.plotting.PlotCollection` to
:func:`quantas.api.rendering.render_plots`.

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
