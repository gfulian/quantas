Rendering API
=============

Scientific modules return neutral :class:`quantas.api.common.ReportTable` and
:class:`quantas.api.plotting.PlotCollection` objects.  The supported functions in
:mod:`quantas.api.rendering` turn those contracts into deterministic text and
static figures without exposing concrete renderer classes.

Plain-text rendering is available in the base installation.  Static figure
rendering requires the ``plot`` extra.

.. code-block:: python

   from pathlib import Path
   from quantas.api import elasticity, rendering

   result = elasticity.run("calcite.dat")
   text = rendering.render_tables(elasticity.build_report(result))
   Path("calcite.log").write_text(text, encoding="utf-8")

   rendered = rendering.render_plots(
       elasticity.build_3d_plots(result),
       output_dir="figures",
       preset="publication",
       close=True,
   )

Rendered plot contracts
-----------------------

.. autoclass:: quantas.api.rendering.RenderedPlot
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.rendering.PlotRenderResult
   :members:
   :show-inheritance:

Table rendering
---------------

.. autofunction:: quantas.api.rendering.render_table

.. autofunction:: quantas.api.rendering.render_tables

Plot rendering
--------------

.. autofunction:: quantas.api.rendering.render_plots

See also
--------

- :doc:`../getting_started/results`
- :doc:`../formats/tabular_outputs`
- :doc:`../workflows/concepts/reporting_and_events`
