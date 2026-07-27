Frontend-neutral plotting contracts
===================================

:mod:`quantas.api.plotting` exposes the passive contracts produced by public
scientific plot builders.  These objects carry prepared numerical values,
scientific axes, units, masks, overlays, uncertainty information, directional
fields, and portable presentation hints without constructing a concrete
figure.

Renderer independence
---------------------

The contracts do not depend on Matplotlib, Plotly, Dash, Rich, or browser
components.  A frontend may translate them into a static figure, an
interactive figure, a notebook object, or another representation while
preserving the scientific data and conventions supplied by Quantas.

Use the public namespace when implementing a renderer:

.. code-block:: python

   from quantas.api import elasticity, plotting

   result = elasticity.run(
       "calcite.dat",
       options=elasticity.Options(calculate_2d=True),
   )
   collection = elasticity.build_2d_plots(result, properties=("young",))

   for specification in collection.plots:
       if isinstance(specification, plotting.PolarPlotSpec):
           print(specification.key, len(specification.panels))

``PlotCollection`` remains available from :mod:`quantas.api.common` for
compatibility with the original common-contract entry point.  Both aliases
refer to the same class.  New renderer implementations should normally import
plot-specific contracts from this namespace.

Axis and Cartesian primitives
-----------------------------

.. autodata:: quantas.api.plotting.AxisLocation

.. autodata:: quantas.api.plotting.AxisOrientation

.. autodata:: quantas.api.plotting.BandOrientation

.. autodata:: quantas.api.plotting.LineStyle

.. autoclass:: quantas.api.plotting.PlotAxis
   :members:

.. autoclass:: quantas.api.plotting.PlotSeriesStyle
   :members:

.. autoclass:: quantas.api.plotting.PlotSeries
   :members:

.. autoclass:: quantas.api.plotting.PlotBandStyle
   :members:

.. autoclass:: quantas.api.plotting.PlotBand
   :members:

.. autoclass:: quantas.api.plotting.ColoredPathStyle
   :members:

.. autoclass:: quantas.api.plotting.ColoredPathSeries
   :members:

.. autoclass:: quantas.api.plotting.SecondaryAxis
   :members:

.. autoclass:: quantas.api.plotting.PlotSpan
   :members:

.. autoclass:: quantas.api.plotting.PlotMask
   :members:

.. autoclass:: quantas.api.plotting.ScalarBackground
   :members:

Cartesian specifications
------------------------

.. autoclass:: quantas.api.plotting.LinePlotSpec
   :members:

.. autoclass:: quantas.api.plotting.ContourPlotSpec
   :members:

Directional and surface specifications
--------------------------------------

.. autoclass:: quantas.api.plotting.PolarPlotPanel
   :members:

.. autoclass:: quantas.api.plotting.PolarPlotSpec
   :members:

.. autoclass:: quantas.api.plotting.SurfaceStyle
   :members:

.. autoclass:: quantas.api.plotting.SurfaceLayer
   :members:

.. autoclass:: quantas.api.plotting.VectorFieldStyle
   :members:

.. autoclass:: quantas.api.plotting.VectorFieldLayer
   :members:

.. autoclass:: quantas.api.plotting.SurfacePlotSpec
   :members:

.. autodata:: quantas.api.plotting.SphericalProjection

.. autoclass:: quantas.api.plotting.SphericalMarker
   :members:

.. autoclass:: quantas.api.plotting.AxisFieldLayer
   :members:

.. autoclass:: quantas.api.plotting.SphericalMapSpec
   :members:

.. autoclass:: quantas.api.plotting.SphericalSummarySpec
   :members:

Composite and collection contracts
----------------------------------

.. autoclass:: quantas.api.plotting.PanelPlotSpec
   :members:

.. autodata:: quantas.api.plotting.PlotSpec

.. autoclass:: quantas.api.plotting.PlotCollection
   :members:

See also
--------

- :doc:`rendering`
- :doc:`../developer/rendering_frontends`
- :doc:`../developer/module_anatomy`
