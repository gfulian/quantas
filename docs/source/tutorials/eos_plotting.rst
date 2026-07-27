EOS plotting
============

EOS plots are generated from immutable HDF5 fit records. Plot generation is a
separate post-processing command: :command:`quantas eos run` never imports
Matplotlib and never creates figures as a side effect.

CLI workflow
------------

Generate every plot supported by the accepted result:

.. code-block:: console

   quantas eos plot quartz_EOS.hdf5 --slot pv/volume

The default output directory is ``quartz_EOS_plots``. A smaller explicit set
can be selected by repeating ``--plot``:

.. code-block:: console

   quantas eos plot quartz_EOS.hdf5 --slot pv/volume \
       --plot fit --plot residuals --plot normalized-pressure \
       --uncertainties --point-size 5 --curve-width 2.0 \
       --no-title --figure-width 7.2 --figure-height 5.1 \
       --axis-label-font-size 16 --legend-font-size 14 \
       --title-font-size 16 --tick-label-font-size 12 \
       --format png --dpi 200

Use ``--list-plots`` to inspect the plots available for a particular record.
The command uses the same public session-aware inventory consumed by Python
clients and Quantas GUI.  When an archive contains no accepted result, or more
than one accepted result slot, the command lists the available slots and
plottable immutable records instead of choosing one silently.  Add ``--slot``
or ``--record-id`` to request the detailed representation inventory.

The list depends on the scientific domain and on the information actually
stored by the solver. Standardized residuals, for example, are offered only
when they exist.

P--V plots include observed data and the fitted curve, residuals, standardized
residuals when available, and finite-strain normalized-pressure diagrams for
Birch--Murnaghan, Natural Strain, and Vinet models. V--T plots include the
fitted thermal-expansion curve and residuals. P--V--T records provide the
sampled pressure-temperature coverage, calculated isotherms, calculated
isobars, and residuals against both pressure and temperature.

Input uncertainty bars can be hidden with ``--no-uncertainties``. Excluded
observations remain visible by default with a distinct marker and can be hidden
with ``--no-excluded``. Input ``GROUP`` identifiers are rendered as separate
series unless ``--no-group-data`` is used. Figure titles are enabled by
default and can be removed with ``--no-title`` for publication layouts that
use external captions. ``--figure-width`` and ``--figure-height`` control the
output size in inches. Axis-label, legend, title, and tick-label sizes are
controlled independently with ``--axis-label-font-size``,
``--legend-font-size``, ``--title-font-size``, and
``--tick-label-font-size``. The first three default to 16, 14, and 16 points;
tick labels preserve the Matplotlib default unless explicitly overridden. The
Matplotlib renderer applies a tight layout after resizing and typography
updates and before saving.

In normalized-pressure plots, finite observations whose measured pressure is
statistically indistinguishable from zero are hidden by default. Such points
carry little information about EOS order and can force the horizontal range
towards the origin because ``F`` contains division by finite strain. This is a
plot-only choice: the observation remains in the fit, diagnostics table, and
HDF5 archive. Use ``--zero-pressure-point`` to display it when desired.

For P--V--T plots, explicit curves can be requested:

.. code-block:: console

   quantas eos plot spessartine_EOS.hdf5 --slot pvt/volume \
       --plot isotherms --isotherm 300 --isotherm 600 --isotherm 900 \
       --plot isobars --isobar 0 --isobar 5 --isobar 10

If no values are supplied, Quantas chooses a compact representative set over
the sampled pressure or temperature range.

Python API
----------

The public inventory keeps the persistent EOS lifecycle explicit.  Archive
datasets, result slots, immutable records, acceptance state, and representation
keys can be inspected without opening a writable session or building figures:

.. code-block:: python

   from quantas.api import eos

   inventory = eos.describe_plots(
       "quartz_EOS.hdf5",
       slot="pv/volume",
   )
   for dataset in inventory.datasets:
       print(dataset.dataset_id, dataset.jobname, dataset.slot_keys)
   for record in inventory.records:
       print(record.record_id, record.disposition.value, record.representation_keys)

   plots = inventory.selected_plots
   if plots is not None:
       for representation in plots.representations:
           print(representation.key, representation.description)

The selected representation keys are accepted directly by
:func:`quantas.api.eos.build_plots`.  The API first prepares frontend-neutral
:class:`~quantas.api.plotting.PlotCollection` objects. Matplotlib is invoked only
through the public rendering bridge:

.. code-block:: python

   from quantas.api import eos, rendering

   collection = eos.build_plots(
       "quartz_EOS.hdf5",
       ("fit", "residuals", "normalized-pressure"),
       slot="pv/volume",
       options=eos.PlotOptions(
           show_uncertainties=True,
           point_size=5.0,
           curve_width=2.0,
           show_title=False,
           show_zero_pressure_point=False,
       ),
   )
   rendered = rendering.render_plots(
       collection,
       output_dir="quartz_plots",
       image_format="png",
       preset="publication",
       dpi=200,
       figure_size=(7.2, 5.1),
       close=True,
   )

   for path in rendered.paths:
       print(path)

The same neutral specifications are consumed by the separate Quantas GUI
package and may be reused by other renderers.

Downloadable example
--------------------

:download:`Download eos_plot_api.py <../_downloads/eos_plot_api.py>`

.. literalinclude:: ../_downloads/eos_plot_api.py
   :language: python
