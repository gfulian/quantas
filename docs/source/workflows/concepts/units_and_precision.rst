Units and numerical precision
=============================

Quantas distinguishes between **stored numerical precision** and **display
precision**.

Stored values
-------------

Native calculations and HDF5 persistence use:

- ``float64`` for real values;
- ``complex128`` for complex values.

There is no runtime precision selector for scientific calculations.

Displayed values
----------------

Formatting and display precision belong to report and plot renderers.  They do
not change the stored data.

Workflow pages document the scientific units relevant to each module.  Where
possible, units are also stored explicitly in metadata associated with results,
reports, or exported tables.
