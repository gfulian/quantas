Inspecting and extracting native HDF5 data
==========================================

Quantas native files are designed to support analysis beyond the built-in
reports and plots.  There are two complementary reading strategies:

#. use :mod:`quantas.api` for schema validation, typed payloads, and historical
   migrations;
#. use h5py for direct access from lightweight scripts or external software.

Use the public API when possible.  Use h5py deliberately, after inspecting the
stored module, schema, paths, units, and validity masks.

Quick inspection from the shell
-------------------------------

Standard HDF5 utilities are useful when installed:

.. code-block:: console

   h5ls -r result.hdf5
   h5dump -A result.hdf5

``h5ls -r`` lists groups and datasets. ``h5dump -A`` reports attributes without
printing all numerical arrays.  These commands do not understand Quantas
schemas; they are structural tools only.

A portable Python inspector is distributed with the documentation:

:download:`Download inspect_hdf5.py <../_downloads/formats/inspect_hdf5.py>`

Run it with:

.. code-block:: console

   python inspect_hdf5.py result.hdf5

The script prints:

- the Quantas module, envelope schema, and writer version;
- every group and dataset path;
- dataset shape and dtype;
- local unit and description attributes;
- compression and chunking;
- scalar attributes attached to groups and datasets.

It does not load large arrays into memory.

Metadata-driven reading with Quantas
------------------------------------

The public registry can identify the responsible module from HDF5 metadata:

.. code-block:: python

   from pathlib import Path
   from quantas.api import registry

   path = Path("result.hdf5")
   descriptor = registry.module_from_result(path)

   print(descriptor.name)
   print(descriptor.title)
   print(descriptor.capabilities)

Open a native result without guessing the module from its filename:

.. code-block:: python

   opened = registry.open_result(path)

For HA, QHA, Elasticity, SEISMIC, and Thermoelasticity, ``opened`` is a
:class:`quantas.api.common.ResultData`.  EOS returns an active archive because
its HDF5 file stores datasets, immutable fit records, and mutable accepted-state
indices rather than one result envelope.

A safe generic dispatcher can therefore be written as:

.. code-block:: python

   from contextlib import closing
   from quantas.api import registry

   descriptor = registry.module_from_result("result.hdf5")
   opened = registry.open_result("result.hdf5")

   if descriptor.name == "eos":
       with closing(opened) as archive:
           print(archive.path)
           # Query datasets, records, and accepted slots through the EOS API.
   else:
       envelope = opened
       print(envelope.metadata)
       print(envelope.options)
       print(envelope.warnings)

Use module-specific readers when the subsequent analysis is known.

Typed HA example
----------------

.. code-block:: python

   import numpy as np
   from quantas.api import ha

   envelope = ha.read_result("mgo_ha.hdf5")
   data = ha.get_result(envelope)

   temperature = np.asarray(data.temperature)
   volume = np.asarray(data.volume)
   cv = np.asarray(data.isochoric_heat_capacity)
   zpe = np.asarray(data.zero_point_energy)

   # Temperature-dependent HA arrays use (temperature, volume).
   print(cv.shape)

   # Zero-point energy is static and retains shape (1, volume).
   print(zpe.shape)

Apply the stored q-point normalization and formula-unit basis only through the
typed result metadata or documented conversion utilities.  Do not infer molar
normalization from array magnitudes.

Typed QHA example
-----------------

.. code-block:: python

   import numpy as np
   from quantas.api import qha

   envelope = qha.read_result("mgo_qha.hdf5")
   data = qha.get_result(envelope)

   temperature = np.asarray(data.temperature)
   pressure = np.asarray(data.pressure)
   volume = np.asarray(data.equilibrium_volume)
   alpha = np.asarray(data.thermal_expansion)

   # QHA grid convention: (temperature, pressure)
   i300 = int(np.argmin(np.abs(temperature - 300.0)))
   i0 = int(np.argmin(np.abs(pressure - 0.0)))

   print(volume[i300, i0])
   print(alpha[i300, i0])
   print(data.valid_mask[i300, i0])

The typed reader reconstructs fit diagnostics, failed points, uncertainties,
metadata, and optional structural arrays.  It also applies supported historical
unit migrations.

Typed Elasticity example
------------------------

.. code-block:: python

   from quantas.api import elasticity

   envelope = elasticity.read_result("calcite_elasticity.hdf5")
   data = elasticity.get_result(envelope)

   print(data.stiffness)
   print(data.compliance)
   print(data.averages.hill.bulk_modulus)
   print(data.stability.is_stable)
   print(data.variations["young_modulus"].maximum)

Three-dimensional surfaces are optional.  Check ``data.properties_3d`` before
using them.

Typed SEISMIC example
---------------------

.. code-block:: python

   import numpy as np
   from quantas.api import seismic

   envelope = seismic.read_result("seismic.hdf5")
   data = seismic.get_result(envelope)

   directions = np.asarray(data.grid.directions)
   phase = np.asarray(data.field.phase.phase_speeds)
   valid = np.asarray(data.field.phase.valid_mask)

   # Last mode index is V_P; mode order is V_S2, V_S1, V_P.
   vp = np.where(valid[..., 2], phase[..., 2], np.nan)
   print(np.nanmin(vp), np.nanmax(vp))

``group`` and ``enhancement`` fields may be absent when a lower sampling level
was requested.  Degeneracy, resolved, finite, and validity masks must be
applied before reducing those arrays.

Typed Thermoelasticity example
------------------------------

.. code-block:: python

   import numpy as np
   from quantas.api import thermoelasticity

   envelope = thermoelasticity.read_result("thermoelastic_grid.hdf5")
   data = thermoelasticity.get_result(envelope)

   temperature = np.asarray(data.temperature)
   pressure = np.asarray(data.pressure)
   c_t = np.asarray(data.stiffness_isothermal)

   # Tensor-grid convention: (temperature, pressure, 6, 6)
   c11 = c_t[..., 0, 0]
   supported = ~np.asarray(data.extrapolation_mask)
   print(np.nanmean(np.where(supported, c11, np.nan)))

For adiabatic analysis, also check ``stiffness_adiabatic`` and
``adiabatic_valid_mask``.  Do not substitute the isothermal tensor at invalid
non-zero-temperature states unless that choice is explicit in the user's own
analysis.

EOS archive example
-------------------

EOS stores a dataset and immutable fit records rather than one
:class:`~quantas.api.common.ResultData` payload.  Open it as a context manager:

.. code-block:: python

   from quantas.api import eos

   with eos.open_archive("quartz_eos.hdf5") as archive:
       summary = archive.inspect()
       print(summary.record_count)
       for slot in summary.slots:
           print(slot.state.slot.key, slot.state.accepted_record_id)

       for record in archive.records():
           print(record.record_id, record.slot.key, record.successful)

       accepted = archive.accepted("pv/volume")
       if accepted is not None:
           print(accepted.result.parameter_values)
           print(accepted.result.fit.covariance)

Use :func:`quantas.api.eos.diagnose` for observation-level residual tables and
:func:`quantas.api.eos.calculate` for derived states.  Directly reading one
``/records`` group without consulting state events can confuse a rejected or
candidate record with the currently accepted result.

Direct h5py inspection
----------------------

A minimal schema-aware script should inspect root metadata before reading a
payload:

.. code-block:: python

   from pathlib import Path
   import h5py

   path = Path("result.hdf5")

   with h5py.File(path, "r") as h5:
       metadata = h5["metadata"].attrs
       module = metadata["module"]
       schema = metadata["schema_version"]
       version = metadata["version"]

       print(module, schema, version)
       print(list(h5.keys()))

h5py may return byte strings or NumPy scalar objects depending on the HDF5
storage type.  Decode bytes and call ``.item()`` on NumPy scalars when a native
Python value is needed.

Reading one dataset safely
--------------------------

.. code-block:: python

   import h5py
   import numpy as np

   with h5py.File("mgo_qha.hdf5", "r") as h5:
       dataset = h5["/results/equilibrium_volume"]
       values = np.asarray(dataset[...], dtype=np.float64)
       unit = dataset.attrs.get("unit")
       description = dataset.attrs.get("description")

       print(values.shape, unit, description)

The absence of a local ``unit`` attribute does not prove that a value is
unitless.  Older payloads and recursively serialized fields may keep units in
``/results/metadata`` or another neighboring metadata group.

Partial reads for large arrays
------------------------------

Do not load a full SEISMIC enhancement field or dense thermoelastic grid when
only one slice is needed.  h5py supports slicing directly on disk:

.. code-block:: python

   import h5py

   with h5py.File("seismic.hdf5", "r") as h5:
       speeds = h5["/results/fields/phase/phase_speeds"]

       # A contiguous block of sampled directions, all three modes.
       block = speeds[1000:2000, :]

       # Only V_P for every sampled direction.
       vp = speeds[:, 2]

SEISMIC arrays are stored with a flattened direction axis.  Read ``ntheta`` and
``nphi`` from ``/results/grid`` attributes before reshaping a scalar mode field
for image-like processing.  For a QHA or Thermoelasticity grid, slice the
temperature and pressure axes before the trailing property or tensor axes.

Validity and support masks
--------------------------

Always search for masks associated with a scientific array.  Common examples
include:

- QHA ``valid_mask`` and structural extrapolation mask;
- SEISMIC phase/group/enhancement validity and resolved masks;
- SEISMIC degeneracy, finite, and caustic-candidate masks;
- Thermoelasticity QHA and elastic extrapolation masks;
- Thermoelasticity adiabatic validity and mechanical-stability masks;
- EOS per-record selection masks and accepted-state references.

A finite value can still be an extrapolation.  A NaN can intentionally preserve
array topology for an invalid state.  A zero can be a legitimate physical value
or an explicit exactly-zero elastic component.  Interpret the value together
with its masks and metadata.

Extracting arrays to NPZ
------------------------

A second distributed script extracts selected dataset paths to a compressed
NumPy archive and writes a JSON manifest with source paths, units, descriptions,
shapes, dtypes, root attributes, and every group and selected-dataset attribute:

:download:`Download extract_hdf5.py <../_downloads/formats/extract_hdf5.py>`

Example:

.. code-block:: console

   python extract_hdf5.py mgo_qha.hdf5 \
      --dataset /results/temperature \
      --dataset /results/pressure \
      --dataset /results/equilibrium_volume \
      --dataset /results/thermal_expansion \
      --output mgo_qha_subset.npz

Read the extracted data with:

.. code-block:: python

   import json
   from pathlib import Path
   import numpy as np

   manifest = json.loads(
       Path("mgo_qha_subset.json").read_text(encoding="utf-8")
   )
   with np.load("mgo_qha_subset.npz", allow_pickle=False) as arrays:
       volume = arrays["results__equilibrium_volume"]
       print(volume.shape)

   print(manifest["datasets"]["results__equilibrium_volume"]["unit"])

``--all`` extracts every dataset, converts variable-length HDF5 text arrays to
portable Unicode arrays, and avoids pickle-dependent NumPy objects.  It may
require substantial memory and may include raw input text, reports, detailed
diagnostics, and large derivative fields.  Prefer an explicit dataset list for
reproducible analysis.  The JSON manifest always records the complete group
attribute tree, including scalar values that are not HDF5 datasets.

Exporting a rectangular table
------------------------------

For a two-dimensional thermodynamic grid, build a tidy table by combining
coordinate axes with selected properties:

.. code-block:: python

   import numpy as np
   import pandas as pd
   from quantas.api import qha

   result = qha.get_result(qha.read_result("mgo_qha.hdf5"))

   temperature, pressure = np.meshgrid(
       result.temperature,
       result.pressure,
       indexing="ij",
   )
   table = pd.DataFrame(
       {
           "temperature_K": temperature.ravel(),
           "pressure_GPa": pressure.ravel(),
           "volume_A3": result.equilibrium_volume.ravel(),
           "KT_GPa": result.isothermal_bulk_modulus.ravel(),
           "alphaV_K-1": result.thermal_expansion.ravel(),
           "valid": result.valid_mask.ravel(),
       }
   )
   table.to_csv("qha_analysis.csv", index=False)

Pandas is not required by Quantas; it is used here only as an example analysis
tool.  Record units in column names or separate metadata rather than relying on
memory.

Extracting all recursive metadata
---------------------------------

Attributes and nested mappings are not captured by a simple dataset export.
When Quantas is installed, the module reader already reconstructs them into:

- ``envelope.metadata``;
- ``envelope.input_data``;
- ``envelope.options``;
- ``envelope.warnings``;
- ``envelope.events``;
- payload-specific ``metadata`` and diagnostics.

For a pure-h5py reader, recursively walk both ``group.attrs`` and
``group.items()``.  The distributed inspector demonstrates the pattern.

Provenance for derived analyses
-------------------------------

A derived file or publication table should record at least:

- source HDF5 filename and SHA-256 checksum;
- Quantas version;
- module and schema version;
- exact HDF5 dataset paths;
- original units and any conversions;
- slices, masks, and exclusions applied;
- equations or transformations used by the external script;
- script version or commit identifier.

Do not overwrite the original native HDF5 file.  Derived results should be new
artifacts whose provenance points back to the immutable source.

Common mistakes
---------------

- Selecting a reader from the filename instead of ``/metadata:module``.
- Reading a large compressed dataset with ``[...]`` when a slice is sufficient.
- Ignoring attributes because only datasets were listed.
- Assuming all scientific arrays use the same axis order.
- Ignoring validity, degeneracy, extrapolation, or acceptance masks.
- Treating display-rounded reports as the numerical source.
- Replacing NaNs or absent optional fields by zero.
- Editing a native file in place and breaking internal provenance.
- Comparing values without reading their normalization basis and units.
