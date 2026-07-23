Tabular exports, reports, and plot files
========================================

Quantas distinguishes the numerical result from the files used to read,
exchange, or visualize that result.  Native HDF5 is the authoritative
machine-readable representation.  Reports, delimited tables, tensor text
files, and figures are derived views produced for a specific purpose.

This distinction matters because exported tables may contain only selected
properties, may use converted units, and use finite display precision.  They
must not be treated as lossless replacements for the native archive.

Common text conventions
-----------------------

All text files written by Quantas use UTF-8 and Unix line endings.  Reports and
plain-text exports are deterministic and contain no ANSI escape sequences, box
drawing, terminal cursor control, or live-progress artefacts.  They therefore
remain suitable for redirection, version control, and archival use.

Unless a format states otherwise:

* numerical calculations and native HDF5 arrays remain ``float64``;
* text formatting changes only the written representation;
* units are written in headers, column names, or leading metadata comments;
* ``NaN`` or an empty field means that a value is unavailable, not zero;
* Boolean masks retain their semantic names, such as ``extrapolated`` or
  ``mechanically_stable``;
* the native HDF5 archive retains the complete provenance and should accompany
  any table used in a publication or downstream workflow.

Plain-text reports
------------------

A report is a human-readable summary of one calculation or archived record.  It
normally contains the resolved options, scientific summary, warnings,
diagnostics, and selected numerical tables.  Reports deliberately prioritize
readability over completeness.

The following rules apply:

* report precision is controlled by the renderer and never changes stored
  values;
* a report may omit large arrays already present in HDF5;
* operational progress events are not written into the report history;
* warnings, errors, meaningful workflow messages, and structured results are
  retained when relevant;
* rerendering the same result with another display profile may change spacing
  or decimal places without changing the scientific result.

For reproducible work, cite or archive the report together with the HDF5 result,
input, and the exact Quantas version.

HA property tables
------------------

``quantas ha export`` writes one selected HA property to a tab-separated
``.dat`` file.  The header records the property, its description, and the
selected unit.  The first column is temperature and every following column
corresponds to one sampled volume::

   # Quantas HA table export
   # Property: Cv - isochoric heat capacity
   # Units: J mol^-1 K^-1
   #
   T / K    V=15.4954 / angstrom^3    V=16.0123 / angstrom^3 ...

For a temperature-dependent property, the data matrix has one row per
requested temperature and one column per volume.  Temperature-independent
properties, notably the zero-point energy, are written as one static row with
``T = 0``.  This does not alter the native zero-point-energy shape
``(1, n_volume)``.

Supported unit conversion applies to energy-like quantities, entropy, and heat
capacity.  Conversion between per-cell and molar units uses the stored
``formula_units_per_cell`` normalization, so one mole means one mole of formula
units.  Volume and temperature are written in their stored units.

QHA tables
----------

``quantas qha export`` can write a human-oriented table or a machine-readable
CSV/TSV table.  The delimited representation is a flat pressure-temperature
table with one row per state.  The row order is deterministic:

1. pressure varies in the outer loop;
2. temperature varies in the inner loop.

The first columns are therefore::

   P (GPa),T (K),...

Selected properties follow with a physical symbol and unit, for example
``V (angstrom^3)``, ``K_T (GPa)``, or ``alphaV (K^-1)``.  When propagated
uncertainties are available, the corresponding ``sigma_<symbol>`` column is
written immediately after the property.  A thermal-expansion export can also
include ``alphaV method``, which records the source actually used at each state,
such as mixed derivative, mode Gruneisen, or numerical fallback.

When structural reconstruction is available, the same row may contain:

* ``a``, ``b``, and ``c`` in angstrom;
* ``alpha``, ``beta``, and ``gamma`` in degrees;
* axial thermal-expansion coefficients;
* the Cartesian thermal-expansion tensor in engineering Voigt order;
* ``trace(alpha)`` and ``trace(alpha)-alphaV`` consistency fields;
* ``structural_extrapolated``.

The human-readable format groups properties by pressure and scientific family.
It is useful for inspection, while CSV/TSV should be preferred for custom
analysis.  The HDF5 payload remains necessary when failed-state records, fit
diagnostics, covariance, or complete structural arrays are required.

Elasticity directional export
-----------------------------

``quantas elasticity export`` writes principal-plane directional data to a
neutral ``.dat`` file.  Each available plane (``xy``, ``xz``, and ``yz``) is a
separate block beginning with the plane name.  Every row contains
``theta_deg`` and ``phi_deg`` followed by the calculated directional fields.

One-dimensional properties occupy one column.  Multi-branch properties are
expanded into numbered columns.  In the current schema this means, for example:

* linear compressibility: positive and negative-magnitude branches;
* shear modulus: minimum and maximum transverse branches;
* Poisson ratio: negative minimum, positive minimum, and maximum branches.

The numbered branch columns preserve the order stored by the elasticity result;
the report and :doc:`../workflows/elasticity` page provide the physical branch
interpretation.  The export contains only the requested 2D planes.  Global
extrema, full compliance and stiffness matrices, VRH averages, and optional 3D
surfaces remain in the native HDF5 payload.

SEISMIC long-form CSV
---------------------

``quantas seismic export`` writes a long-form CSV with one row for every
sampled wave-normal direction and acoustic mode.  The mode order is explicitly
recorded as ``V_S2, V_S1, V_P``.

The phase-level columns include:

* point index and mode symbol;
* polar and azimuthal angles in radians and degrees;
* Cartesian wave normal;
* Christoffel eigenvalue and phase speed;
* local polarization axis;
* absolute and relative eigenvalue gaps;
* validity, clamping, and degeneracy masks.

A result calculated at ``group`` level additionally exports group velocity,
group speed, ray direction, power-flow angle, and validity/resolution masks.  An
``enhancement`` result adds area factor, logarithmic enhancement, finite and
resolved masks, and the caustic-candidate flag.  When polarization tracking is
available, tracked polarization axes and tracking-resolution flags are appended.

The file begins with comment rows giving the job name, density, and mode order,
followed by the CSV header.  Because the table is long-form, it can be filtered
by ``mode`` without reshaping the original spherical grid.  Grid topology,
tracking branch identities, derivatives, and complete metadata remain more
naturally represented in HDF5.

Thermoelastic point exports
---------------------------

``quantas thermoelasticity analysis point`` can write a compact text file for a
single reconstructed state.  The file contains pressure, temperature, volume,
density, one selected 6 x 6 stiffness tensor, and the scientific state flags.
When written in the elasticity/SEISMIC input convention, it can be passed
directly to::

   quantas elasticity run state.dat
   quantas seismic run state.dat

The tensor condition must be recorded and respected.  An isothermal tensor
``C_T`` and an adiabatic tensor ``C_S`` are not interchangeable, particularly
for comparison with acoustic measurements.

Thermoelastic grid and profile tables
-------------------------------------

Thermoelastic grid and profile exports use a wide representation: one row is one
physical state, and symmetry-independent stiffness components occupy separate
columns.  This is the recommended interface for custom plots and statistical
post-processing.

Grid rows begin with::

   pressure_GPa,temperature_K,volume_A3,density_kg_m3

Profile rows add::

   profile,depth_km

before the thermodynamic coordinates.  Depending on the requested tensor
condition, component columns use suffixes such as::

   C11_T_GPa
   C11_S_GPa
   sigma_C11_T_GPa
   sigma_C11_S_GPa

The table also contains, when available:

* ``qha_extrapolated``;
* ``elastic_extrapolated``;
* ``adiabatic_valid``;
* ``mechanically_stable``;
* ``minimum_stiffness_eigenvalue_GPa``.

CSV uses commas; deterministic text tables use tabs.  The profile row order
follows increasing archived depth.  Grid row order is temperature-major and
then pressure-major.  Empty fields denote diagnostics that were not available
for that archive.

Complete thermoelastic tensor export
------------------------------------

The separate tensor exporter writes every selected 6 x 6 tensor from a grid,
profiles, or both.  It supports deterministic text and CSV.

The text representation uses one state block containing:

* state identifier and origin (grid or profile);
* tensor condition and validity;
* depth, pressure, temperature, volume, and density;
* QHA and elastic extrapolation flags;
* mechanical-stability diagnostics;
* the full stiffness matrix in GPa;
* an optional full one-sigma matrix.

The CSV representation flattens the upper Voigt triangle into named columns and
retains the same state metadata.  Use this export when all tensor elements are
needed but HDF5 is inconvenient.  Use the compact wide table instead when only
symmetry-independent components are required.

EOS diagnostic and calculation CSV
----------------------------------

EOS has two post-fit CSV families.

``diagnose`` exports one row per observation and can include the original
coordinates, calculated response, physical and standardized residuals,
inclusion mask, group, finite strain, normalized pressure, and validity fields.
The exact columns depend on the domain and selected model.

``calculate`` exports one row per requested state.  Every calculated property is
written under its stable machine-readable name; a propagated uncertainty is
written as ``sigma_<name>`` immediately after the property.  Units use square
brackets in the header, for example::

   pressure [GPa],volume [angstrom^3],sigma_volume [angstrom^3]

EOS CSV values use up to 16 significant digits.  A non-finite unavailable value
is written as an empty field.  These files are views of one immutable EOS archive
record; they do not contain the full record history, acceptance state, resolved
batch specification, or parameter covariance.

Plot files
----------

Plots are derived visualizations generated from a result or archived record.
Quantas may write formats supported by the selected renderer, commonly PNG,
SVG, or PDF.  Figure preset, size, DPI, colormap, contour levels, marker style,
and polarization stride affect presentation only.

A plot should never be the sole preserved scientific output.  In particular:

* raster pixels do not retain numerical precision;
* contour interpolation does not add physical sampling points;
* a smooth EOS curve does not demonstrate data support outside observations;
* changing DPI or the number of graphical curve points does not rerun the
  scientific calculation.

Recommended archival set
-------------------------

For a calculation that must remain reproducible, preserve at least:

1. the original input or batch/profile specification;
2. the native Quantas HDF5 result or EOS archive;
3. the plain-text report;
4. the Quantas version and execution options;
5. any exported table or figure actually used downstream.

For direct programmatic access to native arrays and metadata, continue with
:doc:`hdf5_inspection`.
