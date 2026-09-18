``quantas eos``
===============

EOS differs intentionally from the other scientific groups.  One dataset is
often fitted repeatedly with different formulations, solvers, selections,
constraints, and initial values.  Quantas therefore stores immutable fit
records in a persistent archive and tracks which record is accepted or marked
as a candidate for each scientific slot.

Model discovery
---------------

Quantas keeps the historical compact EOS tags (for example ``M``, ``BM3``,
``NS4``/``PT4``, ``V3``, ``T3``, and ``SJ``) rather than exposing a separate
family/order option pair in workflows that already consume model tags.  Long
family aliases are accepted by the common resolver and normalized to the
canonical tag used in reports and persisted metadata.

The complete EOS catalogue is available without expanding every command help
page into long choice lists::

   quantas eos show-models
   quantas eos show-models --domain pv
   quantas eos show-models --domain ev
   quantas eos show-models --domain vt
   quantas eos show-models --domain pvt
   quantas eos show-models --domain ev --domain vt

With no filter the command shows the domain-capability matrix and the model
catalogues for P--V, E--V, V--T, and P--V--T.  Repeated ``--domain`` options
select the requested sections; they do not require one model to span unrelated
scientific domains.  The E--V catalogue includes models whose analytical
``P(V) = -dE/dV`` derivative is available even when direct experimental P--V
fitting is not exposed.

EOS-valued CLI options use the same resolver and provide model-aware shell
completion.  Click 8.1 does not provide a native PowerShell completion backend,
so Quantas supplies one.  Enable it for the current PowerShell session with::

   quantas completion powershell | Out-String | Invoke-Expression

Place the same line in ``$PROFILE`` to enable it in future sessions.  Normal
``--help`` output therefore remains compact while TAB completion can expose
canonical tags and established aliases such as ``NS`` and ``SJEOS``.

Energy EOS input generation
---------------------------

``quantas eos inpgen`` collects structure--energy states from electronic-structure
outputs without selecting an EOS formulation or fit settings.  CRYSTAL output
files and independent VASP calculation directories use the same command::

   quantas eos inpgen crystal.out --interface crystal -o energy.dat
   quantas eos inpgen files.txt --interface crystal --list -o energy.dat
   quantas eos inpgen vasp-runs.txt --interface vasp --list -o energy.dat

A CRYSTAL source may contain one static state, one completed geometry
optimization, or a native multi-volume ``EOS`` calculation.  With ``--list``,
each listed output may contribute one or several states.  A VASP source is one
calculation directory (or its ``vasprun.xml``/``OUTCAR`` path) and contributes
exactly one completed ionic state; an optimization history is rejected for this
workflow.  ``vasprun.xml`` is the primary VASP record and a sibling ``OUTCAR``
is used as a consistency check when available.

VASP E--V collection selects ``energy(sigma->0)`` and requires a homogeneous
electronic-energy signature across all listed runs.  This prevents an
optimization performed with one smearing scheme from being inserted silently
into a static EOS series calculated with another scheme.  VASP source cells are
normalized to the primitive cell before collection, with source energy scaled
by the same integer cell multiplicity.  The resulting points from either
backend are checked for compatible composition and energy semantics and sorted
by volume.

The generated table contains ``V A B C ALPHA BETA GAMMA E`` and explicit units.
It also records ``SYSTEM``, the reference space group, ``CRYSTAL_REFERENCE``,
and ``CELL_MULTIPLICITY``.  The default reference is ``primitive``.  Use
``--crystal-reference crystallographic`` to apply one constant conventional-cell
transformation and scale energy and volume by the same cell multiplicity.  The
table remains independent of the later model, solver, constraints, and specfile.

Recommended sequence
--------------------

For a static Energy EOS generated from electronic-structure outputs:

.. code-block:: console

   quantas eos inpgen crystal-files.txt --interface crystal --list -o mgo_ev.dat
   # or: quantas eos inpgen vasp-runs.txt --interface vasp --list -o mgo_ev.dat
   quantas eos run mgo_ev.dat --domain ev --ev-eos BM3 -o mgo_ev.hdf5
   # Optional pressure-form parameterization of the derived axial path:
   quantas eos run mgo_ev.dat --domain ev --ev-eos SJ --axial-eos BM3 \
      -o mgo_ev_axial.hdf5 --force
   quantas eos diagnose mgo_ev.hdf5 --slot ev/energy
   quantas eos calculate mgo_ev.hdf5 --slot ev/energy --pressure-range 0:20:2
   quantas eos plot mgo_ev.hdf5 --slot ev/energy

For one homogeneous P--V fit:

.. code-block:: console

   quantas eos run PV_quartz.dat --domain pv --pv-eos BM --pv-order 3 \
      --solver effective-variance
   quantas eos diagnose PV_quartz_EOS.hdf5 --slot pv/volume
   quantas eos plot PV_quartz_EOS.hdf5 --slot pv/volume
   quantas eos calculate PV_quartz_EOS.hdf5 --slot pv/volume \
      --pressure-range 0 10 1

For a heterogeneous reproducible batch:

.. code-block:: console

   quantas eos spec-template quartz.spec
   quantas eos run PV_quartz.dat --spec quartz.spec --dry-run
   quantas eos run PV_quartz.dat --spec quartz.spec

Why ``run`` is different
------------------------

``run`` is a batch fitting controller rather than a single-result calculator.
It parses and normalizes the dataset, resolves jobs, validates model and solver
requirements, persists every attempted record, and applies explicit acceptance
rules.  A failed job may therefore be present in the archive even though the
command exits unsuccessfully.

``diagnose``, ``plot``, and ``calculate`` select an accepted slot or an explicit
immutable record.  They do not refit the data.  The latest record is not
necessarily the accepted record.

Portable text notation
----------------------

Terminal output and deterministic plain-text reports use ASCII scientific
notation deliberately.  Examples include ``angstrom^3``, ``GPa^-1``, ``K'``,
``K''``, ``alpha0``, and ``eta_a``.  This convention avoids dependence on
terminal fonts or Unicode rendering and does not alter machine-readable names,
stored units, numerical values, or the typographic notation used in plots and
scientific documentation.

Important option families
-------------------------

EOS already follows the self-describing-input rule for unit handling.  Its
``--eunit``, ``--lunit``, ``--punit``, and ``--tunit`` options are input-unit
overrides: when omitted, declarations in the EOS data file are used before the
historical fallbacks.  This differs from HA/QHA, where pressure and temperature
primarily describe the requested calculation/output domain.

* domain and target options define the scientific slot;
* E--V, P--V, V--T, and P--V--T options define the model rather than solver behavior;
* ``--axial-eos`` is an optional secondary pressure-form parameterization for
  an E--V structural path; it does not replace the primary Energy EOS;
* MGD normalization is part of the physical model and must match the volume
  basis;
* ``--fix``, ``--initial``, and ``--bound`` have different scientific meanings;
* OLS, WLS, effective variance, and ODR use different uncertainty assumptions;
* covariance scaling affects parameter uncertainties, not the optimum itself;
* ``--dry-run`` validates the resolved plan without fitting or creating HDF5;
* ``--failure-policy=continue`` continues later independent jobs but does not
  make the batch successful.

See :doc:`../workflows/eos` for implementation and record semantics,
:doc:`../tutorials/eos` for the E--V/P--V/V--T/P--V--T course, and
:doc:`../formats/eos_spec` plus :doc:`../formats/eos_hdf5` for persistence.

Generated command reference
---------------------------

.. click:: quantas.cli.eos:eos
   :prog: quantas eos
   :nested: full
