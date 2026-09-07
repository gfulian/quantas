``quantas ha``
==============

The HA frontend evaluates vibrational thermodynamics independently at every
sampled volume.  It does not minimize the free energy with respect to volume
and does not create a pressure axis.

Recommended sequence
--------------------

.. code-block:: console

   quantas ha inpgen phonon-output.out --output material.yaml
   quantas ha add-kieffer material.yaml elastic.out \
      --output material-kieffer.yaml
   quantas ha run material-kieffer.yaml --kieffer --temperature 0 1000 10
   quantas ha plot material-kieffer_HA.hdf5 --property Cv --property Fvib
   quantas ha plot material-kieffer_HA.hdf5 --property F --axis volume \
      --temperature 300 --temperature 1000 --2d
   quantas ha export material-kieffer_HA.hdf5 --property Cv --unit J/mol

``inpgen`` is optional when a valid Quantas phonon YAML already exists.  ``run``
creates the scientific archive and report.  ``plot`` and ``export`` read that
archive without repeating harmonic sums.

Generating the phonon input
---------------------------

The shared ``inpgen`` command supports three interface routes:

.. code-block:: console

   quantas ha inpgen phonon.out --interface crystal --output material.yaml

   quantas ha inpgen files.txt --list --interface crystal \
      --reference 0 --output material.yaml

   quantas ha inpgen qha.out --interface crystal-qha --output material.yaml

``--reference`` selects the source structure used for reference metadata and
final branch labels in a multi-file series.  It does not change the local
adjacent-volume overlap assignments performed by the mode tracker.

``--formula-units`` records the number of chemical formula units represented by
the normalization cell.  This value participates in later molar conversions
and should not be chosen from the phonon supercell size alone.

``--debug`` prints mode-by-mode continuity diagnostics for multi-volume inputs:
raw source modes, frequencies, selected and competing overlaps, overlap gap,
degenerate-subspace singular value, leave-one-out diagnostics for weak
overlaps, and global frequency-path fit diagnostics.  ``--quiet`` suppresses
normal successful output.  The two options are mutually exclusive.

Adding Kieffer acoustic branches
--------------------------------

``add-kieffer`` is available only for primitive Gamma-only phonons calculated
with the identity supercell. It reads one elastic output at the same volume,
selected with ``--interface crystal``, and writes a new input whose calculated
Gamma frequencies are unchanged. The three Kieffer branches are stored as a
separate additive component.

Use ``--pressure-source auto`` for the usual case. A raw tensor is corrected
using pressure from its unstrained stress; a tensor produced with CRYSTAL's
``PRESSURE`` keyword is recognized as already corrected. Use
``--pressure-source manual --pressure VALUE`` only when no usable output-stress
pressure is available.

The ``energy-eos`` and ``energy-polynomial`` pressure sources require a
multi-volume energy path and are therefore available only through
``quantas qha add-kieffer``. A single HA state must use output-stress, applied
pre-stress, or an explicitly supplied pressure.

The stored block is deliberately opt-in.  Run the enriched file with
``--kieffer`` to add the acoustic contribution.  Running the same file without
the flag performs the ordinary Gamma-phonon calculation, which provides a
controlled comparison without editing either input.  The report identifies
activation and includes the Kieffer reference; the native HDF5 archive retains
the acoustic component separately from the totals.

.. note::

   HA itself does not require phonon branch continuity across volume.  The
   shared generator nevertheless records continuity when it can, so the same
   YAML can later be used safely by QHA.  See
   :doc:`../workflows/phonon_input_generation` for the scientific procedure.

Important distinctions
----------------------

* ``--temperature`` controls the calculated grid and stored arrays.
* unit options on ``run`` describe the YAML values; unit options on ``plot`` or
  ``export`` convert an existing result for presentation.
* plotting during ``run`` is a convenience.  The standalone ``plot`` command is
  preferable when several figure variants are required.
* ``plot --axis temperature`` produces one curve per selected sampled volume;
  ``plot --axis volume`` produces one curve per selected stored temperature.
  ``--volume`` and ``--temperature`` accept only exact native-grid values.
* ``plot --2d`` adds a V--T map when both coordinate axes contain at least two
  points.  It does not interpolate the result grid.
* q-point weights are taken from the input and normalized; the CLI does not
  derive missing symmetry multiplicities.

See :doc:`../workflows/phonon_input_generation` for input-generation science,
:doc:`../workflows/ha` for HA implementation choices,
:doc:`../tutorials/ha` for a complete MgO calculation, and
:doc:`../formats/phonon_yaml` for the input contract.

Generated command reference
---------------------------

.. click:: quantas.cli.ha:ha
   :prog: quantas ha
   :nested: full
