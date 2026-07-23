``quantas eos``
===============

EOS differs intentionally from the other scientific groups.  One dataset is
often fitted repeatedly with different formulations, solvers, selections,
constraints, and initial values.  Quantas therefore stores immutable fit
records in a persistent archive and tracks which record is accepted or marked
as a candidate for each scientific slot.

Recommended sequence
--------------------

For one homogeneous fit:

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

Important option families
-------------------------

* domain and target options define the scientific slot;
* P--V, V--T, and P--V--T options define the model rather than solver behavior;
* MGD normalization is part of the physical model and must match the volume
  basis;
* ``--fix``, ``--initial``, and ``--bound`` have different scientific meanings;
* OLS, WLS, effective variance, and ODR use different uncertainty assumptions;
* covariance scaling affects parameter uncertainties, not the optimum itself;
* ``--dry-run`` validates the resolved plan without fitting or creating HDF5;
* ``--failure-policy=continue`` continues later independent jobs but does not
  make the batch successful.

See :doc:`../workflows/eos` for implementation and record semantics,
:doc:`../tutorials/eos` for the P--V/V--T/P--V--T course, and
:doc:`../formats/eos_spec` plus :doc:`../formats/eos_hdf5` for persistence.

Generated command reference
---------------------------

.. click:: quantas.cli.eos:eos
   :prog: quantas eos
   :nested: full
