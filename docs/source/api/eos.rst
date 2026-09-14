EOS API
=======

:mod:`quantas.api.eos` exposes an archive-oriented equation-of-state lifecycle.
Unlike single-shot workflows, one dataset may generate many immutable fit
records with different domains, formulations, solvers, constraints, and data
selections.  The public API therefore separates passive requests, fitting,
batch execution, archive/session state, and post-fit analysis.

Typical direct fit
------------------

.. code-block:: python

   from quantas.api import eos

   dataset = eos.read_input("PV_quartz.dat")
   request = eos.FitRequest(
       model="BM3",
       domain=eos.FitDomain.PRESSURE_VOLUME,
       options=eos.FitOptions(solver_options=eos.OLSOptions()),
   )
   eos.validate_request(dataset, request)
   result = eos.fit(dataset, request)

Typical Energy EOS fit
----------------------

.. code-block:: python

   from quantas.api import eos

   dataset = eos.read_input("EV_mgo_pbe.dat")
   request = eos.FitRequest(
       model="BM3",
       domain=eos.FitDomain.ENERGY_VOLUME,
       target="energy",
       options=eos.FitOptions(solver_options=eos.OLSOptions()),
   )
   result = eos.fit(dataset, request)
   print(result.parameter_values["V0"], result.parameter_values["K0"])

The public E--V adapter reports ``E0`` in Ha, ``V0`` in angstrom cubed, ``K0``
in GPa, ``KP`` dimensionless, and ``KPP`` in GPa :math:`^{-1}`.

Typical persistent batch
------------------------

.. code-block:: python

   plan = eos.BatchPlan(
       jobs=(eos.BatchJob(request=request, accept=True),),
   )
   batch = eos.run_batch(
       dataset,
       plan,
       "quartz_eos.hdf5",
       overwrite=True,
   )

The archive can then be diagnosed, calculated, and plotted without refitting.

Reference sections
------------------

.. toctree::
   :maxdepth: 1

   eos/contracts
   eos/fitting
   eos/batch_archive
   eos/postfit

See also
--------

- :doc:`../workflows/eos`
- :doc:`../tutorials/eos`
- :doc:`../formats/eos_input`
- :doc:`../formats/eos_spec`
- :doc:`../formats/eos_hdf5`
