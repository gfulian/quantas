Validation matrix
=================

The release-candidate matrix is assembled capability by capability.  The rows
below record validations that are already frozen in the repository; remaining
major workflows will be added during the final validation pass.

.. list-table:: Current scientific validation matrix
   :header-rows: 1
   :widths: 18 22 19 19 12 10

   * - Capability
     - Reference
     - Observable
     - Comparison
     - Test / record
     - Status
   * - Energy EOS E--V
     - Seven-volume CRYSTAL/PBE MgO series
     - ``E0``, ``V0``, ``K0``, ``KP``, ``KPP`` and energy RMSE
     - Public BM3 fit against the frozen real-data regression values
     - :doc:`eos` and ``tests/examples/test_curated_examples.py``
     - validated
   * - Energy EOS analytical consistency
     - Synthetic BM3, modified-Tait, and SJEOS datasets
     - ``E(V)``, ``P(V)``, bulk modulus and pressure derivatives
     - Integrated energy forms against ``P(V) = -dE/dV`` and round-trip tests
     - :doc:`eos`
     - validated

.. admonition:: Work in progress

   The final release-candidate matrix will also identify, for every remaining
   major scientific capability, the reference dataset, comparison target,
   observable, tolerance, implementing test, and validation status.  The
   individual validation pages remain the authoritative detailed records until
   that consolidation is complete.
