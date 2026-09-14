``quantas qha``
===============

Energy-EOS options accept the compact Quantas model tags and aliases.
Use ``quantas eos show-models --domain ev`` for the current integrated
model catalogue and capability summary.

The QHA frontend builds a volume-dependent free-energy representation and
minimizes it at every requested pressure-temperature state.  Its command line
therefore exposes scientific choices that are absent from HA: interpolation
scheme, minimization model, polynomial degrees, derivative strategy,
thermal-expansion route, and local-fit failure policy.

Recommended sequence
--------------------

.. code-block:: console

   quantas qha inpgen qha-outputs.txt --list --output material.yaml
   quantas qha add-kieffer material.yaml --elastic-list elastic-files.txt \
      --interface crystal --pressure-source energy-eos --eos BM3 \
      --output material-kieffer.yaml
   quantas qha inspect material-kieffer.yaml --eos BM3 --degree 3
   quantas qha run material-kieffer.yaml --kieffer --scheme freq \
      --no-mode-gruneisen --minimization poly --temperature 0 1000 25 \
      --pressure 0 10 1
   quantas qha plot material-kieffer_QHA.hdf5 --property VT --property alphaV --2d
   quantas qha plot material-kieffer_QHA.hdf5 --property VT --axis pressure \
      --temperature 300 --temperature 1000
   quantas qha export material-kieffer_QHA.hdf5 --property VT --format csv

Use ``inspect`` before a production run.  It compares the sampled static
energy-volume data with polynomial and EOS previews and reports the implied
pressure support.  A dense requested P--T grid does not extend the support of
the sampled volumes.

Generating and checking the phonon input
----------------------------------------

For independent CRYSTAL phonon calculations at several volumes, place one
output path per line in a text file and run

.. code-block:: console

   quantas qha inpgen files.txt --list --interface crystal \
      --reference 0 --output material.yaml

When printed eigenvectors are available, Quantas validates that all sources use
the same q mesh, q-point weights, supercell, units, and mode count.  It then
tracks modes between adjacent volumes, treats numerical degeneracies as
subspaces, and writes the resulting continuity status and diagnostics to the
YAML.

A monolithic CRYSTAL QHA output uses the source-managed route:

.. code-block:: console

   quantas qha inpgen qha.out --interface crystal-qha \
      --output material.yaml

If CRYSTAL reports that frequency continuity with volume was found, Quantas
records ``mode_continuity: verified`` with ``method: crystal-qha`` rather than
claiming that its own multi-file tracker established the result.

Use ``--debug`` to inspect ambiguous or low-overlap assignments.  Use
``--quiet`` for silent successful batch generation; ``--quiet`` and ``--debug``
are mutually exclusive.

.. warning::

   ``inpgen`` refuses incompatible q meshes and records unresolved mode
   assignments instead of silently reordering a scientifically unsupported
   dataset.  Do not edit the resulting continuity status by hand to bypass the
   ``freq``-scheme preflight check.

The equations and acceptance criteria are documented in
:doc:`../workflows/phonon_input_generation`.

Adding volume-resolved Kieffer branches
---------------------------------------

``add-kieffer`` requires one CRYSTAL ELASTCON or ELAPIEZO output for every QHA
volume. The files may be supplied as positional arguments or through
``--elastic-list``. Select the reader with ``--interface crystal``, following
the same interface naming used by ``inpgen``. Quantas sorts the independently
calculated elastic states, matches them to the QHA volumes under the explicit
matching policy, applies the CRYSTAL finite-pressure correction when necessary,
and writes a separate ``*-kieffer.yaml`` input.

For raw elastic tensors, pressure can be reconstructed directly from the
static ``volume`` and ``energy`` arrays already stored in the QHA input:

.. code-block:: console

   quantas qha add-kieffer material.yaml --elastic-list elastic-files.txt \
      --interface crystal --pressure-source energy-eos --eos BM3 \
      --output material-kieffer-eos.yaml

   quantas qha add-kieffer material.yaml --elastic-list elastic-files.txt \
      --interface crystal --pressure-source energy-polynomial --degree 3 \
      --output material-kieffer-poly.yaml

Both routes evaluate :math:`P(V)=-dE/dV` at the sampled phonon volumes. The
generated Kieffer provenance records the exact EOS or polynomial degree, fit
parameters and diagnostics, input units, evaluated pressures, warnings, and
every elastic-to-phonon volume association. At least three volume-energy points
are required; a selected polynomial also needs enough points for its degree,
while the energy EOS may impose a stricter model-specific minimum. Use
``--pressure-source manual`` when the dataset is insufficient for a fit.
The polynomial route centres and scales the sampled volume coordinate before
fitting; the transform is stored with the coefficients so that conditioning
and the physical :math:`dE/dV` derivative remain independently inspectable.

``energy-eos`` and ``energy-polynomial`` deliberately require raw elastic
tensors. Their parsed output-stress value is not substituted into the fit.
The derived pressures are attached first and the hydrostatic Wallace
CRYSTAL finite-pressure correction is then applied once, with both operations
retained in provenance.
The hydrostatic assumption is not valid for a path carrying substantial
deviatoric stress.

The directional integration defaults to a 12 by 24 coarse quadrature refined
by a factor of two. ``--mu-order``, ``--phi-order``, and
``--refinement-factor`` are convergence controls; changes should be supported
by a sensitivity test for the material under study.

``add-kieffer`` only stores the validated data.  ``qha run`` reads and applies
them when ``--kieffer`` is present; omission of the flag leaves the block
inactive and reproduces the ordinary phonon-only workflow.  Both ``freq`` and
``td`` schemes are supported, and the sampled acoustic component remains
separate in the HDF5 result even though the total properties include it.

For ``--scheme freq``, mode-Gruneisen analysis cannot yet include the continuous
acoustic branches.  Quantas therefore disables its normal CLI default when
``--kieffer`` is selected and records that decision in the run options.  An
explicit ``--mode-gruneisen`` or ``--thermal-expansion mode_gruneisen`` request
is rejected.  Production scripts should state ``--no-mode-gruneisen`` as in the
example above; ``mixed_derivative`` and ``numerical`` thermal expansion remain
available.

Choosing options
----------------

* ``--scheme=freq`` retains mode-resolved information but requires defensible
  mode continuity.  ``--scheme=td`` interpolates integrated harmonic
  properties and is less dependent on branch tracking.
* ``--minimization=poly`` is flexible near a well-sampled minimum;
  ``--minimization=eos`` imposes a selected physical EOS form.
* ``--thermal-expansion`` chooses among mixed-derivative, mode-Gruneisen, and
  numerical-volume routes.  Their agreement is a diagnostic, not an identity
  guaranteed for every dataset.
* ``--poly-grid-points`` and ``--poly-grid-separation`` control local
  derivatives after polynomial minimization.  They do not change the
  equilibrium volume itself.
* ``--failure-policy`` determines whether failed local states terminate,
  accumulate, or raise immediately; it does not convert an unsupported state
  into a valid one.
* ``plot --axis temperature`` selects exact native pressures with
  ``--pressure``; ``plot --axis pressure`` selects exact native temperatures
  with ``--temperature``.  Plot construction never interpolates or snaps the
  requested coordinate.

See :doc:`../workflows/phonon_input_generation` for input-generation science,
:doc:`../workflows/qha` for the decision guide,
:doc:`../tutorials/qha` for reproducible calculations and method comparisons,
and :doc:`../formats/phonon_yaml` for input details.

Generated command reference
---------------------------

.. click:: quantas.cli.qha:qha
   :prog: quantas qha
   :nested: full
