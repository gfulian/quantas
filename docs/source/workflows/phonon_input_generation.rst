Phonon input generation and mode continuity
===========================================

Purpose and scientific scope
----------------------------

The HA and QHA input generator converts external phonon calculations into the
normalized YAML contract consumed by Quantas.  This operation is more than a
change of file syntax.  A usable thermodynamic input must place static energies,
volumes, phonon frequencies, q-point integration weights, structural metadata,
and, when required, phonon branch identities on one consistent physical basis.

For a sampled volume series :math:`V_i`, the normalized dataset contains the
quantities

.. math::

   \left\{V_i,\;U_0(V_i),\;\nu_{qj}(V_i),\;w_q\right\},

with optional primitive-cell structures along the same volume path.  Harmonic
thermodynamics needs the spectra and integration weights to be internally
consistent.  Frequency-based QHA imposes an additional requirement: the mode labelled
:math:`j` at one volume must correspond to the same physical vibrational mode
at the neighbouring volumes.

.. important::

   A YAML file can be syntactically valid while being scientifically unsuitable
   for mode-resolved QHA.  Rectangular frequency arrays do not by themselves
   establish phonon-mode continuity.  Quantas therefore treats input generation,
   structural normalization, and mode-continuity assessment as part of the
   scientific workflow rather than as formatting conveniences.

Supported input routes
----------------------

The public ``ha inpgen`` and ``qha inpgen`` commands, and the corresponding
Python APIs, share one frontend-neutral generator.  The currently supported
interfaces are:

``crystal``
   One CRYSTAL phonon output or a list of independent single-volume CRYSTAL
   phonon outputs.  For a multi-volume list, Quantas verifies q-mesh
   compatibility, parses printed eigenvectors when available, and performs its
   own backend-neutral mode tracking.

``crystal-qha``
   One monolithic CRYSTAL QHA output.  Quantas reads the volume-dependent
   frequencies and preserves the mode-continuity result reported by the native
   CRYSTAL workflow.  The continuity provenance is therefore source-managed,
   not recomputed by the Quantas multi-file tracker.

``phonopy``
   A supported Phonopy input route normalized to the common phonon contract.
   If the selected reader does not expose the eigenvectors required by the
   Quantas tracker, multi-volume continuity is recorded as ``unknown`` rather
   than guessed.

The same generated YAML can be read by HA and QHA.  HA does not require
mode-by-mode continuity because it evaluates each sampled volume independently.
QHA with ``scheme=freq`` does require defensible branch correspondence.

Overall generation pipeline
---------------------------

For independent volume calculations, the scientific pipeline is

.. code-block:: text

   external phonon outputs
          |
          v
   code-specific parsers
          |
          +-- static energy, volume, q mesh, weights
          +-- primitive/supercell structural information
          +-- phonon frequencies
          +-- phonon eigenvectors when available
          |
          v
   structural and unit normalization
          |
          v
   multi-volume compatibility checks
          |
          v
   adjacent-volume phonon-mode tracking
          |
          +-- non-degenerate overlap assignment
          +-- degenerate-subspace matching
          +-- ambiguity diagnostics
          +-- leave-one-out validation of weak overlaps
          |
          v
   normalized HA/QHA YAML

The parser layer reports information extracted from the external code.  The
mode tracker is backend-neutral and operates only on normalized arrays.  This
separation is deliberate: code-specific parsing does not decide QHA policy.

Structural normalization
------------------------

Why a primitive normalization is needed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Phonon calculations are often performed in a supercell even when HA/QHA
thermodynamics are normalized to a smaller primitive cell.  Energies, volumes,
atom counts, mode counts, and structural tensors must all refer to a consistent
cell before the data can be combined safely.

For CRYSTAL, let the primitive direct-lattice matrix be :math:`A`, the
supercell lattice be :math:`B`, and the integer expansion matrix be :math:`D`.
With direct vectors stored by rows, CRYSTAL uses

.. math::

   B = D A.

Quantas therefore reconstructs the primitive lattice as

.. math::

   A = D^{-1} B.

The number of primitive repetitions represented by the supercell is

.. math::

   N_{\mathrm{rep}} = |\det D|,

and the determinant must be a non-zero integer within numerical tolerance.

Folding atomic positions
~~~~~~~~~~~~~~~~~~~~~~~~

If :math:`\mathbf f_s` is a fractional coordinate in the supercell, the
corresponding primitive coordinate is obtained from

.. math::

   \mathbf f_p
   = \left(\mathbf f_s D\right) \bmod 1.

Folded atoms are associated with atoms in the reference primitive structure
only among candidates of the same chemical species.  Quantas uses
minimum-image Cartesian distances to identify the closest reference atom.
Every primitive atom must receive exactly :math:`N_{\mathrm{rep}}`
translational copies.

The copies assigned to one primitive atom are averaged around the reference
position.  If :math:`\Delta\mathbf f_r` are their minimum-image displacement
vectors, the average displacement is

.. math::

   \overline{\Delta\mathbf f}
   = \frac{1}{N_{\mathrm{rep}}}
     \sum_r \Delta\mathbf f_r,

and the residual Cartesian displacement of replica :math:`r` is

.. math::

   \mathbf r_r
   = \left(\Delta\mathbf f_r
     - \overline{\Delta\mathbf f}\right) A.

Quantas records both

.. math::

   r_{\max}=\max_r |\mathbf r_r|

and

.. math::

   r_{\mathrm{RMS}}
   = \left(\frac{1}{N}\sum_r |\mathbf r_r|^2\right)^{1/2}.

The current structural reconstruction labels a path ``exact`` when the maximum
residual is at most :math:`10^{-6}` angstrom, accepts a translationally averaged
path up to :math:`5\times10^{-2}` angstrom, and rejects larger inconsistencies.

.. warning::

   Structural reconstruction is not a symmetry-recovery trick.  Quantas does
   not keep increasing a symmetry tolerance until the desired primitive atom
   count appears.  The reconstruction is constrained by the known supercell
   transformation, species identity, replica count, and translational
   residuals.  A source structure that is not compatible with this physical
   relation is rejected.

The resulting YAML may retain the source-cell provenance in memory while
writing a compact primitive ``structure`` path for downstream QHA and
thermoelastic use.

Units and normalization basis
-----------------------------

Current generated inputs carry explicit units, normally

.. code-block:: yaml

   units:
     energy: Ha
     volume: angstrom^3
     frequency: cm^-1
     length: angstrom

All volume-dependent quantities in one generated dataset must use one common
normalization cell.  The generator also records ``natom`` and ``formula_units``
so that per-cell and molar quantities can be related without inferring the
chemical normalization from file names.

Historical Quantas phonon YAML files without an explicit ``units`` mapping are
still interpreted using the established convention above.  This is a
backward-compatibility rule, not permission to omit units in newly generated
files.

Compatibility checks for a multi-volume series
----------------------------------------------

Before eigenvectors are compared or frequency arrays are stacked, all
single-volume sources must describe the same phonon sampling problem.  Quantas
requires equal:

- number of q-points;
- number of modes at every q-point;
- integer phonon supercell matrix;
- q-point coordinates and ordering, when coordinates are available;
- q-point weights;
- physical units.

For q-point coordinates and weights, the current comparison uses

.. math::

   |x_a-x_b| \le 10^{-12}

with zero relative tolerance.  This is intentionally a strict consistency
check: the mode tracker is not a tool for reconciling different q meshes.

.. warning::

   Do not combine calculations performed on different q meshes, different
   supercells, or different q-point orderings and expect mode tracking to repair
   the mismatch.  Such datasets represent different numerical problems.  The
   generator rejects them before continuity analysis.

CRYSTAL eigenvector representation
----------------------------------

Real and complex modes
~~~~~~~~~~~~~~~~~~~~~~

At Gamma, CRYSTAL prints real normal-mode displacement vectors.  At a general
q-point it prints separate in-phase and anti-phase components.  Quantas
reconstructs the complex displacement vector as

.. math::

   \mathbf u_{qj}
   = \mathbf u^{\mathrm{in}}_{qj}
   + i\,\mathbf u^{\mathrm{anti}}_{qj}.

The in-phase and anti-phase blocks must contain the same frequencies and atom
ordering.  Incomplete or inconsistent blocks are rejected.

Mass weighting and unit normalization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

CRYSTAL prints displacements normalized to classical amplitudes.  Let
:math:`u_{qj,\kappa\alpha}` be the printed displacement of atom :math:`\kappa`
in Cartesian direction :math:`\alpha`.  Quantas restores the mass weighting as

.. math::

   \widetilde e_{qj,\kappa\alpha}
   = \sqrt{m_\kappa}\,
     u_{qj,\kappa\alpha},

and then normalizes each mode:

.. math::

   e_{qj,\kappa\alpha}
   = \frac{\widetilde e_{qj,\kappa\alpha}}
   {\left(\sum_{\kappa\alpha}
   |\widetilde e_{qj,\kappa\alpha}|^2\right)^{1/2}}.

The resulting complex eigenvectors obey

.. math::

   \langle e_{qj}|e_{qj}\rangle = 1.

The classical amplitude and the source displacement length unit are scalar
factors for a given mode and therefore cancel in this final normalization.  The
stored representation is a unit-norm, mass-weighted direction suitable for
backend-neutral overlap comparisons.

Local mode tracking across volume
---------------------------------

Volume ordering
~~~~~~~~~~~~~~~

The supplied files need not be ordered by volume.  Quantas first constructs

.. math::

   V_{(1)} < V_{(2)} < \ldots < V_{(N)}

and performs matching only between adjacent states:

.. math::

   V_{(1)}\leftrightarrow V_{(2)},\quad
   V_{(2)}\leftrightarrow V_{(3)},\quad \ldots

This local construction minimizes the structural change between compared
states and prevents a distant reference volume from controlling all overlap
assignments.

.. important::

   The chosen ``reference_index`` does **not** affect the local overlap
   matrices or the local assignments.  It fixes only the final branch labels.
   Changing the reference is therefore a relabelling operation, not a different
   physical matching calculation.

Overlap matrix and global assignment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For one q-point and two adjacent volumes :math:`V_a` and :math:`V_b`, the
non-degenerate overlap matrix is

.. math::

   O_{ij}^{(a,b)}
   = \left|\left\langle
     e_i(V_a)\middle|e_j(V_b)
     \right\rangle\right|.

The absolute value makes the comparison invariant to an arbitrary complex
phase of either eigenvector.  For normalized vectors,
:math:`0\le O_{ij}\le1`.

Quantas does not independently select the largest value in every row.  A
one-to-one mode mapping is required, so it determines the permutation
:math:`\pi` that maximizes the total overlap

.. math::

   \pi^\star
   = \arg\max_\pi
     \sum_i O_{i,\pi(i)}.

The assignment is solved with the Hungarian algorithm.  This matters near mode
mixing: the best assignment for the complete set can legitimately differ from
the largest scalar product considered for one mode in isolation.

Ambiguous non-degenerate assignments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For a selected pair :math:`i\rightarrow\pi(i)`, Quantas also records the best
competitor in the same row or selected column,

.. math::

   C_i = \max\left[
   \max_{j\ne\pi(i)} O_{ij},
   \max_{k\ne i} O_{k,\pi(i)}
   \right],

and the separation

.. math::

   \Delta_i = O_{i,\pi(i)}-C_i.

The default ambiguity margin is

.. math::

   \Delta_i < 0.4.

This follows the diagnostic separation used by CRYSTAL for potentially unsafe
scalar-product competitors.  An ambiguous assignment is retained but marked as
``caution`` unless another criterion makes it unresolved.

.. caution::

   ``Ambiguous`` does not mean ``wrong``.  It means that the selected overlap is
   not cleanly separated from another plausible assignment.  Near a physical
   mode crossing or strong mixing, this is expected information and should be
   inspected rather than automatically converted into a failure.

Degenerate eigenspaces
----------------------

Why individual eigenvectors cannot be tracked
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Inside an exactly or numerically degenerate manifold, the eigenvectors are not
unique: any unitary rotation of the basis spans the same physical eigenspace.
Comparing one vector to one vector would therefore make the result depend on an
arbitrary diagonalization basis.

Quantas groups modes into a numerical degeneracy when their frequencies satisfy

.. math::

   |\nu_i-\nu_j|
   \le
   \epsilon_{\mathrm{abs}}
   + \epsilon_{\mathrm{rel}}
     \max(|\nu_i|,|\nu_j|,1),

with defaults

.. math::

   \epsilon_{\mathrm{abs}}=0.05\ \mathrm{cm}^{-1},
   \qquad
   \epsilon_{\mathrm{rel}}=10^{-6}.

Subspace overlap
~~~~~~~~~~~~~~~~

Let :math:`E_a` and :math:`E_b` contain orthonormal basis vectors spanning the
matched degenerate manifolds at adjacent volumes.  Quantas forms

.. math::

   S = E_a^\dagger E_b

and evaluates its singular values

.. math::

   \sigma_1,\ldots,\sigma_d.

The conservative measure is

.. math::

   \sigma_{\min}=\min_k \sigma_k.

The default acceptance condition is

.. math::

   \sigma_{\min}\ge0.8.

If the condition fails, all modes belonging to that matched manifold are
``unresolved``.

A useful limiting example is the three acoustic translations at Gamma.  Their
individual eigenvectors may rotate or exchange order between calculations,
but the translational three-dimensional subspace should remain unchanged.  A
subspace criterion tests this physically meaningful object directly.

.. warning::

   Do not interpret a low individual scalar product inside a degenerate
   manifold as evidence of broken continuity.  Individual vector identity is
   not invariant there.  The singular values of the matched subspaces are the
   relevant diagnostic.

Reference-labelled branch composition
-------------------------------------

After every adjacent pair has been matched independently, Quantas composes the
local permutations away from the reference state.  If
:math:`\pi_{12},\pi_{23},\ldots` are the local mappings, a branch label can be
propagated through the chain, schematically,

.. math::

   \Pi_{1\rightarrow4}
   = \pi_{34}\circ\pi_{23}\circ\pi_{12}.

On the compression side the corresponding inverse mappings are composed.  The
resulting ``permutations`` array maps the final branch labels onto the raw mode
indices printed by each source calculation.

The metadata distinguish two counts:

``local_reordered_assignments``
   Number of non-identity assignments in the adjacent-volume problems.

``reordered_assignments``
   Number of final branch-to-raw mappings that differ from identity after the
   local permutations have been composed relative to the reference labels.

The second number can be larger because a local exchange can remain propagated
through several later states.

Global frequency-path diagnostics
---------------------------------

Once branch identities have been constructed from eigenvectors, Quantas checks
how smoothly every tracked frequency path behaves with volume.  These fits are
**diagnostics**; they do not participate in the initial Hungarian assignment.

For :math:`N_V` sampled volumes the global diagnostic degree is

.. math::

   d = \min(3,N_V-2),

provided at least three volumes are available.  The fitted model is

.. math::

   \widehat\nu(V)=\sum_{k=0}^{d}a_k x(V)^k,

where :math:`x(V)` is a centered and scaled volume coordinate.  Quantas reports
for each q-point and branch:

.. math::

   R^2
   = 1-
   \frac{\sum_i[\nu_i-\widehat\nu(V_i)]^2}
        {\sum_i[\nu_i-\overline\nu]^2},

.. math::

   \mathrm{RMSE}
   = \left[
   \frac{1}{N_V}
   \sum_i[\nu_i-\widehat\nu(V_i)]^2
   \right]^{1/2},

and

.. math::

   r_{\max}
   = \max_i|\nu_i-\widehat\nu(V_i)|.

The current descriptive ``supported`` diagnostic uses
:math:`R^2\ge0.98` **or** :math:`\mathrm{RMSE}\le2\ \mathrm{cm}^{-1}`.
This flag summarizes branch smoothness but does not decide whether a weak
overlap is accepted.

.. important::

   A global polynomial fit is not an independent validation of one suspicious
   point because that point also contributes to the fitted coefficients.
   Quantas therefore never uses the global fit alone to rescue a low-overlap
   assignment.  Independent leave-one-out prediction is required instead.

Leave-one-out validation of low overlaps
----------------------------------------

A non-degenerate assignment is classified as low-overlap when

.. math::

   O_{i,\pi(i)} < 0.5.

Such an assignment is not accepted from the eigenvector evidence alone.  For
each endpoint of that adjacent-volume link, Quantas removes the endpoint under
test, fits the tracked branch to all remaining volumes, and predicts the
omitted frequency.

When enough volumes are present, the leave-one-out polynomial degree is

.. math::

   d_{\mathrm{LOO}}=\min(3,N_V-3).

For omitted state :math:`s`, the predictive residual is

.. math::

   r_s
   = \left|
   \nu_s-\widehat\nu_{(-s)}(V_s)
   \right|,

where :math:`\widehat\nu_{(-s)}` is fitted without state :math:`s`.  The
acceptance tolerance is

.. math::

   \tau_s
   = \max\left[
   2.0\ \mathrm{cm}^{-1},
   0.10\left(
   \nu_{\max}^{\mathrm{train}}
   -\nu_{\min}^{\mathrm{train}}
   \right)
   \right].

A low-overlap link remains usable only when both endpoints satisfy

.. math::

   r_a\le\tau_a
   \qquad\mathrm{and}\qquad
   r_b\le\tau_b.

Passing links are retained as ``caution``.  If either endpoint fails, the local
assignment is ``unresolved``.

The two-endpoint test is intentionally symmetric.  Continuity must not depend
on whether a series is mentally traversed toward compression or expansion.

.. warning::

   With only two sampled volumes there is no independent frequency path from
   which to predict either endpoint.  A weak overlap therefore cannot be
   rescued by a trivial two-point fit and remains unresolved.  More generally,
   sparse volume sampling limits what any continuity algorithm can establish.

Continuity status written to the YAML
-------------------------------------

The local diagnostics are summarized into the public QHA continuity contract.
For Quantas eigenvector tracking:

.. list-table:: Local mode-tracking interpretation
   :header-rows: 1
   :widths: 40 24 36

   * - Situation
     - Local status
     - Consequence
   * - Clear non-degenerate assignment
     - ``matched``
     - No continuity warning.
   * - Ambiguous but otherwise acceptable assignment
     - ``caution``
     - Retained and reported for inspection.
   * - Accepted degenerate eigenspace
     - ``degenerate``
     - Individual basis rotation is ignored.
   * - Low overlap passing both LOO predictions
     - ``caution``
     - Retained with predictive diagnostics.
   * - Low overlap failing either LOO prediction
     - ``unresolved``
     - Dataset cannot be verified for mode-resolved QHA.
   * - Degenerate subspace below the singular-value threshold
     - ``unresolved``
     - Dataset cannot be verified for mode-resolved QHA.

The dataset status is then

.. math::

   N_{\mathrm{unresolved}}=0
   \quad\Longrightarrow\quad
   \texttt{mode\_continuity: verified},

while any unresolved assignment gives

.. math::

   N_{\mathrm{unresolved}}>0
   \quad\Longrightarrow\quad
   \texttt{mode\_continuity: unreliable}.

Cautions are intentionally compatible with ``verified``: they report difficult
but still defensible assignments rather than erasing useful scientific
information.

The broader YAML contract also accepts ``assumed`` and ``unknown``.  Their
meaning is described in :doc:`../formats/phonon_yaml`.

Native CRYSTAL QHA continuity
-----------------------------

The ``crystal-qha`` route is different from the independent multi-output route.
If the CRYSTAL QHA output contains the native statement

.. code-block:: text

   FOUND CONTINUITY OF FREQUENCIES WITH VOLUME

Quantas records

.. code-block:: yaml

   mode_continuity: verified
   mode_continuity_metadata:
     method: crystal-qha
     source: crystal

This is a provenance statement: continuity was established by the source QHA
workflow.  Quantas does not relabel it as though its own adjacent-volume
tracker had performed the assessment.

.. note::

   The CRYSTAL eigenvector parser is also characterization-tested on the native
   QHA output, but the generated ``crystal-qha`` YAML deliberately records the
   source-managed method.  Provenance should state what actually established
   the scientific condition.

Reading the generated diagnostics
---------------------------------

The standard CLI output summarizes every adjacent volume interval with:

- number of local reorders;
- number of matched degenerate subspaces;
- low-overlap assignments;
- cautions and unresolved assignments;
- minimum non-degenerate overlap;
- minimum degenerate-subspace singular value.

``--debug`` adds one table for every q-point and adjacent-volume pair.  It shows
raw mode indices, frequencies, selected overlap, best competitor, overlap gap,
subspace singular value, local status, leave-one-out residual and tolerance when
applicable, and the global branch-fit diagnostics.

Debug frequencies are displayed to four decimal places, matching the maximum
precision printed in the characterized CRYSTAL phonon outputs.  This is a
renderer choice only; YAML values and in-memory ``float64`` data are not rounded
for calculation.

When terminal output is redirected to a file, Quantas uses deterministic plain
text rather than Rich table compression.  This prevents Unicode ellipses or
terminal-width truncation from altering numeric diagnostics.

Practical command-line examples
-------------------------------

One CRYSTAL phonon output
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: console

   quantas ha inpgen phonon.out \
       --interface crystal \
       --output material.yaml \
       --jobname "Material harmonic phonons"

A multi-volume CRYSTAL series
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create a plain text file containing one output path per line, then run

.. code-block:: console

   quantas qha inpgen files.txt \
       --list \
       --interface crystal \
       --reference 0 \
       --output material_qha.yaml \
       --jobname "Material QHA phonons"

Use ``--debug`` when the continuity diagnostics need to be inspected in detail.
Use ``--quiet`` for successful batch generation with no normal terminal output.
The two options are intentionally mutually exclusive.

Native CRYSTAL QHA output
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: console

   quantas qha inpgen qha.out \
       --interface crystal-qha \
       --output material_qha.yaml

Only one source file is accepted by the ``crystal-qha`` interface.

Scientific acceptance checklist
-------------------------------

Before using a generated multi-volume file for production QHA, verify all of
the following:

#. the structures represent the same phase and normalization cell;
#. the volume range brackets the thermodynamic states of interest;
#. the q mesh, q-point ordering, weights, and supercell are identical across
   independent volume calculations;
#. structural reconstruction residuals are physically negligible for the
   intended primitive representation;
#. ``mode_continuity`` is appropriate to the selected QHA scheme;
#. if continuity is ``verified`` with cautions, the cautions are localized and
   scientifically interpretable rather than widespread parser or convergence
   failures;
#. no unresolved low-overlap assignment or degenerate subspace is being hidden
   by a downstream fit;
#. the static energies and phonon frequencies were converged at a level
   appropriate to the requested thermodynamic derivatives.

.. warning::

   ``mode_continuity: verified`` is not a statement that the underlying phonon
   calculation is physically correct.  It establishes only that Quantas found a
   defensible correspondence between the supplied vibrational branches.  It
   does not test phase stability, basis-set or k-point convergence, force-
   constant convergence, anharmonicity, or the adequacy of the sampled volume
   range.

Related pages
-------------

- :doc:`ha` explains how the generated spectra are summed harmonically.
- :doc:`qha` explains how the multi-volume dataset is interpolated and minimized.
- :doc:`../formats/phonon_yaml` defines the normalized file contract.
- :doc:`../validation/ha_qha` documents the characterization and real-dataset
  validation used for this implementation.
- :doc:`../developer/interfaces` describes the parser/adapter architecture.
