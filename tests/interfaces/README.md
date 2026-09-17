# Quantum-mechanical interfaces

Interface tests are grouped here. Reference outputs from CRYSTAL, VASP and
Phonopy will be added as validated datasets become available.


## VASP b13 characterization

The generic VASP run parser is characterized first against user-generated MgO
outputs from VASP 5.4.4.  Repository fixtures are reduced directly from the
real files: numerical values are retained rather than synthesized, and POTCAR
data are not stored.  Tests cover directory/source resolution, atom ordering,
float64 fractional coordinates, old and flat ionic-state XML layouts,
VASP-5 energy-tag repair, OUTCAR energy cross-checks, forces, stress, and
explicit convergence markers.

The Energy EOS adapter is characterized separately from the generic parser.
It selects ``e_0_energy`` for the ground-state E--V contract, checks a
conservative electronic-energy signature (including k-point sampling and
pseudopotential labels), rejects optimization histories, and exercises the
shared primitive-cell normalization.  spglib-backed tests also cover an
already primitive structure and a conventional rocksalt MgO cell.

The compact fixtures do not establish scientific validation for every VASP 6
minor version.  Additional real-version fixtures should be added before making
such a claim.

- VASP elasticity: OUTCAR clamped/ionic/total tensor decomposition, explicit
  shear-order mapping, directory-source handling, reference pressure/stress,
  and no implicit CRYSTAL prestress conversion.

- VASP Gamma phonons: primitive-cell q=0 frequencies/eigenvectors, rigid translations, and explicit rejection of unsupported dispersion/supercell folding.
