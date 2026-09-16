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

The compact fixtures do not establish scientific validation for every VASP 6
minor version.  Additional real-version fixtures should be added before making
such a claim.
