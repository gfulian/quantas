Thermoelastic scientific validation
===================================

Validation scope
----------------

The thermoelastic validation matrix separates four questions that must not be
collapsed into one numerical comparison:

#. Can the CRYSTAL interface recover the original ``PRESSURE`` value, final
   structure, energy, volume, density and Wallace stiffness tensor?
#. Does the cold finite-strain implementation reproduce an independently
   assembled linear least-squares problem?
#. Does the quasi-static reconstruction preserve historical Quantas results
   when the same tensor frame and inputs are used?
#. Do the adiabatic tensors satisfy the independent thermodynamic identity,
   tensor symmetry, positive-semidefinite rank-one correction and exact
   zero-temperature limit?

Two real reference systems are used:

``MgO``
   Cubic periclase with eight CRYSTAL SOEC calculations spanning applied
   pressures from -16 to 16 GPa and a QHA grid spanning 0--20 GPa and
   0--2000 K.

``dolomite``
   Low-trigonal CaMg(CO3)2 with eleven CRYSTAL SOEC calculations spanning
   -6 to 20 GPa and a QHA grid spanning 0--10 GPa and 0--1700 K.  This system
   exercises non-zero C14 and C15 terms, tensor-frame normalization and
   anisotropic thermal expansion.

The large electronic-structure and QHA files are not distributed in the test
suite.  Instead, a compact frozen fixture stores the exact source SHA-256
identifiers, volume-dependent observations, fitted targets and selected P--T
states.  This keeps routine tests small while preserving traceability to the
full phase-validation archive.

Validated invariants
--------------------

The real-data tests require all of the following:

* direct CRYSTAL parsing agrees with the generated YAML to the text precision
  of the output files;
* every source contains an explicit ``PRESSURE`` keyword and its value agrees
  with the pressure attached to the corrected elastic tensor;
* an independent NumPy least-squares construction recovers every fitted
  ``C0`` and ``Cprime`` parameter;
* the frozen normalized isothermal reference tensors are reproduced to
  numerical roundoff;
* the adiabatic correction agrees with a separate evaluation of
  ``(T V / C_V) outer(C^T alpha, C^T alpha)``;
* the correction is symmetric, positive semidefinite and of numerical rank no
  greater than one;
* ``C^S = C^T`` exactly at 0 K;
* cubic MgO receives equal corrections to C11 and C12 and no correction to
  C44;
* all states in the declared validation domains remain mechanically stable.

Reference results
-----------------

The frozen reference calculations give the following principal
checkpoints:

.. list-table::
   :header-rows: 1
   :widths: 24 20 20 18 18

   * - System
     - P--T states
     - Minimum Wallace eigenvalue
     - Maximum :math:`|C^S-C^T|`
     - Elastic extrapolations
   * - MgO
     - 4221
     - 102.277 GPa
     - 30.2562 GPa
     - 126
   * - Dolomite, normalized frame
     - 3591
     - 26.7303 GPa
     - 5.14159 GPa
     - 0

The MgO extrapolations occur only where high-temperature QHA expansion places
volumes above the largest static elastic sample.  They are retained and marked;
they are not silently clipped.

Input-frame normalization
-------------------------

Thermoelastic YAML schema 1.0 is the pre-release normalized schema. Every
elastic point must store its source lattice, right-polar-decomposition
co-rotation, removed rotation angle, principal logarithmic strain, and
ordered-atom correspondence diagnostic. Earlier unnormalized schemas are no
longer accepted. Users must regenerate them from the original CRYSTAL outputs
with ``quantas thermoelasticity inpgen``.

The dolomite validation demonstrates why this is required: removing rigid
rotations up to 1.7392 degrees changes some reconstructed stiffness entries by
as much as 1.1254 GPa. No automatic numerical migration is attempted because
the missing frame provenance cannot be reconstructed reliably from an old YAML
alone.

Interpretation
--------------

This matrix validates parsing, finite-strain fitting, tensor reconstruction,
uncertainty-preserving backward compatibility, frame migration, adiabatic
conversion and mechanical stability for two independent symmetry classes.  It
does not establish universal accuracy of the quasi-static approximation.  A
QSA calculation omits explicit vibrational derivatives with respect to strain,
and its physical accuracy must still be assessed against full strain-dependent
QHA, experiment or another finite-temperature method for the material of
interest.

Literature context
------------------

The MgO case is appropriate because its three independent single-crystal
elastic constants and their pressure dependence have been extensively studied
experimentally and by first-principles methods.  The dolomite case provides a
low-symmetry carbonate whose ambient single-crystal tensor has also been
measured by Brillouin spectroscopy.

* Karki, B. B., Stixrude, L., Clark, S. J., Warren, M. C., Ackland, G. J., and
  Crain, J. (1997). Structure and elasticity of MgO at high pressure.
  *American Mineralogist*, **82**, 51--60.
  https://doi.org/10.2138/am-1997-1-207
* Sinogeikin, S. V. and Bass, J. D. (1999). Single-crystal elasticity of MgO
  at high pressure. *Physical Review B*, **59**, R14141--R14144.
  https://doi.org/10.1103/PhysRevB.59.R14141
* Jiang, F., Speziale, S., and Duffy, T. S. (2006). Elasticity of magnesite and
  dolomite from a genetic algorithm for inverting Brillouin spectroscopy
  measurements. *Physics of the Earth and Planetary Interiors*, **155**,
  1--20. https://doi.org/10.1016/j.pepi.2005.08.004
