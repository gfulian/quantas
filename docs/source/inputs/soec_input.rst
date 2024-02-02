.. _soec_input:

=================================================
Input for Second-Order Elastic Constants Analysis
=================================================

  :Last updated: |today|
  :Author: **Gianfranco Ulian**

The input file for the analysis of the second-order elastic constants (SOEC)
can be of any extension, albeit ``.dat`` should be preferred.

.. note::

    No keyword is necessary for this input.

The file is a text file containing the 21 elastic constants in either upper 
triangular, lower triangular and complete form.

The *first line* is intended for a short description of the system under 
analysis, but it can be omitted if the user wants so. Then the six rows of the
second-order elastic constants has to be supplied in Voigt's notation.

After the elastic constants, it can be supplied the density of the crystal 
system, expressed as :math:`kg m^{-3})`. If not supplied, the seismic wave 
velocities will not be calculated as they require this information.

An example of the input file for hydroxylapatite (hexagonal crystal) is the 
following: [1]_

.. code::

   Hydroxylapatite 
   187.208   65.193   84.703    0.000    0.000    0.000 
            187.208   84.703    0.000    0.000    0.000 
                     222.658    0.000    0.000    0.000 
                               39.687    0.000    0.000 
                                        39.687    0.000 
                                                 61.007 
   3178


.. note::

  A command_ for the automatic generation of the input file used by 
  **Quantas** for the second-order elastic constants analysis of crystalline 
  solids is provided to aid the user. However, at the moment, it is compatible 
  only with the output files of CRYSTAL14_ (and above) and VASP_.


.. note::

  .. versionadded:: 0.9.1

    The same input for the analysis of the second-order elastic moduli can be 
    employed for a complete analysis of the propagation of acoustic waves in 
    homogeneous, crystalline solids.
  
.. _CRYSTAL14: http://www.crystal.unito.it/index.php
.. _VASP: https://www.vasp.at/
.. _command: ./input_generator.html

.. rubric:: References

.. [1] Ulian, G., Valdre, G. Second-order elastic constants of hexagonal 
   hydroxylapatite (P63) from ab initio quantum mechanics: comparison between 
   DFT functionals and basis sets, Int. J. Quantum Chem., 118 (2018) e25500.