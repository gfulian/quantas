
=====================
Quantas documentation
=====================

**Quantas** stands for **Quant**\itative **A**\nalysis of **S**\olids. It is an open source 
package in Python_ for the analysis of the thermodynamics, and elastic properties of solid 
phases starting from theoretical or experimental results.

.. _Python: https://www.python.org/

Features
========

- Calculation of thermodynamics of solid systems at harmonic approximation 
  (HA) level

- Calculation of both thermodynamics and thermoelastic properties of solids at selected
  pressure and temperature conditions via quasi-harmonic approximation (QHA):

- Calculation of the equation of state (EoS) from experimental data

- Analysis of the second-order elastic moduli

- Being written in Python 3, Quantas is completely **cross-platform**!

References
==========

If you use Quantas to produce data for a publication, you are kindly requested to cite the 
following work::

  Gianfranco Ulian and Giovanni Valdre'
  'QUANTAS, a Python software for the analysis of solids from ab initio quantum mechanical simulations and experimental data'
  Journal of Applied Crystallography 55, (pages) (2022)
  http://dx.doi.org/10.1107/S1600576722000085
  
Also, the theory behind the different kind of available calculations is discussed in specific
literature, and we kindly ask you to cite them accordingly.

License
=======

New BSD.


Contact
=======

Author: `Gianfranco Ulian <mailto:gianfranco.ulian2@unibo.it>`_

Acknowledgements
================

The development of QUANTAS was supported by the Regione Emilia Romagna project PA2019-11452/RER to Giovanni Valdr |egrave|. 
The authors wish to thank also the beta testers of QUANTAS for their feedback on the program.

.. |egrave| unicode:: U+00E8 .. grave accent on e
   :ltrim:


.. toctree::
   :maxdepth: 1
   :caption: Installation
   :hidden:

   installation/prerequisites
   installation/installation
   installation/releasenotes

.. toctree::
   :maxdepth: 1
   :caption: Theoretical background
   :hidden:

   background/background_qha
   background/background_eos
   background/background_soec

.. toctree::
   :maxdepth: 1
   :caption: Getting started
   :hidden:
   
   usage/commands
   usage/ha_options
   usage/qha_options
   usage/eosfit_options
   usage/soec_options
   
.. toctree::
   :maxdepth: 1
   :caption: Input files
   :hidden:
   
   inputs/input_generator
   inputs/qha_input
   inputs/eos_input
   inputs/soec_input

.. toctree::
   :maxdepth: 1
   :caption: Tutorials
   :hidden: 
   
   tutorials/qha_tutorial
   tutorials/eos_tutorial
   tutorials/soec_tutorial

.. toctree::
   :maxdepth: 1
   :caption: Quantas package
   :hidden:  
   
   modules/modules
