.. _qha_input:

=====================================================
Input for (Quasi-)Harmonic Approximation calculations
=====================================================

  :Last updated: |today|
  :Author: **Gianfranco Ulian**

The input for HA and QHA calculation used by Quantas is a YAML_ file (extension ``.yaml``), which must contains serveral data that will be employed during the software run.

.. _YAML: https://yaml.org/

Each data type in the input is preceded by a ``keyword`` (in lower cases), which describes the
kind of information that will be collected. The format employed is reported in the following::

    keyword1:    value(s)
    keyword2:    [ (array of) values(s) separated by commas ]
    ...

In the following table, the keywords employed by Quantas are reported, alongside a short 
description of their functionality.

===================== ======================================================
Mandatory Keywords    Description
===================== ======================================================
**natom**             Number of atoms in the unit cell
**supercell**         3x3 matrix related to the unit cell expansion
**qpoints**           Number of **q** points used to calculate phonon 
                      properties
**volume**            List of unit cell volumes over over which the phonon 
                      properties were calculated
**energy**            List of energy values related to each unit cell volume
**phonon**            List of phonon values
===================== ======================================================

===================== ======================================================
phonon sub-keywords   Description
===================== ======================================================
**q-position**        Vector related to the specific *q* point
**weight**            Weight of the phonon band
**band**              Phonon band(s)
===================== ======================================================

===================== ======================================================
band sub-keywords     Description
===================== ======================================================
*# n*                 Band number
*frequency*           Array of *natom* phonon frequency values (float)
===================== ======================================================

For the sake of an example, if you calculated the phonon dispersion relations for a FCC lattice
containing 2 atoms in the primitive cell, the input file for the (Q)HA analysis will be 
something like the following:

.. code::

    natom:   2      
    supercell:
    - [     -2,      2,      2 ]
    - [      2,     -2,      2 ]
    - [      2,      2,     -2 ]
    qpoints: 656
    volume: [  V1, V2, V3, ... Vn ]
    energy: [  E1, E2, E3, ... En ]
    phonon:
    - q-position: [ q1_coord_1, q1_coord_2, q1_coord_2 ]
      weight: 1    
      band:
      - # 1
        frequency: [ v1_1, v1_2, v1_3, ..., v1_n ]
      - # 2
        frequency: [ v2_1, v2_2, v2_3, ..., v2_n ]
      - # 3
        frequency: [ v3_1, v3_2, v3_3, ..., v3_n ]
      - # 4
        frequency: [ v4_1, v4_2, v4_3, ..., v4_n ]
      - # 5
        frequency: [ v5_1, v5_2, v5_3, ..., v5_n ]
      - # 6
        frequency: [ v6_1, v6_2, v6_3, ..., v6_n ]
    - q-position: [ q2_coord_1, q2_coord_2, q2_coord_2 ]
      weight: 1    
      band:
      - # 1
        frequency: [ v1_1, v1_2, v1_3, ..., v1_n ]
      - # 2
        frequency: [ v2_1, v2_2, v2_3, ..., v2_n ]
      - # 3
        frequency: [ v3_1, v3_2, v3_3, ..., v3_n ]
      - # 4
        frequency: [ v4_1, v4_2, v4_3, ..., v4_n ]
      - # 5
        frequency: [ v5_1, v5_2, v5_3, ..., v5_n ]
      - # 6
        frequency: [ v6_1, v6_2, v6_3, ..., v6_n ]
    ...

.. note::

  A command_ for the automatic generation of the input file used by **Quantas** for the (quasi-)
  harmonic approximation analysis of crystalline solids is provided to aid the user. However, 
  at the moment, it is compatible only with the output files of CRYSTAL14_ (and above) and 
  phonopy_.
  
.. _CRYSTAL14: http://www.crystal.unito.it/index.php
.. _phonopy: https://atztogo.github.io/phonopy/
.. _command: ./input_generator.html

Keywords description
====================

A detailed description of the data requested for each keyword is reported in the following.


natom
-----

This is the number of atoms in the unit cell that was employed to perform the calculation of 
the phonon properties, either :math:`\Gamma`-point only or phonon dispersion relations.

:Example 1:
  
  You have calculated the phonon properties of NaCl by considering its crystallographic cell.
  In this case, you should write:
  
  .. code::

    natom: 8

:Example 2:

  Same as for Example 1, but instead you considered the primitive cell. In this case you 
  should write:

  .. code::

    natom: 2

.. note::

    ``natom`` is an integer type. Each unit cell volume employed to calculate the phonon 
    properties must contain the same number of atoms.


supercell
---------

The keyword ``supercell`` tells Quantas the expansion matrix used to calculate phonon 
dispersion relations.

:Example 1:

  If you have performed calculations using a :math:`2 \times 2 \times 2` expansion matrix, you 
  should write:

  .. code::
  
    supercell:
    - [      2,      0,      0 ]
    - [      0,      2,      0 ]
    - [      0,      0,      2 ] 

:Example 2:

  If the phonon calculations involved only :math:`\Gamma`-point (*i.e.* for a large unit 
  cell), you should write:

  .. code::
  
    supercell:
    - [      1,      0,      0 ]
    - [      0,      1,      0 ]
    - [      0,      0,      1 ] 

.. note::

    ``supercell`` is an integer type.


qpoints
-------

The keyword ``qpoints`` represents the number of **q** points sampled in the reciprocal space 
and used to calculate phonon properties (*i.e.*, phonon dispersion relations).

.. note::

    ``qpoints`` is an integer value.
    
.. warning::

    The number of **q** points indicated by ``qpoints`` must be consistent with the sum of the
    weights of each phonon band (see below).


volume/energy
-------------

The keywords ``volume`` and ``energy`` are each one followed by an array of :math:`n` values.

.. note::

    The values of ``volume`` and ``energy`` are float type. It is possible to consider a 
    single-volume calculation :math:`(n = 1)`, but at least :math:`(n > 4)` points are 
    required to perform quasi-harmonic approximation calculations.


phonon
------

``phonon`` represents the block of data containing the phonon properties of the material. 
**phonon** has some sub-keywords.

========================= ===================================================
**phonon** *sub-keywords*    
*q-position*              Vector related to the specific **q** point
*weight*                  Weight of the phonon band
*band*                    Phonon band(s)

**band** *sub-keywords*
*# n*                     Band number
*frequency*               Array of phonon frequency values (float)
========================= ===================================================

A brief explanation of each sub-keyword is here presented:

  - **q-position**: it represents one building block containing the phonon band calculated
    at that **q** point. It is a :math:`1 \times 3` array containing the fractional 
    coordinates of the sampled *q** point.
    
    - **weight** (integer number): this sub-keyword is a child of **q-position**, representing
      the number of times that this phonon band is seen when the phonon properties have been 
      calculated. This value is strictly related to the symmetry of the system under 
      consideration.
      
    - **band**: this sub-keyword is another child of **q-position**, and begins the declaration
      of the phonon frequencies calculated at that **q-position**.
      
      - **# n**: child of **band**, it is simply a label of the phonon frequency.
      
      - **frequency**: :math:`(1 \times natoms)` array of float values.

.. warning:
  
  The sum of all the *weights* in the **phonon** block must be equal to the number of **q**-
  points indicated in the **qpoints** section.
