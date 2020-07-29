.. _input_generator:

=======================
Quantas input generator
=======================

  :Last updated: |today|
  :Author: **Gianfranco Ulian**

.. note::

  At the moment, only input for (quasi-)harmonic approximation and second-
  order elastic constants analysis can be generated from from the outputs of
  *ab initio* simulations.

Quantas is shipped with a script intended to aid the generation of input files from different
sources. This script is called from the console via:

.. code-block:: console

  quantas inpgen

To see the available options, type

.. code-block:: console

  quantas inpgen --help
  
or 

.. code-block:: console

  quantas inpgen -h

which outputs:

.. code-block:: console

  Usage: quantas inpgen [OPTIONS] COMMAND [ARGS]...
  
    Generate inputs file for Quantas from input files.
  
  Options:
    -h, --help  Show this message and exit.
  
  Commands:
    ha    Input generator for (Quasi-)Harmonic Approximation calculations.
    soec  Input generator for second-order elastic moduli analisys.

The generator *sub-command* should be chosen according to the data contained in the provided
input file. 

The mandatory argument is an input file that will be read by the software.

General options
===============

:code:`-o out_file`, :code:`--outfile out_file` option
------------------------------------------------------

By default, the generated input file for Quantas is named as ``quantas-`` + 
``generator_type`` + ``generator_extension``. By using this option, it is possible to
set a specific name for the generated input.


:code:`ha` sub-command
======================

By using :code:`quantas inpgen ha`, an input file for both harmonic and quasi-harmonic 
approximation calculations is generated, with ``.yaml`` extension.

Available *options* are listed below, which can be also obtained via prompting
:code:`quantas inpgen ha -h ` or :code:`quantas inpgen ha --help`:

.. code-block:: console

  Usage: quantas inpgen ha [OPTIONS] FILENAME
  
    Input generator for (Quasi-)Harmonic Approximation calculations.
  
    This command requires a file (FILENAME) that will be read to provide the
    input data for the input generation.
  
  Options:
    -o, --outfile out_file          Output file where data will be stored,
                                    without extension.  [default: (quantas_ha)]
    -l, --list                      Input files provided as a file list
    -r, --ref INTEGER               Reference file for (Q)HA input
    -i, --interface [crystal|crystal-qha|phonopy]
                                    Interface for ab initio codes.  [default:
                                    (crystal)]
    -h, --help                      Show this message and exit.

:code:`-l`, :code:`--list` option
---------------------------------
  
Usually, more that a single file will be processed to create an input for (quasi-)harmonic
approximation calculations. In this case, the user can feed :code:`quantas inpgen` with a 
single text file (with any extension, or even without one) containing a list of the required 
file names. For example, if four outputs from *ab initio* codes have to be read, it is possible
to put their file name in a file :code:`list.txt` with the following format:

.. code-block:: console

  file_1
  file_2
  file_3
  file_4
  
and we can give it to :code:`quantas inpgen ha` as:

.. code-block:: console

  > quantas inpgen ha list.txt --list
  
Then, Quantas creates the proper list of *ab initio* simulation output files to be read and 
processed.

.. warning::

  The list must contain **one** file name per line, without any comment or blank line between 
  them.

:code:`-r`, :code:`--ref` option
--------------------------------

Set the index of the file that is considered the reference of all the files provides as input. 
The reference unit cell (and file) should be that corresponding to the equilibrium geometry.
By default, following Python conventions, the first file has index 0.

:code:`-i`, :code:`--interface` option
--------------------------------------

Available interfaces are: 


  - :code:`crystal`: read CRYSTAL14/17 output files.

  - :code:`crystal-qha`: read a CRYSTAL17 output related to QHA calculations.

  - :code:`phonopy`: read phonopy output files, with extension ``.yaml``. 

.. warning:: 

  Phonopy outputs do not contain information of the original unit cell energy. 
  Quantas needs these energy values for the (Q)HA analysis. 
  
At present, two methods can be employed to set the energy of each unit cell:

  - VASP is fully supported, by providing the :code:`vasprun.xml` file of a single-point energy
    of the unit cell. The :code:`vasprun.xml` file **must** be renamed with the same root name 
    of the :code:`phonopy_output.yaml`. For example, is the phonopy output is called
    :code:`file_1.yaml`, the VASP xml file must be called :code:`file_1.xml`;
    
  - for all other *ab initio* software, the energy of the unit cell can be inputed by hand 
    in the generated input. Be careful to place the correct energy value at the same index of 
    the corresponding unit cell volume!



:code:`soec` sub-command
=========================

If selected, Quantas creates an input for the analysis of the second-order elastic constants,
with ``.dat`` extension. 
The generated input contains the elastic moduli expressed in GPa and the crystal density 
expressed in :math:`kg m^{-3}` (if available).

.. code-block:: console
  
  Usage: quantas inpgen soec [OPTIONS] FILENAME
  
    Input generator for second-order elastic moduli analisys.
  
    This command requires a file (FILENAME) that will be read to provide the
    input data for the calculations.
  
  Options:
    -o, --outfile out_file          Output file where data will be stored,
                                    without extension.
    -i, --interface [crystal|vasp]  Interface for ab initio codes.  [default:
                                    (crystal)]
    -h, --help                      Show this message and exit.

:code:`-i`, :code:`--interface` option
--------------------------------------

Available interfaces are: 


  - :code:`crystal`: read CRYSTAL14/17 output files.

  - :code:`vasp`: read vasp output files. 

.. warning::

  At the moment, only the OUTCAR file can be read. 


