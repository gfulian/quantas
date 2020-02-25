.. _commands:

=====================================
Command Line Interface (CLI)
=====================================

  :Last updated: |today|
  :Author: **Gianfranco Ulian**

Quantas is a Command Line Interface (CLI) tool. It was intended in this way mainly for three reasons:

  1. with respect to Graphical User Interfaces (GUIs), it saves computing 
     resources for... computing!
  
  2. since it does not need any graphical environment, it is possible to 
     install the software on High-Performance Computing systems and use it on-site.
     
  3. development time.

After the installation of the software, you can see the different command-line arguments and options by simply writing on your console (UNIX), command-propt (Windows) or terminal (Mac OS).

.. code-block:: console

  quantas -h
  
or 

.. code-block:: console

  quantas --help
  
The output should be this one:

.. code-block:: console

  Usage: quantas [OPTIONS] COMMAND [ARGS]...
  
    ________                       __
    \_____  \  __ _______    _____/  |______    ______
     /  / \  \|  |  \__  \  /    \   __\__  \  /  ___/
    /   \_/.  \  |  // __ \|   |  \  |  / __ \_\___ \
    \_____\ \_/____/(____  /___|  /__| (____  /____  >
           \__>          \/     \/          \/     \/
  
  
  Options:
    -v, --version  Show the software version and exit.
    -h, --help     Show this message and exit.
  
  Commands:
    eosfit  Equation of state (EoS) fitting.
    export  Export results from binary format (HDF5) to text format.
    ha      Harmonic Approximation calculation.
    inpgen  Generate inputs file for Quantas from input files.
    qha     Quasi-Harmonic Approximation calculation.
    soec    Second-order elastic moduli analisys.


General concepts
================

Quantas routines are subdivided in different ``sub-commands``, which can be individually
called. Each command is followed by some parameters.

Parameters
----------

The parameters used by Quantas are of two types:

  - *arguments*: they are positional parameters, for example ``input_file_name.yaml`` when
    ``quantas ha input_file_name.yaml`` launching harmonic approximation 
    calculations;
  - *options*: they are called by using flags (for example ``-f`` or ``--flag-name``), 
    followed by a single value.
    

Multi-value options
-------------------

Some *options* of ``quantas`` are associated with multiple values, for example the temperature
settings employed during (quasi-)harmonic approximation calculations. In this case, the appropriate number of (expected) values are to be provided.


Case-sensitive options
----------------------

Some *options* of ``quantas`` are case-sensitive, in particular those related to the 
measurement unit selection.


Getting help
------------

Each ``quantas`` sub-command has an associated help string that provides guidance to the user 
on the mandatory arguments and options that can be provided. The help strings are called by 
appending ``-h`` or ``--help`` after the selected command, for example 
``quantas soec --help``, or just ``quantas soec -h``.
  

Printing/output options
-----------------------

The following are general options that can be employed in (most) sub-commands.

:code:`-o OUTFILE`, :code:`--outfile OUTFILE`
---------------------------------------------

Specify the output file name. By default, Quantas creates an output file composed as
``input_file_basename`` + ``_COMMAND_NAME.log``.

:code:`--q`, :code:`--quiet`
----------------------------

Suppress streaming information on the console, reporting the output only the log file.
  
:code:`-p`, :code:`--plot`
--------------------------

Activate plotting options for calculators.

.. note::

  At the moment, only the SOEC calculator has the capability of doing plot of the results.
  
:code:`--dpi DPI`
-----------------

Set the resolution (dot-per-inch, DPI) of the output plot figures (default: 80).

:code:`-d`, :code:`--debug`
---------------------------

Activates debugging information on screen. Useful to report Quantas bugs or strange behaviour 
of the code on some systems.


Aborting a calculation
----------------------

If the user needs to stop a currently running calculation, just press :kbd:`Ctrl+C` 
and Quantas will gently stop execution.

