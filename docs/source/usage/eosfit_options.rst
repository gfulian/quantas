.. _eos_options:

=====================================
Options for Equation of State fitting
=====================================

  :Last updated: |today|
  :Author: **Gianfranco Ulian**


The options for the interactive Equation of State fitting procedure are reported below. This
is the output of the :code:`quantas eos -h`.

.. code-block:: console

  Usage: quantas eosfit [OPTIONS] FILENAME
  
    Equation of state (EoS) fitting.
  
    This command requires a file (FILENAME) that will be read to provide the
    input data for the calculations.
  
  Options:
    -o, --outfile out_file  Output file for EoS log, without extension.
    --punit [GPa]           Measurement unit for pressure values.  [default:
                            (GPa)]
    --vunit [A]             Measurement unit for volume values.  [default: (A)]
    -h, --help              Show this message and exit.

