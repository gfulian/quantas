.. _soec_options:

=====================================================
Options for Second-order elastic constants analysis
=====================================================

  :Last updated: |today|
  :Author: **Gianfranco Ulian**


The options for the analysis of the second-order elastic constants (or moduli) tensor
are briefly presented below. This is the output of the :code:`quantas soec -h`.

.. code-block:: console

  Usage: quantas soec [OPTIONS] FILENAME
  
    Second-order elastic moduli analisys.
  
    This command requires a file (FILENAME) that will be read to provide the
    input data for the calculations.
  
  Options:
    -o, --outfile out_file  Output file where data will be stored, without
                            extension.
    --punit [GPa]           Measurement unit for pressure values.  [default:
                            (GPa)]
    -p, --plot              Create plots of the results.
    --dpi INTEGER           Resolution (DPI) of the output figures.  [default:
                            80]
    -q, --quiet             Output will not be printed on screen.
    -d, --debug             Activate debug option.
    -h, --help              Show this message and exit.

