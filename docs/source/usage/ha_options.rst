.. _ha_options:

===================================
Harmonic Approximation calculations
===================================

  :Last updated: |today|
  :Author: **Gianfranco Ulian**

Harmonic approximation calculations have some minor options, mainly related to the choice
of the measurement units involved and temperature settings. These options can be seen from
prompting :code:`quantas ha --help` or :code:`quantas ha -h`.

.. code-block:: console

  Usage: quantas ha [OPTIONS] FILENAME
  
    Harmonic Approximation calculation.
  
    This command requires a file (FILENAME) that will be read to provide the
    input data for the calculations.
  
  Options:
    -o, --outfile out_file          Output file where data will be stored,
                                    without extension.
    -T, --temperature min max step  Temperature range provided as a tuple.
                                    [default: ((298.15, 298.15, 1.))]
    --eunit [Ha|eV|Ry]              Measurement unit for energy values.
                                    [default: (Ha)]
    --vunit [A|bohr]                Measurement unit for volume values.
                                    [default: (A)]
    --funit [cm^-1|THz|Hz]          Measurement unit for phonon frequency
                                    values.  [default: (cm^-1)]
    --tunit [K|C]                   Measurement unit for temperature values.
                                    [default: (K)]
    -q, --quiet                     Output will not be printed on screen.
    -d, --debug                     Activate debug option.
    -h, --help                      Show this message and exit.


Define a temperature range
==========================

:code:`-T min max step`, :code:`--temperature min max step`
-----------------------------------------------------------

Sets the temperature range over which harmonic thermodynamic properties are calculated. This 
range can be set with three numbers after the flag:

  - the first value represents the minimum temperature;
  - the second value represents the maximum temperature;
  - the third value represents the temperature step.

For example, a calculation between 400 K and 1200 K, with an increment of 0.1 K can be 
requested by promting:

  .. code-block:: console
  
    quantas ha -T 400 1200 0.1
    
.. note::

  Temperature values are float numbers.


Measurement units
=================

:code:`--eunit EUNIT`
---------------------

Sets the units for energy values (default :code:`Ha`). Possible choices are:
  
============ ================
EUNIT value  Measurement unit
============ ================
:code:`Ha`   Hartree
:code:`eV`   electronVolt
:code:`Ry`   Rydberg
============ ================


:code:`--vunit VUNIT`
---------------------

Sets the units for unit-cell volume values (default :code:`A^3`). Possible choices are:
  
============== ================
VUNIT value    Measurement unit
============== ================
:code:`A^3`    cubic Angstrom
:code:`bohr^3` cubic bohr
============== ================


:code:`--funit FUNIT`
---------------------

Sets the units for (phonon) frequency values (default :code:`cm^-1`). Possible choices are:
  
============== ================
FUNIT value    Measurement unit
============== ================
:code:`cm^-1`  wavenumber
:code:`THz`    TeraHertz
:code:`Hz`     Hertz
============== ================


:code:`--tunit TUNIT`
---------------------

Sets the units for temperature values (default :code:`K`). Possible choices are:
  
============== ================
FUNIT value    Measurement unit
============== ================
:code:`K`      Kelvin
:code:`C`      Celsius degrees
============== ================


