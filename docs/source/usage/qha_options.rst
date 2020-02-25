.. _qha_options:

=====================================================
Options for Quasi-Harmonic Approximation calculations
=====================================================

  :Last updated: |today|
  :Author: **Gianfranco Ulian**


The same options for the harmonic approximation applies, but there is also
the possibility to change the units used to describe pressure. In addition, 
there are several options to control how the quasi-harmonic properties
are calculated.

.. code-block:: console

  Usage: quantas qha [OPTIONS] FILENAME
  
    Quasi-Harmonic Approximation calculation.
  
    This command requires a file (FILENAME) that will be read to provide the
    input data for the calculations.
  
  Options:
    -o, --outfile out_file          Output file where data will be stored,
                                    without extension.
    -s, --scheme [freq|td]          QHA scheme, select between frequency (freq)
                                    or thermodynamic (td) interpolation.
                                    [default: (td)]
    -m, --minimization [poly|eos]   Volume minimization scheme, select between
                                    polynomial roots (poly) or equation of state
                                    fitting (eos).  [default: (poly)]
    -T, --temperature min max step  Temperature range provided as a tuple.
                                    [default: ((298.15, 298.15, 1.))]
    -P, --pressure min max step     Pressure range provided as a tuple.
                                    [default: ((0., 0., 1.))]
    --fdeg INTEGER RANGE            Set the degree of the polynomials used to
                                    fit phonon frequency values.  [default: 3]
    --edeg INTEGER RANGE            Set the degree of the polynomials used to
                                    fit energy values.  [default: 3]
    --eos [M|BM|PT|V]               Set the equation of state formulation.
                                    [default: (BM)]
    --eunit [Ha|eV|Ry]              Measurement unit for energy values.
                                    [default: (Ha)]
    --vunit [A|bohr]                Measurement unit for volume values.
                                    [default: (A)]
    --funit [cm^-1|THz|Hz]          Measurement unit for phonon frequency
                                    values.  [default: (cm^-1)]
    --tunit [K|C]                   Measurement unit for temperature values.
                                    [default: (K)]
    --punit [GPa|kbar]              Measurement unit for pressure values.
                                    [default: (GPa)]
    -q, --quiet                     Output will not be printed on screen.
    -d, --debug                     Activate debug option.
    -h, --help                      Show this message and exit.


Define a temperature range
==========================

:code:`-T min max step`, :code:`--temperature min max step`
-----------------------------------------------------------

Sets the temperature range over which quasi-harmonic thermodynamic properties are calculated. 
This range can be set with three numbers after the flag:

  - the first value represents the minimum temperature;
  - the second value represents the maximum temperature;
  - the third value represents the temperature step.

For example, a calculation between 400 K and 1200 K, with an increment of 0.1 K can be 
requested by promting:

  .. code-block:: console
  
    quantas ha -T 400 1200 0.1
    
.. note::

  Temperature values are float numbers.


Define a pressure range
=======================

:code:`-P min max step`, :code:`--pressure min max step`
-----------------------------------------------------------

Sets the pressure range over which quasi-harmonic thermodynamic properties are calculated. This
range can be set with three numbers after the flag:

  - the first value represents the minimum pressure;
  - the second value represents the maximum pressure;
  - the third value represents the pressure step.

For example, a calculation between 0 GPa and 10 GPa, with an increment of 2 GPa can be 
requested by promting:

  .. code-block:: console
  
    quantas ha -P 0 10 2
    
.. note::

  Pressure values are float numbers.


Thermodynamic calculation approach
==================================  

Two methods to calculate thermodynamic and derived properties can be employed in Quantas.
They can be selected as with the :code:`-s` or :code:`--scheme` option

:code:`-s td`, :code:`--scheme td`
----------------------------------

Interpolate (harmonic) thermodynamic properties at each considered temperature by polynomial
functions (default). Then, the same functions are employed to calculate these properties at 
specific *P-V-T* conditions.

:code:`-s freq`, :code:`--scheme freq`
--------------------------------------

Use frequency continuity to calculate thermodynamic properties at specific *P-V-T*
conditions. 

With this scheme, all frequency bands are interpolated by polynomial functions, which 
are then employed to calculate thermodynamics at unit-cell volumes at selected *P-T* 
settings.
  

Volume minimization scheme
==========================

Volume minimization at selected *P-T* conditions is performed by minimizing the Helmholtz 
free energy on isothermal :math:`F(V)` curves at target pressures. These curves can be
described by:

  - polynomial functions (*numerical* approach, default) 
  - equation of state (EoS) functions (*phenomenological* approach)
  
The minimization scheme can be selected appending the :code:`-m` or :code:`--minimization`
flag to the :code:`qha` sub-command.

:code:`-m poly`, :code:`--minimization poly`
--------------------------------------------

Use polynomial functions of selected degree to fit and minimize isothermal :math:`F(V)` curves.


:code:`-m eos`, :code:`--minimization eos`
------------------------------------------

Use equation of state functions to fit and minimize isothermal :math:`F(V)` curves.


Equation of state selection
===========================

In the case of volume minimization by EoS, its specific volume-integrated formulation can 
be selected by using the :code:`--eos` flag:

:code:`--eos M`
---------------

Use a volume-integrated Murnaghan equation of state of the form:

.. math::

  E = E_0 + K_0 \frac{V}{K^{\prime}} \Bigg[
  \frac{\big(V_0/V\big)^{K^{\prime}}}{{K^{\prime}} -1} + 1
  \Bigg] - V_0 \frac{K_0}{{K^{\prime}} - 1}

:code:`--eos BM`
----------------

Use a volume-integrated Birch-Murnaghan equation of state of the form (default):

.. math::

  E = E_0 + K_0 V_0 \frac{9}{16} \Big\{ K^{\prime} \Big(
  \eta^2 - 1 \Big)^{3} + \Big[ \big( \eta^{2} -1 \big)^{2}
  \big(6 - 4 \eta^{2} \big) \Big] \Big\}

.. math::

  \eta = \Bigg( \frac{V_0}{V} \Bigg)^{1 / 3}

:code:`--eos V`
---------------

Use a volume-integrated Vinet equation of state, expressed as:

.. math::

  E = E_0 + 2 \frac{ K_0 V_0}{\big(K^{\prime} - 1)}^2 \Big\{
  2 - \bigg[ 5 + 3K^{\prime} \big(\eta - 1\big) \bigg]
  e^{-3\big(K^{\prime}-1\big)\big(\eta -1\big) / 2}

.. math::

  \eta = \Bigg( \frac{V_0}{V} \Bigg)^{1 / 3}

:code:`--eos PT`
----------------

Use a volume-integrated Pourier-Tarantola (Natural Strain) equation of state, formulated
as:

.. math::

  E = E_0 + \frac{B_0 V_0 \rho^2}{6} \bigg(3 + \rho
  \big(K^{\prime}-2\big) \bigg)

.. math::

  \eta = \Bigg( \frac{V_0}{V} \Bigg)^{1 / 3}

.. math::

  \rho = -3 ln(\eta)

Polynomial fitting options
==========================

Polynomial fit options are related to the degree of the employed polynomial. Default values 
should provide adequate results, but the user may want to change them according to the
situation.

:code:`--edeg`
--------------

Set the degree of the polynomial function used to fit energy *vs* volume data. Here, the term
'energy' refers to any energy value calculated by Quantas (internal energy, entropy and so on).
These polynomial functions are employed if it is requested one or both of these options:

  - fit of thermodynamic properties (see :code:`--thermodynamics`)
  - minimization of volume by polynomials (see :code:`--polymin`)

The default value is 3.

:code:`--fdeg`
--------------

Set the degree of the polynomial function used to fit frequency *vs* volume data, and is
employed in conjuction with the frequency interpolation scheme (see :code:`--frequency`). 

The default value is 3.


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

:code:`--punit PUNIT`
---------------------

Sets the units for pressure values (default :code:`GPa`). Possible choices are:
  
============ ================
PUNIT value  Measurement unit
============ ================
:code:`GPa`  Gigapascal
:code:`kbar` kilobar
============ ================
