.. _qha_tutorial:

=======================================
(Quasi-)Harmonic Approximation tutorial
=======================================

  :Last updated: |today|
  :Author: **Gianfranco Ulian**

Preliminary operations
======================

Download the :download:`periclase (MgO) input file <../downloads/mgo_b3lyp_qha.yaml>` and put
it in a folder of your choice. From a console (or command prompt in Windows), go in that 
folder.

This input file was generated from several CRYSTAL17_ simulations performed on the primitive 
cell of MgO:

  - geometry optimization of eleven unit-cell volumes of the mineral, between 82% and 112% the 
    equilibrium cell volume;
  - calculation of the phonon dispersion relations at each unit-cell volume, using a very 
    simple supercell expansion of this form
    
    .. math::
    
      \begin{bmatrix}
      -2 &  2 &  2 \\
       2 & -2 &  2 \\
       2 &  2 & -2
      \end{bmatrix}
    
    that led to the sampling of 32 :math:`k`-points.

.. note::

  CRYSTAL17 units are Hartree for energy, :math:`\require{mediawiki-texvc}\AA^3` for volume and :math:`cm^{-1}` for 
  phonon frequency.
  

.. _CRYSTAL17: http://www.crystal.unito.it/index.php


Harmonic approximation
======================

Let's perform a very simple calculation of the harmonic (constant-volume) thermodynamic properties of this mineral phase. We will set a temperature range from 0 to 2000 K, with an 
increment of 10 K.

To run this calculation, simply write:

.. code-block:: console

  > quantas ha mgo_b3lyp_qha.yaml -T 0 2000 10
  
.. note::

  We could also include the measurement units for the data reported in the input file and 
  for the temperature range with:
  
  .. code-block:: console

    > quantas --ha mgo_b3lyp_qha.yaml -T 0 2000 10 --eunit Ha --vunit A^3 --funit cm^-1 --tunit 
    K 
  
  However, since these are the default units employed by Quantas, we skipped this declaration
  for simplicity.
  
The software prints out on the console the setting we chose:

.. code-block:: console

  ________                       __
  \_____  \  __ _______    _____/  |______    ______
   /  / \  \|  |  \__  \  /    \   __\__  \  /  ___/
  /   \_/.  \  |  // __ \|   |  \  |  / __ \_\___ \
  \_____\ \_/____/(____  /___|  /__| (____  /____  >
         \__>          \/     \/          \/     \/ 
                                              v0.9.0
  Authors: Gianfranco Ulian and Giovanni Valdre'
  Copyright (c) Gianfranco Ulian and Giovanni Valdre'.
  
  
  Calculator: Harmonic Approximation
  
  Temperature settings
  -------------------------------------
   - min:      0.00 K
   - max:   2000.00 K
   - step:    10.00 K
  
  Measurement units
  -------------------------------------
   - energy:      Ha
   - lenght:      A
   - frequency:   cm^-1
   - temperature: K
   
After this setup, Quantas reads the input file and tells us what it does contain:

.. code-block:: console
  
  Job: Quasi-Harmonic analysis of periclase (MgO)
  
  System:
  - Number of volumes                        11
  - Number of atoms                          2
  - Number of sampled k-points               32
  - Number of frequencies                    192

  Volume and energy values:
  
    Volume (A^3)            Energy (Ha)
   --------------- ------------------------------
      15.495427         -2.754467055837e+02
      16.068019         -2.754523424409e+02
      16.654547         -2.754566420093e+02
      17.255176         -2.754597306884e+02
      17.870076         -2.754617227028e+02
      18.499413         -2.754627209751e+02
      18.896862         -2.754628820004e+02
      19.143355         -2.754628184464e+02
      19.802069         -2.754620989310e+02
      20.475723         -2.754606377438e+02
      21.164485         -2.754585023425e+02

Then, the harmonic approximation calculation starts:

.. code-block:: console

  #------------------Harmonic Approximation calculation started------------------#
  
   - Start calculation of zero-point energy
     Finished, elapsed time     0.00 sec
  
   - Start calculation of thermal internal energy
     Finished, elapsed time     0.01 sec
  
   - Start calculation of entropy
     Finished, elapsed time     0.01 sec
  
   - Start calculation of isochoric heat capacity
     Finished, elapsed time     0.01 sec
  
   - Start calculation of vibrational Helmholtz free energy
     Finished, elapsed time     0.01 sec
  
   - Calculate total internal internal energy and Helmholtz free energy
     Finished, elapsed time     0.00 sec
  
  All done!
  
  Total calculation time:   0.04 sec
  #---------------------------HA calculations finished---------------------------#

As it can see, it is usually a very fast procedure, depending on the computing capabilities
of the computer where Quantas is installed.

Finally, the software tells us where the results were saved, in ``HDF5`` binary format.

.. code-block:: console

  Calculated data exported to mgo_b3lyp_qha_HA.hdf5
  
You can analyze these results by accessing the ``HDF5`` file with the method of your choice. For simplicity, Quantas is able to export these data in a human-readable format by employing the following command:

.. code-block:: console

  > quantas export ha mgo_b3lyp_qha_HA.hdf5

In this way, you enter in a interactive shell of Quantas.

.. code-block:: console
  
  ________                       __
  \_____  \  __ _______    _____/  |______    ______
   /  / \  \|  |  \__  \  /    \   __\__  \  /  ___/
  /   \_/.  \  |  // __ \|   |  \  |  / __ \_\___ \
  \_____\ \_/____/(____  /___|  /__| (____  /____  >
         \__>          \/     \/          \/     \/
                                              v0.9.0
  Authors: Gianfranco Ulian and Giovanni Valdre'
  Copyright (c) Gianfranco Ulian and Giovanni Valdre'.
  
  
  This file was created with Quantas.
  
  Job: Quasi-Harmonic analysis of periclase (MgO)
  
  It contains the results of the harmonic approximation (HA) calculations
  performed with the following settings:
    - energy scale: Ha
    - volume scale: A
    - frequency scale: cm^-1
    - temperature scale: K
  
  Thermodynamics properties were calculated as a function of volume and
  temperature.
  
  
  Select data to export (Cv, F, Fvib, S, U0, Uth, Utot, Uzp):

You can select any of the properties listed in the prompt (case-sensitive).
For example, if you want to extract the entropy data, you should write:

.. code-block:: console

  Select data to export (Cv, F, Fvib, S, U0, Uth, Utot, Uzp): S
  
and give Return. Now, the software will ask you the name of the output file:

.. code-block:: console

  Select data to export (Cv, F, Fvib, S, U0, Uth, Utot, Uzp): S
  Data successfully exported to 'mgo_b3lyp_qha_HA_S.dat'
  
In the case the ``mgo_b3lyp_qha_HA_S.dat`` file is present in the work folder, Quantas will 
ask if it can overwrite the file. If you do not want to overwrite the previous file (thus 
prompting "no"), it will ask for a file name:

.. code-block:: console

  File 'mgo_b3lyp_qha_HA_S.dat' exists. Overwrite it? [y/N]: n
  Please, enter a file name:
  
We are fine with the harmonic entropy data, but you are free to extract any other information
that you would like to see.

If you open the entropy file generated, ``mgo_b3lyp_qha_HA_S.dat``, you can see the results 
reported in a table-like format.

.. code-block:: console

  Entropy
  Data in mHa units

  +----------------------------------------------------- [...] ----------------------------------------------+
     T (K)                                            Volume (A^3)                                                                                                              
                  15.49542672          16.06801923       [...]          20.47572337          21.16448527      
  +===================================================== [...] ==============================================+
      0.00     0.000000000000E+00   0.000000000000E+00   [...]       0.000000000000E+00   0.000000000000E+00  
     10.00     1.977421580893E-17   2.747911180059E-17   [...]       2.764828172628E-15   4.551017733552E-14  
     20.00     4.754272242107E-10   5.577188406358E-10   [...]       7.829027654096E-09   2.616193123357E-08  
     30.00     1.189019045599E-07   1.325634760029E-07   [...]       1.042154363874E-06   2.202731892808E-06  
     40.00     1.810739753020E-06   1.988267139183E-06   [...]       1.208114541732E-05   2.087249201741E-05  
     50.00     9.346012248791E-06   1.026785196648E-05   [...]       5.356262336294E-05   8.273904658350E-05  
     ...
    1990.00    3.770624467366E-02   3.847324899680E-02   [...]       4.471287080910E-02   4.579122110104E-02  
    2000.00    3.779881639452E-02   3.856591829512E-02   [...]       4.480604835991E-02   4.588444880471E-02  
  +===================================================== [...] ==============================================+
  
Just some columns and rows were shown above for the sake of clarity.

These results could be used for creating bi-dimensional plot of the harmonic properties at 
selected volume (for isochoric sections) or temperature (isothermal sections). Also, they
could be used even for three-dimensional plots or bi-dimensional *V-T* contour maps of 
selected properties.

------

Quasi-Harmonic approximation
============================

The same input file for :download:`periclase (MgO)<../downloads/mgo_b3lyp_qha.yaml>` can be 
employed to perform a quasi-harmonic approximation (QHA) analysis of the mineral properties 
at different *P-T* conditions.

In this case, we can set also a pressure range over which the thermodynamic and thermomechanic
properties are calculated. 

.. warning::
  Depending on the different QHA scheme that you use during a Quantas run, it may be better 
  to provide a pressure range within the one explored in static (0 K) conditions.

At the beginning of the analysis, Quantas provides the static pressure valued for each unit-
cell volume in the input file, calculated according to the selected volume minimization scheme.

.. seealso::

  `Volume minimization methods coded in Quantas <../usage/qha_options.html#volume-minimization-scheme>`_.
  
For the sake of an example, let's run a QHA calculation using the default values:

.. code-block:: console

  > quantas qha mgo_b3lyp_qha.yaml
  
The output printed on screen will be:

.. code-block:: console
  
  ________                       __
  \_____  \  __ _______    _____/  |______    ______
   /  / \  \|  |  \__  \  /    \   __\__  \  /  ___/
  /   \_/.  \  |  // __ \|   |  \  |  / __ \_\___ \
  \_____\ \_/____/(____  /___|  /__| (____  /____  >
         \__>          \/     \/          \/     \/
                                              v0.9.0
  Authors: Gianfranco Ulian
  Copyright 2020, University of Bologna
  
  
  Calculator: Quasi-Harmonic Approximation
  
  Quasi-Harmonic Approximation approach
  -------------------------------------
   - scheme: thermodynamics interpolation
   - volume minimization: polynomial functions
  
  Polynomial fitting settings (degree)
  -------------------------------------
   - energy:    3
   - frequency: 3
  
  Temperature settings
  -------------------------------------
   - min:    298.15 K
   - max:    298.15 K
   - step:     1.00 K
  
  Pressure settings
  -------------------------------------
   - min:      0.00 GPa
   - max:      0.00 GPa
   - step:     1.00 GPa
  
  Measurement units
  -------------------------------------
   - energy:      Ha
   - lenght:      A
   - frequency:   cm^-1
   - temperature: K
   - pressure:    GPa

This tells us that for this quasi-harmonic approximation calculation:

 - Quantas interpolates thermodynamic properties at :math:`V(T,P)` conditions;
 - the minimum :math:`V(T,P)` is found by minimizing Helmholtz free energy 
   curves, :math:`F(T,V)`, by using polynomial functions.

The polynomial functions are of the third degree for both frequency values
and energy values. By default, only data at 298.15 K and 0 GPa will be calculated.

The following lines are similar to those of the harmonic approximation calculation:

.. code-block::

  Reading input file: mgo_b3lyp_qha.yaml
  
  Job: Quasi-Harmonic analysis of periclase (MgO)
  
  System:
  - Number of volumes                        11
  - Number of atoms                          2
  - Number of sampled k-points               32
  - Number of frequencies                    192
  
  Volume and energy values:
  
   Pressure (GPa)   Volume (A^3)            Energy (Ha)          
   --------------- --------------- ------------------------------
            42.917    15.495427         -2.754467055837e+02      
            37.504    16.068019         -2.754523424409e+02      
            27.245    16.654547         -2.754566420093e+02      
            18.320    17.255176         -2.754597306884e+02      
            10.561    17.870076         -2.754617227028e+02      
             3.759    18.499413         -2.754627209751e+02      
            -0.018    18.896862         -2.754628820004e+02      
            -2.115    19.143355         -2.754628184464e+02      
            -7.083    19.802069         -2.754620989310e+02      
           -11.464    20.475723         -2.754606377438e+02      
           -13.516    21.164485         -2.754585023425e+02

but this time the pressure state of the unit-cell is also reported. This table suggests
that a pressure range 0 -- 42 GPa could be a good choice to remain in the interpolation
regime.

.. note::
  
  You should have noted that a complete run was made, and a ``HDF5`` file is now present 
  in the folder, whose name is ``mgo_b3lyp_qha_QHA.hdf5``. This will be discussed in the 
  following.

.. warning:: 

  The calculated pressures can sensibly change by considering different degrees of the
  polynomial functions or if equation of state is employed.
  
For example, let's run again the QHA calculation, but this time we'll employ a phenomenological
approach based on a volume-integrated equation of state (EoS) fit:

.. code-block:: console

  > quantas qha mgo_b3lyp_qha.yaml -m eos
  
By inserting the ``-m eos`` flag, Quantas changes the volume minimization scheme.
The output is consequently different (similar lines were skipped in the snippet below):

.. code-block:: console

  [...]
  Quasi-Harmonic Approximation approach
  -------------------------------------
   - scheme: thermodynamics interpolation
   - volume minimization: Equation of State (EoS)
   - EoS formulation: 3rd-order Birch-Murnaghan
  [...]
  Volume and energy values:
  
   Pressure (GPa)   Volume (A^3)            Energy (Ha)
   --------------- --------------- ------------------------------
            50.282    15.495427         -2.754467055837e+02
            37.892    16.068019         -2.754523424409e+02
            27.298    16.654547         -2.754566420093e+02
            18.209    17.255176         -2.754597306884e+02
            10.389    17.870076         -2.754617227028e+02
             3.644    18.499413         -2.754627209751e+02
            -0.069    18.896862         -2.754628820004e+02
            -2.186    19.143355         -2.754628184464e+02
            -7.234    19.802069         -2.754620989310e+02
           -11.611    20.475723         -2.754606377438e+02
           -15.407    21.164485         -2.754585023425e+02
  
  EoS fitting parameters for static energy vs volume data:
  E0 = -275.462883(1)
  K0 = 167.8(7)
  K' = 3.897(7)
  V0 = 18.8891(4)
  
It is possible to note that, by using the EoS fitting method, slightly different pressure values 
at each unit-cell volume were obtained. With this minimization approache, it could be possible to
run the QHA calculation by considering pressures between 0 GPa and 50 GPa.

:Excercise:

  Try to calculate the pressure state of periclase by using other EoS formulations. Do you see
  any difference?
  
Now that we have seen how to it is possible to set the pressure range of the crystal under 
analysis, we can set up a complete calculation. We will consider the same temperature range 
previously employed for the harmonic approximation analysis (0 -- 2000 K, increment of 10 K) 
and we will perform the QHA calculation between 0 and 10 GPa, with a pressure step of 2 GPa.

As volume minimization scheme, we'll consider the phenomenological equation of state, using a 
third-order Birch-Murnaghan formulation. Thermodynamic properties will be calculated with the 
frequency interpolation scheme.

In addition, since the phonon frequency values were previously sorted, it is possible to employ the
frequency interpolation scheme to obtain thermodynamic quantities.

We translate these choices in Quantas by typing:

.. code-block::

  > quantas qha mgo_b3lyp_qha.yaml -s freq -m eos --eos BM -T 0 2000 10 -P 0 10 2

After the initial stream of settings and input file data, Quantas reports some information on 
the frequency fitting procedure:

.. code-block::

  #---------------Quasi-Harmonic Approximation calculation started---------------#

   - Fitting frequency using polynomial of order 3
   - Band # 0
     * Frequency     1: R^2 = 0.000000    BAD!
     * Frequency     2: R^2 = 0.000000    BAD!
     * Frequency     3: R^2 = 0.000000    BAD!
     * Frequency     4: R^2 = 0.999877    OK
     * Frequency     5: R^2 = 0.999877    OK
     * Frequency     6: R^2 = 0.999877    OK
     Averaged R^2: 0.499939
  
   - Band # 1
     * Frequency     1: R^2 = 0.999876    OK
     * Frequency     2: R^2 = 0.999876    OK
     * Frequency     3: R^2 = 0.999876    OK
     * Frequency     4: R^2 = 0.999876    OK
     * Frequency     5: R^2 = 0.999876    OK
     * Frequency     6: R^2 = 0.999876    OK
     Averaged R^2: 0.999876
  
   - Band # 2
     * Frequency     1: R^2 = 0.999876    OK
     * Frequency     2: R^2 = 0.999876    OK
     * Frequency     3: R^2 = 0.999876    OK
     * Frequency     4: R^2 = 0.999940    OK
     * Frequency     5: R^2 = 0.999940    OK
     * Frequency     6: R^2 = 0.999940    OK
     Averaged R^2: 0.999908
     
     [... omitted lines ...]
     
   - Band # 31
     * Frequency     1: R^2 = 0.999996    OK
     * Frequency     2: R^2 = 0.999996    OK
     * Frequency     3: R^2 = 0.999996    OK
     * Frequency     4: R^2 = 0.999996    OK
     * Frequency     5: R^2 = 0.999996    OK
     * Frequency     6: R^2 = 0.999996    OK
     Averaged R^2: 0.999996

   Operation time         125.777 msec
        
Besides the first three phonon modes at Band 0 (namely :math:`\Gamma`-point acoustic modes,
which have a value of 0 :math:`cm^{-1}` at each unit-cell volume), the vibrational frequencies
were all well fitted.

After this report, the thermodynamic and thermoelastic properties are calculated:

.. code-block::

   - Calculation of harmonic Helmoltz free energy F(V,T)
     Operation time          10.055 msec
  
   - Volume minimization using EoS:
     Operation time        2056.631 msec
  
   - Calculation of quasi-harmonic thermodynamic properties
     Operation time         907.329 msec
  
   - Calculation of enthalpy H(T,P) and Gibbs free energy G(T,P):
      * enthalpy                           - elapsed time           0.471 msec
      * Gibbs free energy                  - elapsed time           0.072 msec
     Operation time           0.989 msec
  
   - Calculation of volumetric thermal expansion coefficient:
     Operation time           0.322 msec
  
    - Calculation of isobaric heat capacity
     Operation time           0.447 msec
  
    - Calculation of adiabatic bulk modulus K_S(P,T)
     Operation time           0.401 msec
  
    - Calculation of Gruneisen parameters
     Operation time           0.392 msec
  
  
  Total wall time:   3.11 sec
  #----------------------------QHA Calculation ended-----------------------------#

By using the frequency interpolation scheme, Quantas calculate first the Helmholtz free
energy of the mineral as :math:`F(T,V)`, then is finds the unit-cell volume that minimizes
the equation of state function at selected pressure and temperature conditions. In addition,
the bulk modulus, :math:`K_T`, and its pressure derivative, :math:`K'` are contextually 
obtained by using the EoS approach.

As third step, thermodynamic properties are calculated at *P-T* conditions from the 
calculated :math:`V(P,T)` and exploiting the frequency *vs* volume, :math:`\nu (V)`,
continuity.

.. note::

  For long-running operations, a progress bar shows up.

Then, the volumetric thermal expansion coefficient, :math:`\alpha_V(P,T)` and the other
elastic and thermodynamic properties are calculated.

.. warning::

  The volumetric thermal expansion coefficient, and associated properties, can be calculated 
  only if unit-cell volume is calculated in a temperature range comprising at least 50 points.
  If less temperature points are considered, this calculation is skipped and a warning is 
  issued.
  
Finally, a brief report of the most interesting quasi-harmonic results is printed (and also
reported in the log file).

As previously explained for the harmonic approximation calculation, it is possible to export
the properties from the binary ``HDF5`` output file to a table-like text file using:

.. code-block:: console

  > quantas export qha hdf5_file_name.hdf5
  