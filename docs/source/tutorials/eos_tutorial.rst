.. _eos_tutorial:

========================================
Equation of State (EoS) fitting tutorial
========================================

  :Last updated: |today|
  :Author: **Gianfranco Ulian**

Preliminary operations
======================

Download the :download:`topaz input file <../downloads/PV_topaz.dat>`, which 
contains experimental *P-V* structural data on topaz, [1]_ and put it in a folder of your choice.

The data reported in the input file are:

  - unit-cell volume
  - *a* lattice parameter
  - *b* lattice parameter
  - *c* lattice parameter

and related uncertainties, measured between 0.0001 GPa and 45.26 GPa. Note that also uncertainty on pressure value is also reported.


Interactive equation of state fitting
=====================================

The EoS fitting procedure is performed in an interactive way by simply typing on the console:

.. code-block:: console

  > quantas eosfit PV_topaz.dat
  
Quantas shows some information collected from the input file:

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
  
  
  Calculator: Equation of State (EoS) fitting
  
  Measurement units
  -------------------------------------
   - pressure:    GPa
  
  Reading input file: PV_topaz.dat
  
  Job: Topaz (Gatta et al., 2014)
  
   P (GPa)     V (A^3)      a (A)       b (A)       c (A)
  ----------  ----------  ----------  ----------  ----------
    0.000     345.46000    4.66270     8.83430     8.38670
    0.410     344.80000    4.65970     8.83030     8.37990
    1.880     341.78000    4.64510     8.81080     8.35090
    2.570     340.35000    4.63780     8.80220     8.33720
    3.140     339.15000    4.63250     8.79380     8.32540
    3.670     338.05000    4.62670     8.78720     8.31500
    4.650     336.28000    4.61790     8.77590     8.29760
    5.160     335.34000    4.61350     8.76920     8.28870
    5.790     334.08000    4.60710     8.76140     8.27650
    6.520     332.82000    4.60090     8.75330     8.26410
    6.620     331.40000    4.57900     8.75100     8.26900
    7.410     331.23000    4.59300     8.74270     8.24880
    8.370     329.45000    4.58410     8.73080     8.23160
    8.970     328.53000    4.57950     8.72470     8.22240
    9.380     326.80000    4.56800     8.71400     8.20900
    11.340    324.10000    4.55900     8.69000     8.18000
    15.270    316.10000    4.52500     8.63600     8.08900
    19.810    309.90000    4.49300     8.58900     8.02900
    26.020    300.00000    4.44600     8.51200     7.92700
    31.490    294.00000    4.41400     8.46400     7.86900
    37.820    289.20000    4.39900     8.41200     7.81400
    41.000    284.00000    4.36400     8.37700     7.76700
    43.840    283.00000    4.36300     8.35600     7.76300
    45.260    281.60000    4.35200     8.34600     7.75300
    45.330    281.70000    4.36000     8.34800     7.73900
  
   - Found weights for pressure
   - Found weights for volume
   - Found weights for a-axis
   - Found weights for b-axis
   - Found weights for c-axis
   
The first thing Quantas asks is what structural property you would like to fit:

.. code-block:: console

  DATA SELECTION
    v. Volume
    a. a-axis
    b. b-axis
    c. c-axis
  
  Select the values to fit:
  
.. note::

  If you select a property not available in input data, Quantas asks you to choose another
  one.
  
Let's start by fitting unit-cell volumes, so select ``v``. Then, we should choose one of the available equation 
of state formulations:

.. code-block:: console

  EOS FORMULATION
    1. Murnaghan
    2. Birch-Murnaghan
    3. Natural Strain
    4. Vinet
    5. Tait
  
  Select an EoS formulation [1]:
  
We will perform a third-order Birch-Murnaghan equation of state fit, thus select ``2`` and, at the next question:

.. code-block:: console

  EOS ORDER
    - 2
    - 3
    - 4
  
  Select the order of the EoS [2]:

choose ``3``. At this point, an initial guess of the EoS parameters is reported:

.. code-block:: console

  EOS PARAMETERS:
    K0       =    187.43910
    K'       =      4.00000
    K''      =      0.00000
    V0       =    345.68229
  
  Would you like to modify any of the EoS parameters? [y/N]:
  
.. note::

  Upper-case letters in the yes/no questions are the default answers. For example, in this 
  question, the default answer is no.
  
Keep these parameters as provided, by selecting 'N' or just pressing Return. At this point, we 
could select some parameters to be kept fixed during the fitting procedure. Let's initially keep 
the :math:`K'` value fixed at 4 (as it would be in a second-order Birch-Murnaghan EoS), so
select ``y`` at the next question and fix the parameters accordingly, as shown also below:

.. code-block:: console

  Would you like to fix any of the EoS parameters? [y/N]: y
    Fix K0     [y/N]: n
    Fix K'     [y/N]: y
    Fix K''    [y/N]: n
    Fix V0     [y/N]: n
    
In the next question, Quantas asks if the uncertainties on both the independent (pressure) and 
dependent (volume, in this case) variables should be used as weights. Type ``y`` here and in the 
subsequent two questions:

.. code-block:: console

  Would you like to weight the data during fitting? [Y/n]: y
   Use weights for pressure? [Y/n]: y
   Use weights for Volume  ? [Y/n]: y

All the questions have been answered. Now, Quantas perform the fitting procedure using the Orthogonal Distance Regression (ODR), then report the results:

.. code-block:: console

  #------------------------------------------------------------------------------#
  #                        EOS results for volume fitting                        #
  #------------------------------------------------------------------------------#
  
  EOS formulation: Birch-Murnaghan
        EOS order: 3
  
  Results weighted for pressure
  Results weighted for volume
  
  Fitting parameters:
    K0     =  149.57006 +/-    1.83762
    K'     =    4.00000 +/-    0.00000 (FIXED)
    K''    =   -0.02600 +/-    0.00000 [IMPLIED]
    V0     =  345.65545 +/-    0.12009
  
  Covariance matrix:
     7.79573e-02  0.00000e+00  0.00000e+00 -2.24409e-03
     0.00000e+00  0.00000e+00  0.00000e+00  0.00000e+00
     0.00000e+00  0.00000e+00  0.00000e+00  0.00000e+00
    -2.24409e-03  0.00000e+00  0.00000e+00  3.32935e-04
  
  Residual variance: 43.3167
  
  Exited ODR loop because:
    1. Problem is not full rank at solution
    2. Sum of squares convergence
  
  Calculated data:
  
      volume           eps_v      Pressure (GPa)    P_obs (GPa)        eps_P
  --------------- --------------- --------------- --------------- ---------------
     345.46000        0.19522         0.00010         0.00010         0.00000
     344.80000       -0.00255         0.41000         0.37358        -0.03642
     341.78000       -0.01099         1.88000         1.72993        -0.15007
     340.35000       -0.01328         2.57000         2.39245        -0.17755
     339.15000       -0.01381         3.14000         2.95848        -0.18152
     338.05000       -0.01422         3.67000         3.48606        -0.18394
     336.28000       -0.02325         4.65000         4.35705        -0.29295
     335.34000       -0.02669         5.16000         4.82819        -0.33181
     334.08000       -0.02643         5.79000         5.46747        -0.32253
     332.82000       -0.03320         6.52000         6.12228        -0.39772
     331.40000        0.23255         6.62000         6.72957         0.10957
     331.23000       -0.08133         7.41000         6.98722        -0.42278
     329.45000       -0.08262         8.37000         7.95160        -0.41840
     328.53000       -0.16541         8.97000         8.50568        -0.46432
     326.80000        0.00658         9.38000         9.38289         0.00289
     324.10000       -0.37938        11.34000        11.18085        -0.15915
     316.10000        0.68924        15.27000        15.53092         0.26092
     309.90000        0.47752        19.81000        19.97426         0.16426
     300.00000        1.85690        26.02000        26.58167         0.56167
     294.00000        1.66534        31.49000        31.94829         0.45829
     289.20000        0.24297        37.82000        37.88074         0.06074
     284.00000        1.95661        41.00000        41.46332         0.46332
     283.00000        0.61572        43.84000        43.98055         0.14055
     281.60000        0.71484        45.26000        45.41988         0.15988
     281.70000        0.57940        45.33000        45.45951         0.12951
  
  #------------------------------------------------------------------------------#
  
Fitting parameters are reported together with their estimated standard deviation and the
covariance matrix.

:code:`Exited ODR loop because` tells us how the fitting procedure proceeded and
exited. The first sentence is related to the fact that the EoS formulation has four parameters,
but we fitted just two of them (:math:`K_0` and :math:`V_0`), keeping one fixed 
(:math:`K'`) and one implied (:math:`K''`) because of the third-order formulation). 
The second sentence, 
'Sum of squares convergence' means that the fitting procedure converged.

.. note::

  As for any fitting procedure, having the initial parameters to values close to the optimized 
  ones may speed up convergence and ensure not diverging (or unphysical) results.

The residual variance is quite high, suggesting that a fixed :math:`K'` (equivalent to a 
second-order Birch-Murnaghan) is not sufficient to fit the data.

After the report, Quantas asks if we want to update the parameters:

.. code-block:: console

  Update the parameters [y/N]: y
  
Let's try improving the fitting, so answer ``y`` and, at the following question:

.. code-block:: console

  Would you like to continue fitting on the same data? [y/N]: y
  
answer ``y``, too.

We could change the EoS formulation in the subsequent point, but there is no reason to do that
at the moment, so answer ``n``:

.. code-block:: console

  Would you like to change the EoS? [y/N]: n
  
Then, Quantas asks us if we want to change the initial parameters, to keep some of them fixed, 
and so on. This time, let everything free to be optimized.

The fitting results now become:

.. code-block:: console

  #------------------------------------------------------------------------------#
  #                        EOS results for volume fitting                        #
  #------------------------------------------------------------------------------#
  
  EOS formulation: Birch-Murnaghan
        EOS order: 3
  
  Results weighted for pressure
  Results weighted for volume
  
  Fitting parameters:
    K0     =  161.71746 +/-    2.53414
    K'     =    3.00184 +/-    0.16220
    K''    =   -0.02404 +/-    0.00000 [IMPLIED]
    V0     =  345.50950 +/-    0.08264
  
  Covariance matrix:
     3.44825e-01 -1.93179e-02  0.00000e+00 -5.61609e-03
    -1.93179e-02  1.41261e-03  0.00000e+00  2.47375e-04
     0.00000e+00  0.00000e+00  0.00000e+00  0.00000e+00
    -5.61609e-03  2.47375e-04  0.00000e+00  3.66721e-04
  
  Residual variance: 18.6235
  
  Exited ODR loop because:
    1. Problem is not full rank at solution
    2. Sum of squares convergence
  
  Calculated data:
  
      volume           eps_v      Pressure (GPa)    P_obs (GPa)        eps_P
  --------------- --------------- --------------- --------------- ---------------
     345.46000        0.04929         0.00010         0.00010         0.00000
     344.80000       -0.00558         0.41000         0.33609        -0.07391
     341.78000       -0.00724         1.88000         1.78741        -0.09259
     340.35000       -0.00623         2.57000         2.49165        -0.07835
     339.15000       -0.00395         3.14000         3.09098        -0.04902
     338.05000       -0.00182         3.67000         3.64769        -0.02231
     336.28000       -0.00724         4.65000         4.56312        -0.08688
     335.34000       -0.00874         5.16000         5.05626        -0.10374
     334.08000       -0.00572         5.79000         5.72310        -0.06690
     332.82000       -0.01009         6.52000         6.40374        -0.11626
     331.40000        0.55134         6.62000         6.87166         0.25166
     331.23000       -0.02537         7.41000         7.28242        -0.12758
     329.45000       -0.01972         8.37000         8.27287        -0.09713
     328.53000       -0.05641         8.97000         8.81543        -0.15457
     326.80000        0.39442         9.38000         9.55036         0.17036
     324.10000        0.03346        11.34000        11.35395         0.01395
     316.10000        1.09195        15.27000        15.68985         0.41985
     309.90000        0.78060        19.81000        20.08821         0.27821
     300.00000        1.84821        26.02000        26.61418         0.59418
     294.00000        1.32555        31.49000        31.88485         0.39485
     289.20000       -0.51291        37.82000        37.67865        -0.14135
     284.00000        0.86232        41.00000        41.22723         0.22723
     283.00000       -0.67228        43.84000        43.66804        -0.17196
     281.60000       -0.70700        45.26000        45.08218        -0.17782
     281.70000       -0.84307        45.33000        45.11804        -0.21196
  
  #------------------------------------------------------------------------------#

The residual variance is better than before, and the :math:`K'` value is less than 4. The 
calculated third-order Birch-Murnaghan EoS parameters are very close to those reported in literature: [1]_

  - :math:`K_0` = 158(4) GPa
  - :math:`K'` = 3.3(3) GPa
  
If you are satisfied with these results, you can exit by answering 'N' to all the subsequent
questions asked by Quantas.

.. note::

  All the fitting results reported on screen will be also available in the output file 
  ``PV_topaz_EOS.log``, generated during the run.
  
:Excercise 1:

  The example file contains also the *a*, *b* and *c* lattice parameters at different
  pressures. Try to fit those data with a third-order Birch-Murnaghan equation of state.
  
:Excercise 2:

  The different EoS formalisms may some important impacts on the final parameters. Try to fit 
  the unit-cell volume, *a*, *b* and *c* lattice parameters using:
  
  - a (second-order) Murnaghan equation of state
  - a third-order Natural Strain equation of state
  - a third-order Tait equation of state
  
  and compare the results.


.. rubric:: References

.. [1] Gatta, G.D., Morgenroth, W., Dera, P., Petitgirard, S., Liermann, H.P., 2014. Elastic
       behavior and pressure-induced structure evolution of topaz up to 45 GPa. Phys. Chem. 
       Miner. 41, 569-577