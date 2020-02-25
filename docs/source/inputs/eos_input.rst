.. _eos_input:

=========================================
Input for Equation of State (EoS) fitting
=========================================

  :Last updated: |today|
  :Author: **Gianfranco Ulian**

The input file for the equation of state fitting can be of any extension, albeit ``.dat`` 
should be preferred. 

The input file format for equation of state fitting closely resembles that of the EOSFit7 
software. [1]_ An example is reported in the following.

.. code::

    #
    # Comment lines start with a '#' character and are ignored by the file parser.
    #
    JOB
    Quartz (Angel et al., 2000)
    FORMAT
    P V sigmaP sigmaV
    DATA
    0.001  112.981   1.000E-6 0.002
    0.430  111.725   0.008    0.014
    0.794  110.711   0.008    0.006
    1.651  108.597   0.007    0.007
    1.845  108.150   0.007    0.019
    1.934  107.974   0.007    0.008
    2.629  106.467   0.012    0.008
    3.300  105.141   0.007    0.009
    3.480  104.831   0.012    0.012
    3.780  104.285   0.012    0.010
    4.026  103.810   0.009    0.013
    4.554  102.939   0.009    0.007
    4.819  102.534   0.012    0.012
    5.214  101.933   0.010    0.010
    5.416  101.592   0.009    0.008
    5.737  101.143   0.008    0.011
    6.478  100.051   0.012    0.015
    6.751   99.677   0.010    0.009
    7.192   99.064   0.014    0.012
    7.899   98.204   0.009    0.014
    8.449   97.545   0.012    0.016
    8.905   96.989   0.009    0.017

In the following table, the keywords employed by for EoS analysis are reported, alongside a
short description of their functionality.

===================== ===================================================
*Mandatory Keywords*

**FORMAT**            Specify what the different data columns represent
**DATA**              Actual data that will be collected

*Optional Keywords*
**JOB**               Short description of the provided data
===================== ===================================================


Keywords description
====================

A detailed description of the data requested for each keyword is reported in the following.


FORMAT
------

This is keyword is used to specify how the data is reported in the file, namely what each 
column represents. **FORMAT** is followed in the next line by *sub-keywords* indicating the 
type of data in the input file. Each sub-keyword is separated by at least one white space.

The sub-keywords of **FORMAT** are showed in the following table.

========================= ===================================================
**FORMAT sub-keywords**   **Related quantity**

*P*                       Pressure
*V*                       Unit cell volume
*a*                       **a**-axis values
*b*                       **b**-axis values
*c*                       **c**-axis values

*sigmaP*                  Uncertainties of pressure data
*sigmaV*                  Uncertainties of volume data
*sigmaa*                  Uncertainties of **a**-axis data
*sigmab*                  Uncertainties of **b**-axis data
*sigmac*                  Uncertainties of **c**-axis data
========================= ===================================================

It is expected that pressure data is present, thus the *P* sub-keyword is mandatory. Also, at
least one of the structural information, either volume or axis lenghts, must be present, otherwise **Quantas** will throw an error.

It is not important the order of the sub-keywords, it is sufficient that the structural data 
are reported in the same order as in **FORMAT**.

.. note::

    Please, do not insert comments between **FORMAT** and the next line.
    
.. warning::

    **FORMAT** must be inserted in the input file *before* **DATA**, otherwise the software is 
    not able to collect the quantities.
    
:Example:
    
    If your input data consists of unit cell volume and **a**-axis lenghts, and you have the
    uncertainties only for the latter, you should write:
    
    .. code::
    
      [...]
      FORMAT
      P V a sigmaa
      [...]


DATA
----

This keyword informs Quantas that the next block is composed by pressure and unit cell data, 
as specified by the **FORMAT** keyword. 

.. warning::

    **FORMAT** must be inserted in the input file *after* **FORMAT**, otherwise the software is 
    not able to collect the quantities.
   

JOB (optional)
--------------

This keyword sets some information on the input file, that will be printed by **Quantas** 
during the run. This description have to be inserted in the line following this keyword. 

.. note::

  Please, do not insert comments between **JOB** and the next line, otherwise the comment 
  will be employed as data description!


.. rubric:: References

.. [1] Angel, R.J., 2010. EOSFit7. Computer Program (http://www.rossangel.com).