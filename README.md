# README
#
# Quantas (QUANTitative Analysis of Solids) is a software 
# developed to aid the analysis of data obtained from (ab initio) 
# quantum mechanical simulations and/or experiments. 
# It is a cross-platform program coded in Python 3.x. 
#
# Copyright (c) 2020, Gianfranco Ulian
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
# 
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
# * Neither the name of the Quantas project nor the
#   names of its contributors may be used to endorse or promote products
#   derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL GIANFRANCO ULIAN 
# BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
Requirements
------------
- Python_ 3.6+
- Cython_ 0.29+
- click_ 7.0+
- NumPy_ (base N-dimensional array package)
- SciPy_ (library for scientific computing)
- HDF5_ (library for HDF5)
- YAML_ (library for data serialization)

Optional packages
-----------------

- Matplotlib_ (used for plotting results)

.. _Python: http://www.python.org/
.. _Cython: https://cython.org/
.. _click: https://click.palletsprojects.com/en/7.x/
.. _NumPy: http://docs.scipy.org/doc/numpy/reference/
.. _SciPy: http://docs.scipy.org/doc/scipy/reference/
.. _HDF5: https://www.h5py.org/
.. _YAML: https://pypi.org/project/PyYAML/
.. _Matplotlib: https://matplotlib.org/

.. note::

  On Windows, you could benefit from the 
  `Unofficial Windows Binaries <http://www.lfd.uci.edu/~gohlke/pythonlibs/>`_
  made available Christoph Gohlke, in particular:

  - Numpy+MKL_
  - Scipy_Numpy+MKL_

  After downloading the package wheels, you can install them using ``pip``.

.. _Numpy+MKL: http://www.lfd.uci.edu/~gohlke/pythonlibs/#numpy
.. _Scipy_Numpy+MKL: http://www.lfd.uci.edu/~gohlke/pythonlibs/#scipy

Python virtual environment
==========================

Since Quantas depends on several third-party python packages, we strongly recommend using a 
`virtual environment <https://docs.python.org/3/tutorial/venv.html>`_ to prevent any 
interference with you current working environment. The following instructions are focused on
the `virtualenv <https://virtualenv.pypa.io/en/latest/>`_ tool, but you can freely use any 
other environment manager of your choice.

Create the virtual environment
------------------------------

The creation of a virtual environment using virtualenv is slightly different on Windows and 
Unix machines.

:Linux and Mac OS:

  On a shell (or terminal), issue the following commands::
  
    $ pip install --user --upgrade virtualenv
    $ python3 -m venv ~/.virtualenvs/quantas
  
  Using the above commands, you installed/upgraded virtualenv and created home directory of 
  the virtual environment. If you look in the home folder under ``.virtualenv``, you can find
  a directory named ``quantas``.
  
  To activate the virtual environment, just issue::
  
    $ source ~/.virtualenvs/quantas/bin/activate
  
  After activation, your prompt should appear as::
  
    (quantas) $
  
  indicating that you are working on a virtual environment. Now, the python executable of the 
  virtualenv is the first in PATH, and that python programs have access only to packages 
  installed inside the virtualenv.
  
  To leave or deactivate the ``quantas`` virtual environment, run::
  
    (quantas) $ deactivate
  
:Windows:

  On a command prompt, issue the following command to install/upgrade the virtualenv package::
  
    > pip install --user --upgrade virtualenv

  
  Then, in a folder of your choice, run the venv command as::
  
    > python -m venv quantas
    
  or::
  
    > virtualenv quantas
    
  to create the home folder of the virtual environment. To activate it, issue::
  
    > quantas\Scripts\activate.bat
  
  Differently from Linux and Mac OS, no clear indication on the successful activation of the 
  environment is provided. However, you can run the following::
  
    > python
    >>> import sys
    >>> sys.path
  
  and check if in the output the path to the virtual environment directory is the first one 
  provided.
  
  To leave or deactivate the ``quantas`` virtual environment, run::
  
    > deactivate


Installation from source
------------------------

At the moment, it is only possible to install Quantas from its 
:download:`source code <./downloads/quantas-0.9.0.tar.gz>`. Common alternatives such 
as ``pip`` will be available in the next future.

:Linux and Mac OS:

  On a shell (or terminal), unpack the compressed file and enter in the decompressed directory:
  
  .. code-block:: console
   
    $ tar -xf quantas-0.9.0.tar.gz
    $ cd quantas-0.9.0
    
  Then, install the Quantas package:
  
  .. code-block:: console
   
    $ python3 setup.py install


:Windows:

  Use a software as 7zip or WinRar to decompress the file. On a command prompt, enter the 
  directory:
    
  .. code-block:: console
 
    > cd quantas-0.9.0
    
  Then, install the Quantas package:
  
  .. code-block:: console
   
    > python setup.py install
    

.. note::

    On Linux, you may need root privileges to install the package in the ``/usr/local``
    directory (default). If you prefer a local (user) installation, you could use the 
    :envvar:`--prefix` directive to specify a different location.
    In this case, the launching scripts will be installed in the ``~/.local/bin`` directory of 
    the user.


Environment variables
---------------------

If you installed Quantas in a system-wide fashion, please ensure that the following variables 
are set.

.. envvar:: PATH

    Colon-separated paths where programs can be found.

.. envvar:: PYTHONPATH

    Colon-separated paths where Python modules can be found.

Under Linux, you can set these permanently in your :file:`~/.bashrc` file::

    $ export PYTHONPATH=<path-to-Quantas-package>:$PYTHONPATH
    $ export PATH=<path-to-Quantas-command-line-tools>:$PATH

or your :file:`~/.cshrc` file::

    $ setenv PYTHONPATH <path-to-Quantas-package>:${PYTHONPATH}
    $ setenv PATH <path-to-Quantas-command-line-tools>:${PATH}

.. note::

   If running on Mac OSX: be aware that terminal sessions will
   source :file:`~/.bash_profile` by default and not
   :file:`~/.bashrc`. Either put any ``export`` commands into
   :file:`~/.bash_profile` or source :file:`~/.bashrc` in all Bash
   sessions by adding

   ::

      if [ -f ${HOME}/.bashrc ]; then
      source ${HOME}/.bashrc
      fi

   to your :file:`~/.bash_profile`.

.. note::

   Under Windows, the environmental variables should have been set during/after the 
   installation of the Python 3.x package.


Test Quantas installation
-------------------------

Some :download:`input examples <./downloads/examples.zip>` are provided to test 
the installation of Quantas. Unpack them in any folder you like.

Before running the tests, make sure you have set your :envvar:`PATH`.

Run the tests like:

.. code-block:: console

    $ quantas ha examples\mgo_b3lyp_qha.yaml
    $ quantas qha examples\mgo_b3lyp_qha.yaml
    $ quantas eos examples\PV_topaz.dat
    $ quantas soec examples\hydroxylapatite.dat

If something goes wrong, please send us an output log of the failing test.