Requirements
============
- Python 3.6+ (http://www.python.org/)
- Cython 0.29+ (https://cython.org/)
- click 7.0+ (https://click.palletsprojects.com/en/7.x/)
- NumPy (http://docs.scipy.org/doc/numpy/reference/)
- SciPy (http://docs.scipy.org/doc/scipy/reference/)
- HDF5 (https://www.h5py.org/)
- YAML (https://pypi.org/project/PyYAML/)

Optional packages
-----------------

- Matplotlib (https://matplotlib.org/)

On Windows, you could benefit from the Unofficial Windows Binaries 
(http://www.lfd.uci.edu/~gohlke/pythonlibs/) made available by Christoph 
Gohlke, in particular:

  - Numpy+MKL (http://www.lfd.uci.edu/~gohlke/pythonlibs/#numpy)
  - Scipy_Numpy+MKL (http://www.lfd.uci.edu/~gohlke/pythonlibs/#scipy)

After downloading the package wheels, you can install them using ``pip``.

Python virtual environment
==========================

Since Quantas depends on several third-party python packages, we strongly recommend using a 
virtual environment (https://docs.python.org/3/tutorial/venv.html) to prevent any 
interference with you current working environment. The following instructions are focused on
the virtualenv (https://virtualenv.pypa.io/en/latest/) tool, but you can freely use any 
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
========================

At the moment, it is only possible to install Quantas from its source code. Common 
alternatives such as ``pip`` will be available in the next future.

Linux and Mac OS
----------------

  On a shell (or terminal), unpack the compressed file and enter in the decompressed directory::
   
    $ tar -xf quantas-0.9.0.tar.gz
    $ cd quantas-0.9.0
    
  Then, install the Quantas package::
   
    $ python3 setup.py install

On Linux, you may need root privileges to install the package in the ``/usr/local``
directory (default). If you prefer a local (user) installation, you could use the 
`--prefix` directive to specify a different location.
In this case, the launching scripts will be installed in the ``~/.local/bin`` directory of 
the user.

Windows
-------

  Use a software as 7zip or WinRar to decompress the file. On a command prompt, enter the 
  directory:
    
  .. code-block:: console
 
    > cd quantas-0.9.0
    
  Then, install the Quantas package:
  
  .. code-block:: console
   
    > python setup.py install

	
Test Quantas installation
=========================

Some input examples are provided to test the installation of Quantas. 
Unpack them in any folder you like.

Before running the tests, make sure you have set your `PATH`, so that the Quantas 
executable can be accessed.

Run the tests like::

    $ quantas ha examples\mgo_b3lyp_qha.yaml
    $ quantas qha examples\mgo_b3lyp_qha.yaml
    $ quantas eos examples\PV_topaz.dat
    $ quantas soec examples\hydroxylapatite.dat

If something goes wrong, please send us an output log of the failing test.