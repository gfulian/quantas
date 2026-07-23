Installation
============

Requirements
------------

Quantas requires Python 3.10 or newer. The scientific runtime depends on
NumPy, SciPy, h5py, PyYAML, spglib, odrpack, Click, and Rich. Both spglib and
odrpack are installed with the base package.

Install from the source tree
----------------------------

From the root directory of a Quantas source checkout:

.. code-block:: console

   python -m pip install .

Optional plotting support is installed with:

.. code-block:: console

   python -m pip install ".[plot]"

For development, testing, and documentation:

.. code-block:: console

   python -m pip install -e ".[dev]"
   python -m pip install -e ".[docs]"

Verify the installation
-----------------------

.. code-block:: console

   quantas -v
   quantas --help

The Python API can be checked with:

.. code-block:: python

   import quantas
   from quantas.api import elasticity, registry

   print(quantas.__version__)
   print(elasticity.run)
   print(registry.get("elasticity").capabilities)

Building the documentation
--------------------------

On Linux and macOS:

.. code-block:: console

   make -C docs html

On Windows Command Prompt or PowerShell:

.. code-block:: doscon

   .\docs\make.bat html

The generated site is written to ``docs/_build/html``.
