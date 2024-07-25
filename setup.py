# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c), Gianfranco Ulian and Giovanni Valdre'.                      #
# All rights reserved.                                                       #
#                                                                            #
# This file is part of the Quantas code.                                     #
#                                                                            #
# For further information on the license, see the LICENSE file               #
##############################################################################

import os
import sys
from setuptools import setup
from distutils.core import Extension

try:
    from Cython.Build import cythonize
    from Cython.Distutils import build_ext
except ImportError:
    print('\nError: Cython package not found')
    print('\nCython is required to prooceed with the installation of Quantas')
    print("\nPlease, install Cython via 'pip install cython'")
    print('before installing this package')
    print('\nWill now exit')
    sys.exit(0)

# ----------------------------------------------------------------
# Setup utilities
# ----------------------------------------------------------------
def convert_path(string):
    """ Convert a string containing a path into a string containing a module."""
    if sys.platform == 'win32':
        new_string = string.replace('\\','.')
    else:
        new_string = string.replace('/','.')
    return new_string

def get_quantas_package(root_dir, cwd):
    """ Collect all the modules within the Quantas package. """
    package = {}
    modules = []
    modules.append(root_dir)
    paths = []
    paths.append(os.path.join(cwd, root_dir))
    for root, dirs, files in os.walk(root_dir):
        for directory in dirs:
            modules.append(os.path.join(root, directory).replace('\\','.'))
            paths.append(os.path.join(cwd, root, directory))
    modules.sort()
    paths.sort()
    for i in range(len(modules)):
        package[modules[i]] = paths[i]
    return package, modules, paths


def get_extension_modules(root_dir):
    """ Collect all the cython extensions. """
    extensions = []
    # ------------------------------------------------------------
    # Compiler flags for 'cythonized' modules
    # ------------------------------------------------------------
    if sys.platform == 'win32':
        extra_compile_args = ['/openmp']
        extra_link_args = []
    elif sys.platform == 'linux':
        extra_compile_args = ['-fopenmp']
        extra_link_args = ['-fopenmp']
    else:
        extra_compile_args = []
        extra_link_args = []

    # ------------------------------------------------------------
    # Start collection
    # ------------------------------------------------------------
    for root, dirs, files in os.walk(root_dir):
        for filename in files:
            if os.path.splitext(filename)[1] == '.pyx':
                module = convert_path(
                    os.path.splitext(os.path.join(root, filename))[0]
                    )
                path = [os.path.join(root, filename)]
                extensions.append(
                    Extension(module, path,
                              extra_compile_args=extra_compile_args,
                              extra_link_args=extra_link_args
                              )
                    )
    return extensions

# ----------------------------------------------------------------
# Quantas description
# ----------------------------------------------------------------
with open("README.md", "r") as fh:
    long_description = fh.read()

# ----------------------------------------------------------------
# Package requirements
# ----------------------------------------------------------------
requirements = [
    'cython>=0.29',
    'click>=7.0',
    'numpy>=1.18',
    'scipy>=1.4',
    'pyyaml>=5.3',
    'h5py>=2.10'
    ]

# ----------------------------------------------------------------
# Package modules and folders
# ----------------------------------------------------------------
cwd = os.path.join(os.path.dirname(__file__))
package, modules, paths = get_quantas_package('quantas', '.')

# ----------------------------------------------------------------
# Compiler flags for 'cythonized' modules
# ----------------------------------------------------------------
if sys.platform == 'win32':
    extra_compile_args = ['/openmp']
    extra_link_args = []
elif sys.platform == 'linux':
    extra_compile_args = ['-fopenmp']
    extra_link_args = ['-fopenmp']
else:
    extra_compile_args = []
    extra_link_args = []

# ----------------------------------------------------------------
# Setup
# ----------------------------------------------------------------
setup(name='quantas',
      version='0.9.1',
      description='QUANtistic Thermomechanical Analysis of Solids',
      long_description=long_description,
      classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: BSD License',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Operating System :: OS Independent',
        'Environment :: Console',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',
        'Intended Audience :: Science/Research'
      ],
      url='https://github.com/gfulian/quantas',
      author='Gianfranco Ulian',
      author_email='gianfranco.ulian2@unibo.it',
      license='BSD',
      package_dir=package,
      packages=modules,
      entry_points={
          'console_scripts': [
              'quantas = quantas.cmdline.commands.cmd_quantas:cli'
          ]
      },
      ext_modules=cythonize(get_extension_modules('quantas')),
      python_requires='>=3.5',
      install_requires=requirements,
      include_package_data=True,
      zip_safe=False)
