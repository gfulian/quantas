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


with open("README.md", "r") as fh:
    long_description = fh.read()


requirements = [
    'cython>=0.29',
    'click>=7.0',
    'numpy>=1.18',
    'scipy>=1.4',
    'pyyaml>=5.3',
    'h5py>=2.10'
    ]

packages = [
    'quantas',
    'quantas.cmdline',
    'quantas.cmdline.commands',
    'quantas.cmdline.utils',
    'quantas.core',
    'quantas.eosfit',
    'quantas.eosfit.commands',
    'quantas.eosfit.utils',
    'quantas.harmonic',
    'quantas.harmonic.commands',
    'quantas.harmonic.utils',
    'quantas.interfaces',
    'quantas.IO',
    'quantas.qha',
    'quantas.qha.commands',
    'quantas.qha.utils',
    'quantas.soec',
    'quantas.soec.commands',
    'quantas.soec.utils',
    'quantas.utils',
    'quantas.utils.chemistry',
    'quantas.utils.math',
    'quantas.utils.physics',
    ]

cwd = os.path.join(os.path.dirname(__file__),'quantas')
directories = [
    os.path.join(cwd, ''),
    os.path.join(cwd, 'cmdline'),
    os.path.join(cwd, 'cmdline', 'commands'),
    os.path.join(cwd, 'cmdline', 'utils'),
    os.path.join(cwd, 'core'),
    os.path.join(cwd, 'eosfit'),
    os.path.join(cwd, 'eosfit', 'commands'),
    os.path.join(cwd, 'eosfit', 'utils'),
    os.path.join(cwd, 'harmonic'),
    os.path.join(cwd, 'harmonic', 'commands'),
    os.path.join(cwd, 'harmonic', 'utils'),
    os.path.join(cwd, 'interfaces'),
    os.path.join(cwd, 'IO'),
    os.path.join(cwd, 'qha'),
    os.path.join(cwd, 'qha', 'commands'),
    os.path.join(cwd, 'qha', 'utils'),
    os.path.join(cwd, 'soec'),
    os.path.join(cwd, 'soec', 'commands'),
    os.path.join(cwd, 'soec', 'utils'),
    os.path.join(cwd, 'utils'),
    os.path.join(cwd, 'utils', 'chemistry'),
    os.path.join(cwd, 'utils', 'math'),
    os.path.join(cwd, 'utils', 'physics'),
    ]

dirs = {}
for i in range(len(packages)):
    dirs[packages[i]] = directories[i]

# Compiler flags for 'cythonized' modules
if sys.platform == 'win32':
    extra_compile_args = ['/openmp']
    extra_link_args = []
elif sys.platform == 'linux':
    extra_compile_args = ['-fopenmp']
    extra_link_args = ['-fopenmp']
else:
    extra_compile_args = ['-fopenmp']
    extra_link_args = ['-lomp']

setup(name='quantas',
      version='0.9.0',
      description='QUANtistic Thermomechanical Analysis of Solids',
      long_description=long_description,
      classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: BSD License',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.6',
        'Operating System :: OS Independent',
        'Environment :: Console',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',
        'Intended Audience :: Science/Research'
      ],
      url='',
      author='Gianfranco Ulian',
      author_email='gianfranco.ulian2@unibo.it',
      license='MIT',
      package_dir=dirs,
      packages=packages,
      entry_points={
          'console_scripts': [
              'quantas = quantas.cmdline.commands.cmd_quantas:cli'
          ]
      },
      ext_modules=cythonize([
          Extension('quantas.utils.physics.statistical_mechanics',
                    ['quantas/utils/physics/statistical_mechanics.pyx'],
                    extra_compile_args=extra_compile_args,
                    extra_link_args=extra_link_args,
                    ),
          Extension('quantas.utils.physics.thermodynamics',
                    ['quantas/utils/physics/thermodynamics.pyx'],
                    extra_compile_args=extra_compile_args,
                    extra_link_args=extra_link_args,
                    ),
          Extension('quantas.utils.math.fast_math',
                    ['quantas/utils/math/fast_math.pyx'],
                    extra_compile_args=extra_compile_args,
                    extra_link_args=extra_link_args,
                    ),
          ]),
      python_requires='>=3.6',
      install_requires=requirements,
      include_package_data=True,
      zip_safe=False)
