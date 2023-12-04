#!/usr/bin/env python
# Copyright (c) University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

from distutils.core import setup
from distutils.extension import Extension
import numpy

from setuptools import find_namespace_packages
from Cython.Build import cythonize
from setuptools import Extension, setup, command

extensions = [
       Extension("narupa.lammps.accelerated_bonds",
                ["src/narupa/lammps/accelerated_bonds.pyx"],
             extra_compile_args=['-fopenmp', '-O3', '-ffast-math'],
             extra_link_args=['-fopenmp'],),
]

setup(name='narupa-lammps',
      version='0.1.0',
      description='LAMMPS integration for Narupa',
      author='Intangible Realities Lab',
      author_email='simonbennie@gmail.com',
      url='https://gitlab.com/intangiblerealities/',
      packages=find_namespace_packages('src', include='narupa.*'),
      package_dir={'': 'src'},
      package_data={
          '': ['py.typed']
      },
      install_requires=(
            'narupa',
            'mpi4py',
            'numpy',
            'cython'
      ),
      ext_modules=cythonize(extensions), include_dirs=[numpy.get_include()],

      )
