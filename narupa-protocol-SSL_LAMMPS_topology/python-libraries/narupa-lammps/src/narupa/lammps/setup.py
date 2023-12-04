#!/usr/bin/env python
# Copyright (c) University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
import numpy

from setuptools import find_namespace_packages
from Cython.Build import cythonize
from setuptools import Extension, setup

extensions = [
    Extension("accelerated_bonds",
              ["accelerated_bonds.pyx"],
              extra_compile_args=['-fopenmp','-O3', '-ffast-math'],
              extra_link_args=['-fopenmp'],),

]

setup(
      ext_modules=cythonize(extensions), include_dirs=[numpy.get_include()],

      )
