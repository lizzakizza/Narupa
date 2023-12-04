# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
#!/usr/bin/env python

from distutils.core import setup
from setuptools import find_namespace_packages

setup(
    name='narupa-openmm',
    version='0.1.0',
    description='OpenMM server for Narupa',
    author='Intangible Realities Lab',
    author_email='m.oconnor@bristol.ac.uk',
    url='https://gitlab.com/intangiblerealities/',
    packages=find_namespace_packages('src', include='narupa.*'),
    package_dir={'': 'src'},
    package_data={
        '': ['py.typed']
    },
    install_requires=(
        'narupa',
        'openmm',
    ),
    entry_points={
        'console_scripts': ['narupa-omm-server=narupa.openmm.cli:main'],
    },
)
