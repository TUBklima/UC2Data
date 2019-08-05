# -*- coding: utf-8 -*-

from setuptools import setup, find_packages


with open('README.rst') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='uc2data',
    version='0.1.0',
    description='Package for working with netCDF files following the [UC]² data standard',
    long_description=readme,
    author='Achim Holtmann, Tom Grassmann',
    author_email='achim.holtmann@tu-berlin.de, grassmann@tu-berlin.de',
    url='https://gitlab.klima.tu-berlin.de/klima/uc2data.git',
    license=license,
    packages=find_packages(exclude=('tests', 'docs'))
)
