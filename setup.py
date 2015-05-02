#!/usr/bin/env python
from setuptools import setup, find_packages
import itertools

options = dict(
    name='gamma_index',
    version='0.1',
    packages=find_packages(),
    license='MIT',
    include_package_data = True,
    description='gamma_index - calculation of gamma index on multi-dimensional distributions',
    long_description=open('README.rst').read(),
    author='Jan Pipek',
    author_email='jan.pipek@gmail.com',
    url='https://github.com/janpipek/gamma_index',
    install_requires = ['numpy']
)
setup(**options)