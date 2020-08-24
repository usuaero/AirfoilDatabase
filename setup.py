"""AirfoilDatabase: A Python module for modeling airfoils using Xfoil."""

from setuptools import setup
import os

setup(name = 'airfoil_db',
    version = '1.4.2',
    description = "AirfoilDatabase: A Python module for modeling airfoils using Xfoil.",
    url = 'https://github.com/usuaero/AirfoilDatabase',
    author = 'usuaero',
    author_email = 'doug.hunsaker@usu.edu',
    install_requires = ['numpy', 'scipy', 'matplotlib'],
    python_requires ='>=3.6.0',
    license = 'MIT',
    packages = ['airfoil_db'],
    zip_safe = False)
