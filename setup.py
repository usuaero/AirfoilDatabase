"""AirfoilDatabase: A python module for creating databases of complex airfoils."""

from setuptools import setup
import os

setup(name = 'airfoil_db',
    version = '1.0',
    description = "AirfoilDatabase: A Python module for creating databases of complex airfoils",
    url = 'https://github.com/usuaero/AirfoilDatabase',
    author = 'usuaero',
    author_email = 'doug.hunsaker@usu.edu',
    install_requires = ['numpy', 'scipy', 'matplotlib'],
    python_requires ='>=3.6.0',
    license = 'MIT',
    packages = ['airfoil_db'],
    zip_safe = False)
