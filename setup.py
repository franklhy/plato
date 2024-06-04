from setuptools import setup, find_packages
import subprocess

subprocess.run(["make"])

setup(
    name='plato',
    version='1.0.0',
    packages=['plato'],
    package_data={
        'plato': ['_plato.so', 'plato.py', '__init__.py'],
    },
    install_requires=[
        #'swig',
    ],
)
