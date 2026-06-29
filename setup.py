from setuptools import setup, find_packages
import subprocess
import sys
import os
import shutil
import urllib.request
import numpy as np

def get_numpy_i():
    """Download numpy.i from GitHub matching the installed numpy version."""
    dst = os.path.join('code', 'helper', 'numpy.i')
    version = np.__version__

    # try exact version tag first
    url = f"https://raw.githubusercontent.com/numpy/numpy/v{version}/tools/swig/numpy.i"
    try:
        print(f"Downloading numpy.i for numpy {version}...")
        urllib.request.urlretrieve(url, dst)
        print(f"Downloaded numpy.i from {url}")
        return
    except urllib.error.HTTPError:
        print(f"Exact version v{version} not found, falling back to main branch...")

    # fallback to main branch
    url = "https://raw.githubusercontent.com/numpy/numpy/main/tools/swig/numpy.i"
    try:
        urllib.request.urlretrieve(url, dst)
        print("Downloaded numpy.i from main branch")
        return
    except urllib.error.HTTPError:
        raise FileNotFoundError(
            "Could not download numpy.i from GitHub. "
            "Check your internet connection or manually place numpy.i in code/helper/"
        )

# download numpy.i if not present or if numpy version changed
get_numpy_i()

# run make
subprocess.run(["make"], stdout=sys.stdout, stderr=sys.stderr, check=True)

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
