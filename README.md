# plato
Python LAMMPS Analysis Tools

## Install this package

Run the following in the package root directory
```
pip install .
```

During the installation, it will run the Makefile by `make`. To see the output of `make`, use the `-v` option of pip:
```
pip install . -v
```

**Note:** The `numpy.i` file in `code/helper` might not be compatible with SWIG versions newer than 4.0. If you run into SWIG-related installation issues, try downloading a compatible or updated version of `numpy.i` and replacing the existing file in `code/helper`.
