# code/helper/

This directory contains helper files for the SWIG interface.

## numpy.i

`numpy.i` is **not tracked by git** — it is automatically downloaded from the
numpy GitHub repository at build time (`pip install -e .`).

The correct version is fetched to match the installed numpy version.
To manually download it:

```python
import urllib.request, numpy as np
url = f"https://raw.githubusercontent.com/numpy/numpy/v{np.__version__}/tools/swig/numpy.i"
urllib.request.urlretrieve(url, "code/helper/numpy.i")
```
