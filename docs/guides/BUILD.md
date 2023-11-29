# Build guide
The following document provides a guide for building $pi$-PIC in order to make it available through [PyPI](https://pypi.org/).

### Project structure
When building the project, the C-module is compiled into a binary `_pipic*.so` (linkable object). This module can be imported in a python script by writing e.g. `import _pipic`. The source code for the C-module is located under [src/](../../src) and the python-package is found separately under [pipic/](../../pipic), wherein the `_pipic` module is imported (this remains hidden for the end-user). Because python generally treat files, packages and variables with a leading underscore as "private", extra care is needed when importing such modules, and this also presents some complications with certain build-scripts.

### Stubs
Because `_pipic` is a binary package, IDEs are unable to provide tab completion, type hinting and other information except at runtime, _e.g._ in a jupyter notebook. To fix this, `stubs` are needed to expose the interface. This is done with [\_\_init\_\_.pyi](../../pipic/__init__.pyi), which is a stub file (`.pyi`), that can be automatically generated using [mypy](https://mypy-lang.org/) (`pip install mypy`). However, when providing a stub file it entirely replaces that which an IDE could otherwise infer (_e.g._ from the python-package). To address this we need combine the stubs for `_pipic` and `pipic` when creating [\_\_init\_\_.pyi](../../pipic/__init__.pyi).

To generate the stub file, simply write the following (assuming the `_pipic*.so` binary is available):
```
stubgen -m _pipic
```

## Build
### Overview
The build process primarily revolves around three files: [setup.py](../../setup.py), [pyproject.toml](../../pyproject.toml) and [MANIFEST.in](../../MANIFEST.in), where we use [`setuptools`](https://pypi.org/project/setuptools/) to compile and install the package. This can be done by writing:
```
python setup.py install
```
which generates several folders and installs the package in your environment in the form of an `*.egg` (generates local folders: `build/`, `dist/`, `pipic.egg-info/`, `var/`). To instead install the package using _pip_, one can write:
```
pip install .
```
This also installs the package in your local environment, but now generates no excess local folders. The installation is now split between the compiled binary on one hand and the python package on the other.

A quick summary of the three files mentioned above:
- [setup.py](../../setup.py) is the primary script used to actually build the project, using `setuptools`. It defines the build parameters, such as flags and installation dependencies, and grabs project metadata from [\_\_init\_\_.py](../../pipic/__init__.py).

- [pyproject.toml](../../pyproject.toml) is used to define the build dependencies. This is needed because _e.g._ pybind11 is required for building before `setup()` is run in [setup.py](../../setup.py).

- [MANIFEST.in](../../MANIFEST.in) is used to define which files and folders will be used in the build process. When building, these files are moved to a separate folder. Failure to specify them leads to errors, _e.g._ due to missing `.h` files.

### Deployment
To build the project for deployment requires the python package `build` (`pip install build`). To build, simply write:
```
CC=gcc-13 CXX=g++-13 python -m build
```
This creates two local folders, `dist/` and `pipic.egg-info/`, where the former contains the distribution file(s), often in terms of a binary wheel (`.wh`) and source archive (`.tar.gz`). Distribution can be done with either (or both) but source is generally recommended as guaranteeing compatibility of binaries is likely to require more work.

Finally, distribution of the package is done using `twine` (`pip install twine`). To test the distribution you can write:
```
twine check dist/*
```

In order to actually distribute it, here via [TestPyPI](https://test.pypi.org/), you write (using _only_ source):
> [!CAUTION]  
> Only run the following command if you actually want to distribute the current build.
```
twine upload --repository testpypi dist/*.tar.gz
```

## Links
For more information see the following PEPs:
- [PEP440](https://peps.python.org/pep-0440/)
- [PEP518](https://peps.python.org/pep-0518/)
- [PEP621](https://peps.python.org/pep-0621/)

For guides, see:
- https://packaging.python.org/en/latest/tutorials/packaging-projects/
- https://scikit-hep.org/developer

On version capping, read:
- https://iscinumpy.dev/post/bound-version-constraints/