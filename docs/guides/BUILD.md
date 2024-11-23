# Build guide
The following document provides a guide for building $\pi$-PIC in order to make it available through [PyPI](https://pypi.org/).

### Project structure
When building the project, the C-module is compiled into a binary `_pipic*.so` (linkable object). This module can be imported in a python script by writing e.g. `import _pipic`. The source code for the C-module is located under [src/](../../src) and the python-package is found separately under [pipic/](../../pipic), wherein the binary `_pipic` module is imported (this remains hidden for the end-user). Because python generally treat files, packages and variables with a leading underscore as "private", extra care is needed when importing such modules, and this also presents some complications with certain build-scripts.

### Stubs
Because `_pipic` is a binary package, IDEs are unable to provide tab completion, type hinting and other information except at runtime, _e.g._ in a jupyter notebook. To fix this, `stubs` are needed to expose the interface. This is done with [\_\_init\_\_.pyi](../../pipic/__init__.pyi), which is a stub file (`.pyi`), that can be automatically generated using [mypy](https://mypy-lang.org/) (`pip install mypy`). However, when providing a stub file it entirely replaces that which an IDE could otherwise infer (_e.g._ from the python-package). To address this we need to combine the stubs for `_pipic` and `pipic` when creating [\_\_init\_\_.pyi](../../pipic/__init__.pyi).

To generate the stub file, simply write the following (assuming the `_pipic*.so` binary is available):
```
stubgen -m _pipic
```

## Build
### Overview
The build process primarily revolves around two files: [setup.py](../../setup.py) and [pyproject.toml](../../pyproject.toml) (previously also a [MANIFEST.in](../../MANIFEST.in)), where we use [`setuptools`](https://pypi.org/project/setuptools/) to compile and install the package. This can be done by writing:
```
python setup.py install
```
which generates several folders and installs the package in your environment in the form of an `*.egg` (generates local folders: `build/`, `dist/`, `pipic.egg-info/`, `var/`). To instead install the package using _pip_, one can write:
```
pip install .
```
This also installs the package in your local environment (for example `somewhere/site-packages/`), but now generates no excess local folders. The installation is now split between the compiled binary on one hand and the python package on the other.

A quick summary of the three files mentioned above:
- [setup.py](../../setup.py) is the primary script used to actually build the project, using `setuptools`. It defines the build parameters, such as flags and installation dependencies. (Previously grabbed project metadata from [\_\_init\_\_.py](../../pipic/__init__.py)).

- [pyproject.toml](../../pyproject.toml) is used to define the build dependencies. This is needed because _e.g._ pybind11 is required for building before `setup()` is run in [setup.py](../../setup.py). Project metadata is defined statically in this file.

- [MANIFEST.in](../../MANIFEST.in) can be used to define which files and folders will be used in the build process. When building, these files are moved to a separate folder. Failure to specify them may lead to errors, _e.g._ due to missing `.h` files. This file can also be used to exclude files from the source distribution (`sdist`), such as excluding the `examples/`, `docs/` or `tests/` folders. To do so, one can for example add the following to the [MANIFEST.in](../../MANIFEST.in):
```
prune docs
prune examples
prune tests
exclude MANIFEST.in .*ignore
```

### Deployment
The publishing process can be automated using [cibuildwheel](https://cibuildwheel.pypa.io), in order to build wheels for multiple platforms, architectures and python versions automatically. This is done in the `deploy.yml` workflow, which is automatically triggered when a new version tag ([see below](#Versioning)) is pushed to the main branch.

#### Manually
To build the project for manual deployment requires the python package `build` (`pip install build`). To build, simply write:
```
python -m build
```
This creates two local folders, `dist/` and `pipic.egg-info/`, where the former contains the distribution file(s), often in terms of a binary wheel (`.wh`) and source distribution (`.tar.gz`). Distribution can be done with either (or both) but source is generally recommended as guaranteeing compatibility of binaries is likely to require more work.

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

In order to distribute it via [PyPI](https://pypi.org/), you instead write
> [!CAUTION]
> Only run the following command if you actually want to distribute the current build.
> You may have several `.tar.gz` in your `dist/` folder. Make sure to only distribute the one you want.
```
twine upload dist/*.tar.gz
```

### Versioning
Versioning is done automatically using [`setuptools-scm`](https://setuptools-scm.readthedocs.io/en/latest/config/#api-reference). When building the project, a version file `_version.py` containing the version number is created, based on the latest `git` tag. This file should not be tracked by git, but must be shipped with source distribution (this is already done automatically). All versioning should rely on this file.

To check the presumptive version, simply write
```bash
python setup.py --version
```

To bump the version, create a `git` tag beginning with `v` and following the [_SemVer_](https://semver.org/) scheme. Although not enforced, it is strongly encouraged to use an annotated tag (with the `-a` flag) as this creates a record of when the tag was created, by whom, and allow for attaching a comment.

For example:
```bash
git tag -a v1.3 -m "pi-PIC v1.3"
```
will create the tag `v1.3`, which can then be pushed to the remote repo using:
```bash
git push origin --follow-tags
```
or
```bash
git push origin v1.3
```


## Links
For more information see the following PEPs:
- [PEP440](https://peps.python.org/pep-0440/)
- [PEP518](https://peps.python.org/pep-0518/)
- [PEP621](https://peps.python.org/pep-0621/)

For guides, see:
- https://packaging.python.org/en/latest/tutorials/packaging-projects/
- https://scikit-hep.org/developer
- https://pypi.org/project/setuptools-scm/4.1.2/
- https://learn.scientific-python.org/development/guides/packaging-compiled/
- https://setuptools.pypa.io/en/latest/userguide/miscellaneous.html

On version capping, read:
- https://iscinumpy.dev/post/bound-version-constraints/
