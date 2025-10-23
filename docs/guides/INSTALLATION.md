# Installation

$\pi$-PIC has been tested and confirmed to work Linux / macOS / Windows under [`WSL`](learn.microsoft.com/windows/wsl/about). Installation requirements are:
- Python 3
- [`gcc`](https://gcc.gnu.org/)
- [`fftw3`](http://www.fftw.org/) (`gcc` and `fftw3` must be both of 64-bit format)
- [`openmp`](https://www.openmp.org/)
- [`pybind11`](https://github.com/pybind/pybind11) (handled automatically when installed through `pip`)

Furthermore, the following python packages are generally needed for running the code:
- [`numba`](https://numba.pydata.org/)
- [`numpy`](https://numpy.org/)



## PyPI
$\pi$-PIC is available through PyPI, which is the recommended installation method:
```
pip install pipic
```
If using a non-standard compiler (for example using `gcc` installed through [`homebrew`](https://brew.sh/) on macOS) this can be specified as:
```
CC=<c-compiler> CXX=<c++-compiler> pip install pipic
```

## Conda
Follow the instructions at: https://github.com/hi-chi/pipic/blob/pipic_tutorial/tutorials/installation.ipynb
```


### macOS
Installation on macOS depends on which compiler is being used.

#### gcc
Assuming that `gcc` is installed through homebrew it will be named _e.g._ `gcc-14`. Note that the built-in `gcc` on macOS is in fact an Apple-compiled `clang` (this can be confirmed with the `gcc --version`), which lacks `openmp`.

First, make sure you have `fftw` installed
```
brew install fftw
```
Then, assuming `gcc-14`, `pipic` can be installed with
```
export LDFLAGS="-L$(brew --prefix)/opt/fftw/lib"
export CPPFLAGS="-I$(brew --prefix)/opt/fftw/include"
CC=gcc-14 CXX=g++-14 pip install pipic
```

The export flag statements can be added to `.bashrc`/`.bash_profile` in order to not have to manually specify it in every shell.

#### clang
`openmp` can be installed for `clang` through homebrew, using `libomp`.
```
brew install fftw libomp
```
Then, assuming that the default compiler (`gcc`) has not been changed (_i.e._ it still invokes `clang`), `pipic` can be installed with
```
export LDFLAGS="-L$(brew --prefix)/opt/fftw/lib -L$(brew --prefix)/opt/libomp/lib"
export CPPFLAGS="-I$(brew --prefix)/opt/fftw/include -I$(brew --prefix)/opt/libomp/include"
pip install pipic
```


## Compiling latest or specific version from GitHub

1. Create a local folder:
	`mkdir <name>`
	`cd <name>`
2. Clone the project from there and cd to a folder named `pipic`:
	`git clone https://github.com/hi-chi/pipic.git`
	`cd pipic`
3. (optional) Switch to a needed branch:
	`git checkout <name>`
4. (optional) Uninstall previous version of pipic:
	`pip uninstall pipic`
5. Install pipic (you might need to use python instead of python3):
	`python3 -m pip install .`
6. (optional) Check that it works:
    `python3 examples/basic_example.py`

Note that the install command also here can be prefixed with `CC=<c-compiler> CXX=<c++-compiler>` in order to use a non-default compiler, and that `LDFLAGS`/`CPPFLAGS` may need to be supplied as shown earlier.

<!---
## CMake
As an alternative, manual compilation is also available using `cmake`. This extends the above requirements to also include:
- [`CMake`](https://cmake.org/) 3.1 or higher

To compile:
- Clone the repository
    ```
    git clone https://github.com/hi-chi/pipic.git
    ```
- Go to `src` folder and fetch `pybind11`:
    ```
    cd pipic/src
    git clone https://github.com/pybind/pybind11
    ```
- Generate binary .so-file by running sequentially:
    ```
    cmake .
    make
    ```
- To use compilers other than default, set the `CC` and `CXX` environment variables prior to running `cmake`, or pass the compilers as arguments to `cmake` using the `-DCMAKE_C_COMPILER=` and `-DCMAKE_CXX_COMPILER=` flags.

To use $\pi$-PIC requires that both the `pipic/` subfolder (python package) and the `_pipic.*.so` binary (C/C++ package) are present in your project folder or otherwise made available, e.g. through your `$PATH`.
-->
