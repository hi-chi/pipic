<p align="center">
<img src="https://raw.githubusercontent.com/hi-chi/pipic/testpypi/docs/logo/pipic_logo.png" width="300">
</p>

$\pi$-PIC (PIPIC, Python-controlled Interactive PIC) is an open-source collection of relativistic particle-in-cell solvers featuring
- exact energy conservation;
- absence of numerical dispersion.

The solvers provide a way to either suppress or eliminate numerical artefacts (instabilities, heating, numerical Cherenkov radiation, etc.) permitting larger space and time steps, as well as lower number of particles per cell.
Because of reduced computational demands, the solvers can be found useful for quick tests of ideas, as well as for scanning parameter spaces. For a description of the underlying methods see [References](#References).

---

# Overview
$\pi$-PIC provides all tools necessary for designing simulations and arbitrary outputs directly from Python. In addition, it has interfaces for incorporating extensions (read/modify field and particles, add/remove particles) that can be developed in Python, C/C++, Fortran or any other language that generate callable functions. The project and its development is hosted on [GitHub](https://github.com/hi-chi/pipic). 

To get started, it should for most cases be sufficient to install $\pi$-PIC via _pip_ (this requires: `gcc`, `openmp` and `fftw3`, see [installation instructions](https://github.com/hi-chi/pipic/blob/testpypi/INSTALLATION.md) for details):
```
pip install pipic
```

The basic layout of a simulation includes five elements:
- creating a container with cells with given parameters
- adding particles of all necessary types
- setting initial electromagnetic field
- defining output (via loops over particles and grid values of field)
- advance and read the state of the defined physical system

We demonstrate the use of these elements in the [tutorial](https://github.com/hi-chi/pipic/blob/testpypi/TUTORIAL.md). A complete list of supported interfaces can be found [here](https://github.com/hi-chi/pipic/blob/testpypi/docs/interfaces.md). The development of extensions is detailed and exemplified [here](https://github.com/hi-chi/pipic/blob/testpypi/docs/making_extentions.md). 


# References
A. Gonoskov, Explicit energy-conserving modification of relativistic PIC method, [arXiv:2302.01893](https://arxiv.org/abs/2302.01893) (2023).
