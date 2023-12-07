<p align="center">
<img src="https://raw.githubusercontent.com/hi-chi/pipic/main/docs/logo/pipic_logo.png" width="300">
</p>

---

[![PyPI version][pypi-version]][pypi-link]
[![PyPI platforms][pypi-platforms]][pypi-link]

$\pi$-PIC (PIPIC, Python-controlled Interactive PIC) is an open-source collection of relativistic particle-in-cell solvers featuring
- exact energy conservation;
- absence of numerical dispersion.

The solvers provide a way to either suppress or eliminate numerical artefacts (instabilities, heating, numerical Cherenkov radiation, etc.) permitting larger space and time steps, as well as lower number of particles per cell.
Because of reduced computational demands, the solvers can be found useful for quick tests of ideas, as well as for scanning parameter spaces. For a description of the underlying methods see [References](#References).

---

# Overview
$\pi$-PIC provides all tools necessary for designing 1D/2D/3D simulations and arbitrary outputs directly from Python. In addition, it has interfaces for incorporating extensions (read/modify field and particles, add/remove particles) that can be developed in Python, C/C++, Fortran or any other language that generate callable functions. The project and its development are hosted on [GitHub][]. 

To get started, it should for most cases be sufficient to install $\pi$-PIC via _pip_ (this requires: `gcc`, `openmp` and `fftw3`, see [installation instructions][installation] for details):
```
pip install pipic
```

The basic layout of a simulation includes five elements:
- creating a container with cells with given parameters
- adding particles of all necessary types
- setting initial electromagnetic field
- defining output (via loops over particles and grid values of field)
- advance and read the state of the defined physical system

We demonstrate the use of these elements in the [tutorial][]. A complete list of supported interfaces can be found [here][interfaces]. The development of extensions is detailed and exemplified [here][extensions]. 


# References
A. Gonoskov, Explicit energy-conserving modification of relativistic PIC method, [arXiv:2302.01893][] (2023).


<!-- prettier-ignore-start -->
[pypi-link]:                https://pypi.org/project/pipic/
[pypi-platforms]:           https://img.shields.io/pypi/pyversions/pipic
[pypi-version]:             https://badge.fury.io/py/pipic.svg

[GitHub]: https://github.com/hi-chi/pipic
[installation]: https://github.com/hi-chi/pipic/blob/main/docs/guides/INSTALLATION.md
[tutorial]: https://github.com/hi-chi/pipic/blob/main/docs/guides/TUTORIAL.md
[interfaces]: https://github.com/hi-chi/pipic/blob/main/docs/guides/INTERFACES.md
[extensions]: https://github.com/hi-chi/pipic/blob/main/docs/guides/EXTENSIONS.md
[arXiv:2302.01893]: https://arxiv.org/abs/2302.01893
<!-- prettier-ignore-end -->