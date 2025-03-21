<p align="center">
<img src="https://raw.githubusercontent.com/hi-chi/pipic/main/docs/logo/pipic_logo.png" width="300">
</p>

---
[![Actions Status][actions-badge]][actions-link]
[![PyPI version][pypi-version]][pypi-link]
[![PyPI platforms][pypi-platforms]][pypi-link]

$\pi$-PIC (PIPIC, Python-controlled Interactive PIC) is an open-source collection of relativistic particle-in-cell solvers featuring
- exact energy conservation;
- absence of numerical dispersion.

The solvers provide a way to either suppress or eliminate numerical artefacts (instabilities, heating, numerical Cherenkov radiation, etc.) permitting larger space and time steps, as well as lower number of particles per cell.
Because of reduced computational demands, the solvers can be found useful for quick tests of ideas, as well as for scanning parameter spaces. For a description of the underlying methods see [Reference](#Reference) or a [presentation at PIF24](https://plymouth.cloud.panopto.eu/Panopto/Pages/Viewer.aspx?id=86a17791-6db2-4009-be2f-b1dc00c846d0).

---

# Overview
$\pi$-PIC provides all tools necessary for designing 1D/2D/3D simulations and arbitrary outputs directly from Python. In addition, it has interfaces for incorporating extensions (read/modify field and particles, add/remove particles) that can be developed in Python, C/C++, Fortran or any other language that generate callable functions (see [extensions][extensions] for a list of extensions included in $\pi$-PIC installation). The project and its development are hosted on [GitHub][].

To get started, it should for most cases be sufficient to install $\pi$-PIC via _pip_ (this requires: `gcc`, `openmp` and `fftw3`; for details and information on compilation via CMake see [installation instructions][installation]):
```
pip install pipic
```

The basic layout of a simulation includes five elements:
- creating a container with cells with given parameters
- adding particles of all necessary types
- setting initial electromagnetic field
- defining output (via loops over particles and grid values of field)
- advance and read the state of the defined physical system

We demonstrate the use of these elements in the [tutorial][]. A complete list of supported interfaces can be found [here][interfaces]. The development of extensions is detailed and exemplified [here][extension_development].

# New in $\pi$-PIC v1.2
- options for energy correction routine in `ec` and `ec2` (see `docs/guides/INTERFACES.md`)
- [`extension`][extensions] for initializing arbitrary tightly focused pulses, e.g. dipole waves (`focused_pulse`)

# New in $\pi$-PIC v1.1
- [`extensions`][extensions] for QED-PIC simulations (`qed_volokitin2023`, `qed_gonoskov2015`)
- [`extension`][extensions] for ensemble down-sampling (`downsampler_gonoskov2022`)

See all releases [here][releases].

# Reference
A. Gonoskov, Explicit energy-conserving modification of relativistic PIC method, [J. Comput. Phys., 502, 112820](https://doi.org/10.1016/j.jcp.2024.112820); [arXiv:2302.01893][] (2024).



<!-- prettier-ignore-start -->
[actions-badge]:            https://github.com/hi-chi/pipic/workflows/CI/badge.svg
[actions-link]:             https://github.com/hi-chi/pipic/actions
[pypi-link]:                https://pypi.org/project/pipic/
[pypi-platforms]:           https://img.shields.io/pypi/pyversions/pipic
[pypi-version]:             https://badge.fury.io/py/pipic.svg

[GitHub]: https://github.com/hi-chi/pipic
[installation]: https://github.com/hi-chi/pipic/blob/main/docs/guides/INSTALLATION.md
[tutorial]: https://github.com/hi-chi/pipic/blob/main/docs/guides/TUTORIAL.md
[interfaces]: https://github.com/hi-chi/pipic/blob/main/docs/guides/INTERFACES.md
[extensions]: https://github.com/hi-chi/pipic/blob/main/docs/EXTENSIONS.md
[extension_development]: https://github.com/hi-chi/pipic/blob/main/docs/guides/EXTENSION_DEVELOPMENT.md
[releases]: https://github.com/hi-chi/pipic/blob/main/docs/RELEASES.md
[arXiv:2302.01893]: https://arxiv.org/abs/2302.01893
<!-- prettier-ignore-end -->
