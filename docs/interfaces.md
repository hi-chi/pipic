# User interfaces

Interfaces and functions of $\pi$-PIC containers 
--
(whenever applicable default values of parameters are specified after `=`)

- `init(solver, nx, XMin, XMax, ny=1, YMin=-0.5, YMax=0.5, nz=1, ZMin=-0.5, ZMax=0.5)` allocates data for the field grid and creates a container for particles
    - `nx, ny, nz ` are the grid sizes along $x$, $y$ and $z$; must be powers of two (by default `ny=nz=1,`, i.e. 1D simulation, for 2D case keep `nz=1`). 
    - `XMin, XMax, YMin, YMax, ZMin, ZMax` are the physical limits of the simulations region (by default `YMax=ZMax=-YMin=-ZMin=0.5`)
    - `solver` provides a way to select one of the following solvers to be used for the simulation:
        - `fourier_boris` is a combination of standard Boris pusher and Fourier-based field solver
        - `ec` is a first-order energy-conserving solver with Fourier-based field solver
        - `ec2` is a second-order energy-conserving solver with Fourier-based field solver

- `addParticles(name, number, charge, mass, temperature, density, dataDouble = 0, dataInt = 0)` adds particles of new or existing type (in case of an exact match)

- `fieldLoop(handler, dataDouble = 0, dataInt = 0, useOmp = False)` makes a loop over grid nodes

- `customFieldLoop(numberOfIterations, it2r, field2data, dataDouble = 0, dataInt = 0)` makes a loop over a subset of location, at which electromagnetic field is interpolated and sent to callback `field2data` for output

- `advance(timeStep, numberOfIterations = 1)` advances the state of field and particles

- `fourierSolverSettings(divergenceCleaning=-1, sin2_kFilter=-1)` provides additional setting of Fourier-based field solver, value `-1` means keeping unchanged, `0` and `1` mean disable and enable, respectively

- `setRngGenSeed(seed)` provides a way to change the seed used for pseudo-random number generation

- `logPolicy(logToFile=True, logToScreen=False)` defines how log messages are reported

- `getTypeIndex(typeName)` returns the iteger index of the type with name `typeName`

- `addHandler(name, subject, handler, dataDouble = 0, dataInt)` incorporates an extension, which can add and remove particles during state advancement (see [Making extensions](making_extensions.md))

Convensions, assumptions and properties
--

- CGS units are used for all dimensional quantities
- The components of vectors are enumerated: $x$, $y$, $z$ correspond to `[0]`, `[1]`, `[2]`, respectively
- The sizes of the grid `nx`, `ny` and `nz` must be powers of two (to facilitate FFT)
- 2D simulation is enabled by setting `nz=1`, 1D simulations is enabled by `nz=1` and `ny=1`
- The topology of space is toroidal (periodic boundary conditions)
- Before adding particles the field can be advanced over an arbitrary time, however particles should not traverse a distance larger then one spatial step over a single time step
- All solvers are deterministic (results are exactly reproducible)