# User interfaces

$\pi$-PIC containers can be created by calling (whenever applicable default values of parameters are specified after `=`):

- `init(solver, nx, XMin, XMax, ny=1, YMin=-0.5, YMax=0.5, nz=1, ZMin=-0.5, ZMax=0.5)` allocates data for the field grid and creates a container for particles
    - `nx, ny, nz ` are the grid sizes along $x$, $y$ and $z$; must be powers of two (by default `ny=nz=1,`, i.e. 1D simulation, for 2D case keep `nz=1`). 
    - `XMin, XMax, YMin, YMax, ZMin, ZMax` are the physical limits of the simulations region (by default `YMax=ZMax=-YMin=-ZMin=0.5`)
    - `solver` provides a way to select one of the following solvers to be used for the simulation:
        - `fourier_boris` is a combination of standard Boris pusher and Fourier-based field solver
        - `ec` is a first-order energy-conserving solver with Fourier-based field solver
        - `ec2` is a second-order energy-conserving solver with Fourier-based field solver

Interfaces and functions of containers 
--

- `addParticles(name, charge, mass, number, temperature, density, dataDouble = 0, dataInt = 0)` adds particles of new or existing type (in case of an exact match)
    - `name`, `charge` and `mass` are attributes of the type of particles to be added
    - `number` is the number of macroparticles to be distributed, the weight is determined automatically
    - `density` is the address of a callback that gives density as a function of coordinates: `density_callback(r, dataDouble, dataInt)` (decorator `@cfunc(type_addParticles)` must be placed before function definition)
        - `r[0]`, `r[1]` and `r[2]` are $x$, $y$ and $z$ coordinates (read only)
        - `dataDouble` and `dataInt` are pointers to data of double and int type (read/write) 
    - `dataDouble` and `dataInt` are addresses for exchanging data between `density` function and the remaining Python script

- `particleLoop(name, handler, dataDouble, dataInt)` makes a loop over particles of type `name` and applies action `handler` (non-OMP loop, fast OMP loops can be arranged via extension)
    - `name` defines the type of particles to be processed
    - `handler` is the address of a callback that defines actions to be performed with each particle: `handler(r, p, w, id, dataDouble, dataInt)` (decorator `@cfunc(type_particleLoop)` must be placed before function definition)
        - `r[0]`, `r[1]` and `r[2]` are $x$, $y$ and $z$ coordinates of the particle being processed (read only)
        - `p[0]`, `p[1]` and `p[2]` are momentum components along $x$, $y$ and $z$ (read/write)
        - `w[0]` is the particle's weight, can be set to zero to remove the particle (read/write)
        - `id[0]` is the particle's id in double representation (read only)
        - `dataDouble` and `dataInt` are pointers to data of double and int type (read/write)
    - `dataDouble` and `dataInt` are addresses for exchanging data between `handler` function and the remaining Python script

- `fieldLoop(handler, dataDouble = 0, dataInt = 0, useOmp = False)` makes a loop over grid nodes
    - `handler` is the address of a callback that defines the action on field components: `handler(ind, r, E, B, dataDouble, dataInt)` (decorator `@cfunc(type_fieldLoop)` must be placed before function definition)
        - `ind[0]`, `ind[1]` and `ind[2]` are indices of the gird node along $x$, $y$ and $z$, respectively (read only)
        - `r[0]`, `r[1]` and `r[2]` are the coordinates along $x$, $y$ and $z$, respectively (read only, collocated grid is used in all solvers)
        - `E[0]`, `E[1]`, `E[2]`, `B[0]`, `B[1]` and `B[2]` are the components of electric and magnetic fields (read/write)
        - `dataDouble` and `dataInt` are pointers to data of double and int type (read/write)
    - `dataDouble` and `dataInt` are addresses for exchanging data between `handler` function and the remaining Python script
    - `useOmp` defines whether to make parallel loop (via openMP), if `useOmp = True` the `handler` must be thread-safe

- `customFieldLoop(numberOfIterations, it2r, field2data, dataDouble = 0, dataInt = 0)` makes a loop over a subset of locations, at which electromagnetic field is interpolated and sent to callback `field2data` for output
    - `numberOfIterations` defines the number of calls, i.e. the size of the subset
    - `it2r` is the address of a callback that defines the locations, at which the field interpolated is needed: `handler(it, r, dataDouble, dataInt)` (decorator `@cfunc(type_it2r)` must be placed before function definition)
        - `it[0]` is the index of the subset element (read only)
        - `r[0]`, `r[1]` and `r[2]` are the coordinates along $x$, $y$ and $z$, respectively (write)
        - `dataDouble` and `dataInt` are pointers to data of double and int type (read/write)
    - `field2data` is the address of a callback that defines how the interpolated field contributes to data output via `dataDouble` and `dataInt`: `handler(it, r, E, B, dataDouble, dataInt)` (decorator `@cfunc(type_field2data)` must be placed before function definition)
        - `it[0]` is the index of the subset element (read only)
        - `r[0]`, `r[1]` and `r[2]` are the coordinates along $x$, $y$ and $z$, respectively (read only)
        - `E[0]`, `E[1]`, `E[2]`, `B[0]`, `B[1]` and `B[2]` are the components of interpolated electric and magnetic fields (read only)
        - `dataDouble` and `dataInt` are pointers to data of double and int type (read/write)
    - `dataDouble` and `dataInt` are addresses for exchanging data between `handler` function and the remaining Python script

- `advance(timeStep, numberOfIterations = 1)` advances the state of field and particles
    - `timeStep` is time step
    - `numberOfIterations` is the number of iterations to perform

- `fourierSolverSettings(divergenceCleaning=-1, sin2_kFilter=-1)` provides additional setting of Fourier-based field solver, value `-1` means keeping unchanged, `0` and `1` mean disable and enable, respectively
    - `divergenceCleaning` defines is divergence cleaning is applied (`True` by default for `solver=fourier_boris`)
    - `sin2_kFilter` defines if high frequency field components are suppressed (see [fourier_solver.h](src/) lines 219-224)
    
- `setRngGenSeed(seed)` provides a way to change the seed used for pseudo-random number generation

- `logPolicy(logToFile=True, logToScreen=False)` defines how log messages are reported

- `getTypeIndex(typeName)` returns the integer index of the type with name `typeName`

- `addHandler(name, subject, handler, dataDouble = 0, dataInt)` incorporates an extension, which can add and remove particles during state advancement (see [Making extensions](making_extensions.md))
    - `name` is the name of the extension
    - `subject` a string list of particle types to be affected, `cell` indicates that the handler is to be called for all cells (even if empty), `cell` can be used to affect all particles
    - `handler` is the address of the function to act on particles (see [Making extensions](making_extensions.md))
    - `dataDouble` and `dataInt` are addresses for exchanging data between `handler` function and the remaining Python script

Convensions, assumptions and properties
--

- CGS units are used for all dimensional quantities
- The components of vectors are enumerated: $x$, $y$, $z$ correspond to `[0]`, `[1]`, `[2]`, respectively
- The sizes of the grid `nx`, `ny` and `nz` must be powers of two (to facilitate FFT)
- 2D simulation is enabled by setting `nz=1`, 1D simulations is enabled by `nz=1` and `ny=1`
- The topology of space is toroidal (periodic boundary conditions)
- Before adding particles the field can be advanced over an arbitrary time, however particles should not traverse a distance larger then one spatial step over a single time step
- All solvers are deterministic (results are exactly reproducible)
