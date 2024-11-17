# User interfaces

$\pi$-PIC containers can be created by calling (whenever applicable default values of parameters are specified after `=`):

- `init(solver, nx, xmin, xmax, ny=1, ymin=-0.5, ymax=0.5, nz=1, zmin=-0.5, zmax=0.5)` allocates data for the field grid and creates a container for particles
    - `nx, ny, nz ` are the grid sizes along $x$, $y$ and $z$; must be powers of two (by default `ny=nz=1,`, i.e. 1D simulation, for 2D case keep `nz=1`).
    - `xmin, xmax, ymin, ymax, zmin, zmax` are the physical limits of the simulations region (by default `ymax=zmax=-ymin=-zmin=0.5`)
    - `solver` provides a way to select one of the following solvers to be used for the simulation:
        - `fourier_boris` is a combination of standard Boris pusher and Fourier-based field solver
        - `ec` is a first-order energy-conserving solver with Fourier-based field solver
        - `ec2` is a second-order energy-conserving solver with Fourier-based field solver

Interfaces and functions of containers
--

- `add_particles(name, charge, mass, number, temperature, density, data_double = 0, data_int = 0)` adds particles of new or existing type (in case of an exact match)
    - `name`, `charge` and `mass` are attributes of the type of particles to be added
    - `number` is the number of macroparticles to be distributed, the weight is determined automatically
    - `density` is the address of a callback that gives density as a function of coordinates: `density_callback(r, data_double, data_int)` (decorator `@cfunc(types.add_particles_callback)` must be placed before function definition)
        - `r[0]`, `r[1]` and `r[2]` are $x$, $y$ and $z$ coordinates (read only)
        - `data_double` and `data_int` are pointers to data of double and int type (read/write)
    - `data_double` and `data_int` are addresses for exchanging data between `density` function and the remaining Python script

- `particle_loop(name, handler, data_double, data_int)` makes a loop over particles of type `name` and applies action `handler` (non-OMP loop, fast OMP loops can be arranged via extension)
    - `name` defines the type of particles to be processed
    - `handler` is the address of a callback that defines actions to be performed with each particle: `handler(r, p, w, id, data_double, data_int)` (decorator `@cfunc(types.particle_loop_callback)` must be placed before function definition)
        - `r[0]`, `r[1]` and `r[2]` are $x$, $y$ and $z$ coordinates of the particle being processed (read only)
        - `p[0]`, `p[1]` and `p[2]` are momentum components along $x$, $y$ and $z$ (read/write)
        - `w[0]` is the particle's weight, can be set to zero to remove the particle (read/write)
        - `id[0]` is the particle's id in double representation (read only)
        - `data_double` and `data_int` are pointers to data of double and int type (read/write)
    - `data_double` and `data_int` are addresses for exchanging data between `handler` function and the remaining Python script

- `field_loop(handler, data_double = 0, data_int = 0, use_omp = False)` makes a loop over grid nodes
    - `handler` is the address of a callback that defines the action on field components: `handler(ind, r, E, B, data_double, data_int)` (decorator `@cfunc(types.field_loop_callback)` must be placed before function definition)
        - `ind[0]`, `ind[1]` and `ind[2]` are indices of the grid node along $x$, $y$ and $z$, respectively (read only)
        - `r[0]`, `r[1]` and `r[2]` are the coordinates along $x$, $y$ and $z$, respectively (read only, collocated grid is used in all solvers)
        - `E[0]`, `E[1]`, `E[2]`, `B[0]`, `B[1]` and `B[2]` are the components of electric and magnetic fields (read/write)
        - `data_double` and `data_int` are pointers to data of double and int type (read/write)
    - `data_double` and `data_int` are addresses for exchanging data between `handler` function and the remaining Python script
    - `use_omp` defines whether to make parallel loop (via openMP), if `use_omp = True` the `handler` must be thread-safe

- `custom_field_loop(number_of_iterations, it2r, field2data, data_double = 0, data_int = 0)` makes a loop over a subset of locations, at which electromagnetic field is interpolated and sent to callback `field2data` for output
    - `number_of_iterations` defines the number of calls, i.e. the size of the subset
    - `it2r` is the address of a callback that defines the locations, at which the field interpolated is needed: `handler(it, r, data_double, data_int)` (decorator `@cfunc(types.it2r_callback)` must be placed before function definition)
        - `it[0]` is the index of the subset element (read only)
        - `r[0]`, `r[1]` and `r[2]` are the coordinates along $x$, $y$ and $z$, respectively (write)
        - `data_double` and `data_int` are pointers to data of double and int type (read/write)
    - `field2data` is the address of a callback that defines how the interpolated field contributes to data output via `data_double` and `data_int`: `handler(it, r, E, B, data_double, data_int)` (decorator `@cfunc(types.field2data_callback)` must be placed before function definition)
        - `it[0]` is the index of the subset element (read only)
        - `r[0]`, `r[1]` and `r[2]` are the coordinates along $x$, $y$ and $z$, respectively (read only)
        - `E[0]`, `E[1]`, `E[2]`, `B[0]`, `B[1]` and `B[2]` are the components of interpolated electric and magnetic fields (read only)
        - `data_double` and `data_int` are pointers to data of double and int type (read/write)
    - `data_double` and `data_int` are addresses for exchanging data between `handler` function and the remaining Python script

- `advance(time_step, number_of_iterations = 1, use_omp = True)` advances the state of field and particles
    - `time_step` is time step
    - `number_of_iterations` is the number of iterations to perform
    - `use_omp` defines whether to perform parallel processing using threads (via openMP); `use_omp = False` should only be used for testing purposes

- `fourier_solver_settings(divergence_cleaning=-1, sin2_kfilter=-1)` provides additional setting of Fourier-based field solver, value `-1` means keeping unchanged, `0` and `1` mean disable and enable, respectively
    - `divergence_cleaning` defines is divergence cleaning is applied (`True` by default for `solver=fourier_boris`)
    - `sin2_kfilter` defines if high frequency field components are suppressed (see [fourier_solver.h](../../src/fourier_solver.h) lines 219-224)

- `set_rng_seed(seed)` provides a way to change the seed used for pseudo-random number generation

- `log_policy(log_to_file = True, log_to_screen = False)` defines how log messages are reported

- `get_type_index(type_name)` returns the integer index of the type with name `type_name`

- `add_handler(name, subject, handler, data_double = 0, data_int)` incorporates an extension, which can add and remove particles during state advancement (see [Making extensions](EXTENSIONS.md))
    - `name` is the name of the extension
    - `subject` a string list of particle types to be affected, `cell` indicates that the handler is to be called for all cells (even if empty), `cell` can be used to affect all particles
    - `handler` is the address of the function to act on particles (see [Making extensions](EXTENSIONS.md))
    - `data_double` and `data_int` are addresses for exchanging data between `handler` function and the remaining Python script

- `ensemble_data()` returns a global pointer to data of all particles to be used in extensions if standard way of accessing particles is insufficient (see example of use in [downsampler_gonoskov2022](src/extensions/downsampler_gonoskov2022/downsampler_gonoskov2022.cpp))

- `en_corr_type(correction_type = 2)` can be used to change the way of correcting energy in `ec` and `ec2` solvers; the argument is an integer that enumerates the following options:
    - `0` no correction: the feedback to energy exchange is still accounted for (energy is preserved to the next order accuracy as compared to `boris` pusher) but the energy is not preserved to machine accuracy; the inaccuracy becomes larger if a particle goes from non-relativistic to relativistic motion (or vice versa) within single time step (5 - 10 % faster than other options);
    - `1` achieves machine accuracy for energy conservation by multiplying particle's momentum by a number close to 1, while keeping the precomputed change of electric field unchanged;
    - `2` (default) achieves machine accuracy for energy conservation by multiplying precomputed change of electric field by a number close to 1, while keeping particle's momentum unchanged.

Conventions, assumptions and properties
--

- CGS units are used for all dimensional quantities
- The components of vectors are enumerated: $x$, $y$, $z$ correspond to `[0]`, `[1]`, `[2]`, respectively
- The sizes of the grid `nx`, `ny` and `nz` must be powers of two (to facilitate FFT)
- 2D simulation is enabled by setting `nz=1`, 1D simulations is enabled by `nz=1` and `ny=1`
- The topology of space is toroidal (periodic boundary conditions)
- Before adding particles the field can be advanced over an arbitrary time, however particles should not traverse a distance larger then one spatial step over a single time step
- All solvers are deterministic (results are exactly reproducible)
