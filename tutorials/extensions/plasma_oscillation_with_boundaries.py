''' Description: this file demonstrates the use of the absorbing boundaries extension
to simulate plasma oscillations in a 1D plasma with open boundaries. The absorbing
boundaries extension removes particles leaving the simulation box and adds new
particles according to a specified density profile, simulating an open plasma system.
The script uses the absorbing boundary extension installable with pipic.
For more details on installtion and usage of the extension, please refer to:
https://'''
import pipic
from pipic import consts, types
import numpy as np
from numba import cfunc, carray
# Note this script uses the absorbing boundaries extension installed with pipic
from pipic.extensions import absorbing_boundaries
import matplotlib.pyplot as plt
import os

# =============================================================================
# SIMULATION SETUP
# =============================================================================

# Plasma parameters (CGS units)
# Temperature in units of [erg] (T=T_k * k_B, where T_k is temperature in Kelvin 
# and k_B is the Boltzmann constant [erg/K])
temperature = 1e-6 * consts.electron_mass * consts.light_velocity**2
density = 1e18 # particle number density in units of [1/cm^3]
debye_length = np.sqrt(
    temperature / (4 * np.pi * density * consts.electron_charge**2)
) # in units of [cm]
plasma_period = np.sqrt(
    np.pi * consts.electron_mass / (density * consts.electron_charge**2)
) # in units of [s]

# Simulation box parameters
l = 128 * debye_length # simulation box length
xmin, xmax = -l / 2, l / 2
nx = 128 # number of cells
timestep = plasma_period / 64 

# Electric field parameters
field_amplitude = (
    0.01 * 4 * np.pi * (xmax - xmin) * consts.electron_charge * density
)
dmin, dmax = -l / 4, l / 4  # region to apply field

# Initialize simulation with energy-conserving solver
sim = pipic.init(solver="ec2", xmin=xmin, xmax=xmax, nx=nx)

# =============================================================================
# INITIAL FIELD
# =============================================================================
@cfunc(types.field_loop_callback)
def initial_field(ind, r, E, B, data_double, data_int):
    """Applies a sinusoidal initial electric field in the region dmin < x < dmax."""
    if dmin < r[0] < dmax:
        E[0] = (
            field_amplitude
            * np.sin(4 * np.pi * r[0] / (xmax - xmin))
            * np.exp(-r[0] ** 2 / (2 * (0.2 * l) ** 2))
        )

# Apply initial field
sim.field_loop(handler=initial_field.address)

# =============================================================================
# PARTICLES
# =============================================================================
@cfunc(types.add_particles_callback)
def density_profile(r, data_double, data_int):
    """Uniform density profile."""
    return density

# Add particles according to density_profile
sim.add_particles(
    name="particle_name", # name of the particle species
    number=nx*100,  # total number of particles to add
    density=density_profile.address, # density profile function
    charge=consts.electron_charge, # particle charge
    mass=consts.electron_mass, # particle mass
    temperature=temperature, # particle temperature
)

# =============================================================================
# ABSORBING BOUNDARIES EXTENSION
# =============================================================================
data_int = np.zeros((1,), dtype=np.intc)  # used to pass iteration number
boundary_size = xmax / 2

sim.add_handler(
    name=absorbing_boundaries.name,
    subject="particle_name,cells",  # apply to both particles and cells
    # particle handler
    handler=absorbing_boundaries.handler(
        # pass address to cell and particle data
        sim.ensemble_data(),
        # pass address to simulation box geometry
        sim.simulation_box(),
        # pass adress to density profile function
        density_profile=density_profile.address,
        # size of boundary (in cm)
        boundary_size=boundary_size,
        # size of temperature of particles to be added
        temperature=temperature,
        # number of particles to be added per cell 
        particles_per_cell=100,
    ),
    field_handler=absorbing_boundaries.field_handler(
        # pass address to simulation box geometry
        sim.simulation_box(), 
        # pass size of time step
        timestep=timestep
    ),
    # pass address to data_int for passing iteration number 
    # (see RUN SIMULATION section)
    data_int=pipic.addressof(data_int),
)

# =============================================================================
# DIAGNOSTICS
# =============================================================================

# --- Field diagnostic ---
# array for saving Ex field
field_dd = np.zeros((nx,), dtype=np.double)  # array for saving Ex field
# callback function for field diagnostic
@cfunc(types.field_loop_callback)
def field_callback(ind, r, E, B, data_double, data_int):
    """Store Ex."""
    data = carray(data_double, field_dd.shape, dtype=np.double)
    data[ind[0]] = E[0]

# --- Particle phase space diagnostic ---
# array for saving particle (integrated) phase-space
particle_dd = np.zeros((64, nx), dtype=np.double)
# momentum range for phase space
pmin = -5 * np.sqrt(consts.electron_mass * temperature)
pmax = 5 * np.sqrt(consts.electron_mass * temperature)
# momentum and position steps
dp = (pmax - pmin) / particle_dd.shape[0]
dx = (xmax - xmin) / particle_dd.shape[1]
# callback function for particle phase space diagnostic
@cfunc(types.particle_loop_callback)
def particle_callback(r, p, w, id, data_double, data_int):
    """Calculate particle momentum-position phase space density."""
    data = carray(data_double, particle_dd.shape, dtype=np.double)
    ip = int(particle_dd.shape[0] * (p[0] - pmin) / (pmax - pmin))
    ix = int(particle_dd.shape[1] * (r[0] - xmin) / (xmax - xmin))
    if 0 <= ip < particle_dd.shape[0] and 0 <= ix < particle_dd.shape[1]:
        data[ip, ix] += w[0] / (dx * dp)  # normalize

# =============================================================================
# PLOTTING SETUP
# =============================================================================
fig, ax = plt.subplots(2, 1, constrained_layout=True)

# Field plot
x_axis = np.linspace(xmin, xmax, nx)
Ex_plot = ax[1].plot(x_axis, field_dd)[0]
ax[1].set_ylim(field_amplitude, -field_amplitude)

# Phase space plot
xpx_plot = ax[0].imshow(
    particle_dd,
    extent=[xmin, xmax, pmin, pmax],
    aspect="auto",
    origin="lower",
    cmap="YlOrBr",
    vmin=0,
    vmax=6 * density / (2 * pmax),
    interpolation="none",
)

# Labels
ax[0].set_title("Plasma oscillations")
ax[1].set_xlabel("x (cm)")
ax[0].set_ylabel("$p_x$ (cm g/s)")
ax[1].set_ylabel("$E_x$ (StatV/cm)")

# =============================================================================
# RUN SIMULATION
# =============================================================================
simulation_steps = int(8 * plasma_period / timestep)
# Number of figures to save (every 8 steps)
figures = simulation_steps // 8
# Directory to save figures
save_to = "./output/"

# Create output directory
if not os.path.exists(save_to):
    os.makedirs(save_to)

for i in range(figures):
    # Update iteration counter for absorbing boundary extension
    data_int[0] = i * 8

    # Advance simulation (use_omp=True to enable OpenMP parallelization)
    sim.advance(time_step=timestep, number_of_iterations=8, use_omp=True)

    # Collect diagnostics
    sim.field_loop(
        handler=field_callback.address,
        # pass address to field_dd array for storing Ex
        data_double=pipic.addressof(field_dd),
        # enable OpenMP parallelization
        use_omp=True,
    )
    particle_dd.fill(0)
    sim.particle_loop(
        # name of particle species to be processed
        name="particle_name",
        handler=particle_callback.address,
        # pass address to particle_dd array for storing phase space
        data_double=pipic.addressof(particle_dd),
    )

    # Update plots
    Ex_plot.set_ydata(field_dd)
    xpx_plot.set_data(particle_dd)

    # Save figure
    plt.savefig(save_to + f"plasma_oscillation_{i:03d}.png", dpi=150)
