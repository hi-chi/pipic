import pipic
from pipic import consts, types
import numpy as np
from numba import cfunc, carray
import matplotlib.pyplot as plt
import os

# Simulation variables (cgs units)
temperature = 1e-6 * consts.electron_mass * consts.light_velocity**2
density = 1e18
debye_length = np.sqrt(temperature / (4 * np.pi * density * consts.electron_charge**2))
plasma_period = np.sqrt(np.pi * consts.electron_mass / (density * consts.electron_charge**2))
l = 128 * debye_length
xmin, xmax = -l / 2, l / 2
field_amplitude = 0.01 * 4 * np.pi * (xmax - xmin) * consts.electron_charge * density
nx = 128
timestep = plasma_period / 64

# initialize simulation 
sim=pipic.init(solver='ec2', # using energy-conserving (ec) solver, other options include 'fourier_boris'
               xmin=xmin,xmax=xmax,
               nx=nx)

@cfunc(types.field_loop_callback)
def initial_field(ind, r, E, B, data_double, data_int):
    E[0] = field_amplitude * np.sin(4*np.pi * r[0]/ (xmax-xmin))

# add field according to initial_field
sim.field_loop(handler=initial_field.address) # setting initial field

# Define functions for initiating the simulation
@cfunc(types.add_particles_callback)
def density_profile(r, data_double, data_int):
    return density

# add particles according to density_profile
sim.add_particles(name='particle_name',
                  number= nx*200, # total number of particles to add 
                  density=density_profile.address,
                  charge=consts.electron_charge,
                  mass=consts.electron_mass,
                  temperature=temperature,)



# define functions and arrays for reading and saving field and particle phase space 
# For reading Ex field 
field_dd = np.zeros((nx,), dtype=np.double)  # array for saving Ez-field
@cfunc(types.field_loop_callback)
def field_callback(ind, r, E, B, data_double, data_int):
    # read Ez in the xz plane at y=0
    data = carray(data_double, field_dd.shape, dtype=np.double)
    data[ind[0]] = E[0]

# For reading particle phase space
particle_dd = np.zeros((64, nx), dtype=np.double)  # array for saving particle (integrated) phase-space
pmin = -np.sqrt(consts.electron_mass * temperature)*5 # minimum momentum
pmax = np.sqrt(consts.electron_mass * temperature)*5 # maximum momentum
dp = (pmax - pmin) / particle_dd.shape[0] # momentum step
dx = (xmax - xmin) / particle_dd.shape[1] # position step
@cfunc(types.particle_loop_callback)
def particle_callback(r, p, w, id, data_double, data_int):
    # save particle momentum and position
    data = carray(data_double, particle_dd.shape, dtype=np.double)
    ip = int(particle_dd.shape[0] * (p[0] - pmin) / (pmax - pmin))
    ix = int(particle_dd.shape[1] * (r[0] - xmin) / (xmax - xmin))
    if ip >= 0 and ip < particle_dd.shape[0] and ix < particle_dd.shape[1] and ix >= 0:
        data[ip, ix] += w[0] / (dx * dp)  # normalize by dx, dp


# initialize plot
fig, ax = plt.subplots(2, 1, constrained_layout=True)

x_axis = np.linspace(xmin, xmax, nx)
Ex_plot = ax[1].plot(x_axis,field_dd)[0]
ax[1].set_ylim(field_amplitude, -field_amplitude)
xpx_plot = ax[0].imshow(particle_dd,  
             extent=[xmin, xmax, pmin, pmax], 
             aspect='auto', origin='lower', 
             cmap='YlOrBr',vmin=0, vmax= 6 * density / ( 2 * pmax ),
             interpolation = 'none')


# set titles
ax[0].set_title('Plasma oscillations')
ax[1].set_xlabel('x (cm)')
ax[0].set_ylabel('$p_x$ (cm g/s)')
ax[1].set_ylabel('$E_x$ (StatV/cm)')


# ===============================SIMULATION======================================
simulation_steps = int(8 * plasma_period / timestep)
frames = simulation_steps // 8 # number of plots to generate
save_to = "./output/"
# mkdir save_to directory if not exist
if not os.path.exists(save_to):
    os.makedirs(save_to)

for i in range(frames):
    sim.advance(time_step=timestep, number_of_iterations=8,use_omp=True)
    sim.field_loop(handler=field_callback.address, 
                   data_double=pipic.addressof(field_dd),
                   use_omp=True)
    
    particle_dd.fill(0)
    sim.particle_loop(name='particle_name', 
                      handler=particle_callback.address, 
                      data_double=pipic.addressof(particle_dd))
    Ex_plot.set_ydata(field_dd)
    xpx_plot.set_data(particle_dd)
    plt.savefig(save_to + f'plasma_oscillation_{i:03d}.png', dpi=150)

