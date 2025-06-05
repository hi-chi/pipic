import pipic
from pipic import consts, types
import numpy as np
from numba import cfunc, carray
from pipic.extensions import moving_window
import matplotlib.pyplot as plt

# --------------------- simulation variables ----------------------
temperature = 1e-6 * consts.electron_mass * consts.light_velocity**2
density = 1e18
debye_length = np.sqrt(temperature / (4 * np.pi * density * consts.electron_charge**2))
plasma_period = np.sqrt(np.pi * consts.electron_mass / (density * consts.electron_charge**2))
l = 128 * debye_length
zmin, zmax = -l / 2, l / 2
ymin, ymax = -l / 2, l / 2
xmin, xmax = -l / 2, l / 2
field_amplitude = 0.01 * 4 * np.pi * (xmax - xmin) * consts.electron_charge * density
nx, ny, nz = 16, 16, 128
timestep = plasma_period / 64
simulation_steps = int(2 * plasma_period / timestep)

# --------------------- initial state functions ----------------------
@cfunc(types.add_particles_callback)
def density_profile(r, data_double, data_int):
   return density

@cfunc(types.field_loop_callback)
def initial_field(ind, r, E, B, data_double, data_int):
    E[2] = field_amplitude * np.sin(np.pi * r[2]/ zmax)

# --------------------- initialize simulation ----------------------
sim=pipic.init(solver='fourier_boris', # using fourier_boris solver
               xmin=xmin,xmax=xmax,
               ymin=ymin,ymax=ymax,
               zmin=zmin,zmax=zmax,
               nx=nx,ny=ny,nz=nz)

# add field according to initial_field
sim.field_loop(handler=initial_field.address) # setting initial field

# add particles according to density_profile
sim.add_particles(name='particle_name',
                  number= nx*ny*nz*1, # total number of particles to add (1 particle per cell)
                  density=density_profile.address,
                  charge=-consts.electron_charge,
                  mass=consts.electron_mass,
                  temperature=temperature,)



# --------------------- add handler for absorbing boundaries ----------------------
data_int = np.zeros((1, ), dtype=np.intc) 
'''
handler = moving_window.handler(sim.ensemble_data(),
                                thickness=10,
                                particles_per_cell=1,
                                temperature=temperature,
                                density=density_profile.address,)
sim.add_handler(name=moving_window.name, 
                subject='particle_name,cells',  #apply handler to both cells and particles
                handler=handler,
                data_int=pipic.addressof(data_int),) # using data_int to pass the iteration number
'''
# --------------------- add handler for reading field and particle loops ----------------------
field_dd = np.zeros((nx, nz), dtype=np.double)  # array for saving Ez-field
@cfunc(types.field_loop_callback)
def field_callback(ind, r, E, B, data_double, data_int):
    # read Ez in the xz plane at y=0
    data = carray(data_double, field_dd.shape, dtype=np.double)
    if ind[1] == ny // 2:
        data[ind[0],ind[2]] = E[2]

particle_dd = np.zeros((100, 100), dtype=np.double)  # array for saving particle (integrated) phase-space
pmin = -np.sqrt(consts.electron_mass * temperature)*5 # minimum momentum
pmax = np.sqrt(consts.electron_mass * temperature)*5 # maximum momentum
@cfunc(types.particle_loop_callback)
def particle_callback(r, p, w, id, data_double, data_int):
    # save particle momentum and position
    data = carray(data_double, particle_dd.shape, dtype=np.double)
    ip = int(particle_dd.shape[0] * (p[2] - pmin) / (pmax - pmin))
    iz = int(particle_dd.shape[1] * (r[2] - zmin) / (zmax - zmin))
    if ip >= 0 and ip < particle_dd.shape[0] and iz < particle_dd.shape[1]:
        data[ip, iz] += w[0] / (nx * ny * nz)



# --------------------- run simulation ----------------------
for i in range(simulation_steps):    
    
    data_int[0] = i # updating data int so the handler gets the right iteration number
    sim.advance(time_step=timestep,number_of_iterations=1)    

     # read and plot Ez-field and particle phase-space every 10 iterations
    if i % 10 == 0:
        sim.field_loop(handler=field_callback.address, 
                        data_double=pipic.addressof(field_dd),
                        use_omp=True)
        
        sim.particle_loop(name='particle_name', 
                            handler=particle_callback.address, 
                                data_double=pipic.addressof(particle_dd))
        
        # plot Ez-field and particle phase-space        
        fig,ax = plt.subplots(2, 1, constrained_layout=True)
        ax[0].imshow(field_dd, 
                     extent=[zmin, zmax,xmin, xmax], 
                     aspect='auto', origin='lower', 
                     cmap='coolwarm')
        ax[1].imshow(particle_dd, 
                     extent=[zmin, zmax,pmin, pmax], 
                     aspect='auto', origin='lower', 
                     cmap='Reds')
        fig.savefig(f'output_{i//10:04d}.png')
