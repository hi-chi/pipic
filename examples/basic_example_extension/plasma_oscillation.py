import pipic
from pipic import consts, types
import numpy as np
from numba import cfunc, carray
import extension
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
nx, ny, nz = 2, 2, 128
timestep = plasma_period / 64
simulation_steps = int(8 * plasma_period / timestep)

dmin = - l/4
dmax = + l/4

# --------------------- initial state functions ----------------------
@cfunc(types.add_particles_callback)
def density_profile(r, data_double, data_int):
    #if dmin < r[2] < dmax:
    return density
    #else:
    #    return 0.0

@cfunc(types.field_loop_callback)
def initial_field(ind, r, E, B, data_double, data_int):
    if dmin < r[2] < dmax:
        E[2] = field_amplitude * np.sin(4*np.pi * r[2]/ (zmax-zmin))

# --------------------- initialize simulation ----------------------
sim=pipic.init(solver='ec', # using fourier_boris solver
               xmin=xmin,xmax=xmax,
               ymin=ymin,ymax=ymax,
               zmin=zmin,zmax=zmax,
               nx=nx,ny=ny,nz=nz)

# add field according to initial_field
sim.field_loop(handler=initial_field.address) # setting initial field

# add particles according to density_profile
sim.add_particles(name='particle_name',
                  number= nx*ny*nz*500, # total number of particles to add (1 particle per cell)
                  density=density_profile.address,
                  charge=-consts.electron_charge,
                  mass=consts.electron_mass,
                  temperature=temperature,)

# --------------------- add handler ----------------------------------------------------------
data_int = np.zeros((1, ), dtype=np.intc) 
boundary_size = zmax/2 # size of the absorbing boundary 

sim.add_handler(name=extension.name, 
                subject='particle_name,cells',  #apply handler to both cells and particles
                handler=extension.handler(sim.ensemble_data(),
                                          sim.simulation_box(),
                                          density,
                                          boundary_size,
                                          temperature=temperature), # particle handler
                field_handler=extension.field_handler(sim.simulation_box(),boundary_size), # field handler
                data_int=pipic.addressof(data_int),) # using data_int to pass the iteration number

# --------------------- add handler for reading field and particle loops ----------------------
field_dd = np.zeros((nx, nz), dtype=np.double)  # array for saving Ez-field
@cfunc(types.field_loop_callback)
def field_callback(ind, r, E, B, data_double, data_int):
    # read Ez in the xz plane at y=0
    data = carray(data_double, field_dd.shape, dtype=np.double)
    if ind[1] == ny // 2:
        data[ind[0],ind[2]] = E[2]

particle_dd = np.zeros((64, nz), dtype=np.double)  # array for saving particle (integrated) phase-space
pmin = -np.sqrt(consts.electron_mass * temperature)*5 # minimum momentum
pmax = np.sqrt(consts.electron_mass * temperature)*5 # maximum momentum
dp = (pmax - pmin) / particle_dd.shape[0] # momentum step
dz = (zmax - zmin) / particle_dd.shape[1] # position step
@cfunc(types.particle_loop_callback)
def particle_callback(r, p, w, id, data_double, data_int):
    # save particle momentum and position
    data = carray(data_double, particle_dd.shape, dtype=np.double)
    ip = int(particle_dd.shape[0] * (p[2] - pmin) / (pmax - pmin))
    iz = int(particle_dd.shape[1] * (r[2] - zmin) / (zmax - zmin))
    if ip >= 0 and ip < particle_dd.shape[0] and iz < particle_dd.shape[1]:
        #data[ip, iz] += w[0] / (dz * dp) / (3*density/pmax) / (xmax - xmin) / (ymax - ymin)  # normalize by dz, dp and density
        data[ip, iz] += w[0] / (dz * dp) / density #/ (3*density/pmax) / (xmax - xmin) / (ymax - ymin)  # normalize by dz, dp and density

# initial plot
fig,ax = plt.subplots(2, 1, constrained_layout=True)
z_axis = np.linspace(zmin, zmax, nz)
izb = np.argwhere(z_axis > zmin + boundary_size)[0][0]
Ez_plot = ax[1].plot(z_axis[izb:-izb],field_dd[nx//2, izb:-izb])[0]
ax[1].set_ylim(-field_amplitude, field_amplitude)
zpz_plot = ax[0].imshow(particle_dd / (3/pmax) / (xmax - xmin) / (ymax - ymin), 
             extent=[zmin, zmax,pmin, pmax], 
             aspect='auto', origin='lower', 
             cmap='YlOrBr',vmin=0, vmax=1,
             interpolation = 'none')


# --------------------- run simulation ----------------------
for i in range(simulation_steps):    
    
    data_int[0] = i # updating data int so the handler gets the right iteration number
    sim.advance(time_step=timestep,number_of_iterations=1, use_omp=False)    

     # read and plot Ez-field and particle phase-space every 10 iterations
    if i % 10 == 0:
        sim.field_loop(handler=field_callback.address, 
                        data_double=pipic.addressof(field_dd),
                        use_omp=True)
        
        particle_dd.fill(0)
        sim.particle_loop(name='particle_name', 
                            handler=particle_callback.address, 
                                data_double=pipic.addressof(particle_dd))
        dx = (xmax - xmin) / nx
        dy = (ymax - ymin) / ny
        nnp = (2*pmax) / dp # number of momentum bins

        # plot Ez-field and particle phase-space        
        z_axis = np.linspace(zmin, zmax, nz)
        Ez_plot.set_ydata(field_dd[nx//2, izb:-izb])
        zpz_plot.set_data(particle_dd / (3/pmax) / (xmax - xmin) / (ymax - ymin))
        fig.savefig(f'./output_{i//10:04d}.png')
