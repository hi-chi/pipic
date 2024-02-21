import pipic
from pipic.tools import *
from pipic.extensions import downsampler_gonoskov2022
from numba import cfunc, carray
import matplotlib.pyplot as plt
import math, numpy, os, time


#===========================SIMULATION INITIALIZATION===========================
temperature = 1e-6*electron_mass*light_velocity**2
density = 1e+18
debye_length = math.sqrt(temperature/(4*math.pi*density*electron_charge**2))
plasma_period = math.sqrt(math.pi*electron_mass/(density*electron_charge**2))
L = 128*debye_length
xmin, xmax = -L, L
nx, ny, ppc = 64, 32, 100
ymin, ymax = -L/2, L/2
stepx, stepy = (xmax - xmin)/nx, (ymax - ymin)/ny
time_step = plasma_period/64

#---------------------setting solver and simulation region---------------------- 
sim = pipic.init(solver='ec', nx=nx, xmin=xmin, xmax=xmax, ny=ny, ymin=ymin, ymax=ymax)
sim.set_rng_seed(1)
#------------------------------adding electrons---------------------------------
@cfunc(add_particles_callback)
def density_callback(r, data_double, data_int):  # callback function
    return density*(abs(r[0] - r[1]) < L/4)#*2*abs(r[1])/L
sim.add_particles(name='electron', number=int(ppc*nx*ny/8.0),
                  charge=electron_charge, mass=electron_mass,
                  temperature=temperature, density=density_callback.address)

#------------------------------adding extension---------------------------------
downsampler = downsampler_gonoskov2022.handler(sim.ensemble_data(), sim.get_type_index('electron'),
                                               preserve_energy=False, preserve_momentum=False, cap=6)
sim.add_handler(name=downsampler_gonoskov2022.name, subject='cells',
                handler=downsampler)

#=================================OUTPUT========================================
fig, axs = plt.subplots(2, constrained_layout=True)

#-------------preparing output for electron distrribution f(x, px)--------------
Ne = numpy.zeros((ny, nx), dtype=numpy.double)

dxdy = (xmax - xmin)*(ymax - ymin)/(nx*ny)

@cfunc(particle_loop_callback)
def Ne_cb(r, p, w, id, data_double, data_int):   
    data = carray(data_double, Ne.shape, dtype=numpy.double)
    # CIC cell-wise counting (2D):
    ix = int(Ne.shape[1]*(r[0] - 0.5*stepx - xmin)/(xmax - xmin) + 2) - 2
    iy = int(Ne.shape[0]*(r[1] - 0.5*stepy - ymin)/(ymax - ymin) + 2) - 2
    gx = (((r[0] - 0.5*stepx) - xmin) - stepx*ix)/stepx
    gy = (((r[1] - 0.5*stepy) - ymin) - stepy*iy)/stepy
    ix = (ix + nx) % nx
    iy = (iy + ny) % ny
    data[ny - 1 - iy, ix] += gx*gy*w[0]/dxdy
    ix += 1
    ix = (ix + nx) % nx
    data[ny - 1 - iy, ix] += (1 - gx)*gy*w[0]/dxdy
    iy += 1
    iy = (iy + ny) % ny
    data[ny - 1 - iy, ix] += (1 - gx)*(1 - gy)*w[0]/dxdy
    ix -= 1
    ix = (ix + nx) % nx
    data[ny - 1 - iy, ix] += gx*(1 - gy)*w[0]/dxdy


def plot_initial_density():
    Ne.fill(0)
    sim.particle_loop(name='electron', handler=Ne_cb.address,
                      data_double=addressof(Ne))
    axs[0].set_title('Initial density')
    axs[0].set(ylabel='$y$ (cm)')
    axs[0].xaxis.set_ticklabels([])
    plot0 = axs[0].imshow(Ne, vmin=0, vmax=2*density,
                      extent=[xmin, xmax, ymin, ymax], interpolation='none',
                      aspect='auto', cmap='YlOrBr')
    fig.colorbar(plot0, ax=axs[0], location='right')
    
axs[1].set_title('Density after downsampling')
axs[1].set(ylabel='$y$ (cm)')
axs[1].xaxis.set_ticklabels([])
plot1 = axs[1].imshow(Ne, vmin=0, vmax=2*density,
                    extent=[xmin, xmax, ymin, ymax], interpolation='none',
                    aspect='auto', cmap='YlOrBr')
fig.colorbar(plot1, ax=axs[1], location='right')
def plot_density():
    Ne.fill(0)
    sim.particle_loop(name='electron', handler=Ne_cb.address,
                      data_double=addressof(Ne))
    plot1.set_data(Ne)
    
#===============================SIMULATION======================================
outputFolder = 'output_downsampler_gonoskov2022'
if not os.path.exists(outputFolder):
    os.makedirs(outputFolder)
plot_initial_density()
print('Before downsampling N_macro = ', sim.get_number_of_particles())
print('CIC density at x = y = 0: ', Ne[ny - 1 - int(ny/2), int(nx/2)])
sim.advance(time_step=0.0, number_of_iterations=1)
print('After downsampling N_macro = ', sim.get_number_of_particles())
sim.advance(time_step=0.0, number_of_iterations=1)
print('After extra iteration (to remove zero-weight particles in processed cells) N_macro = ', sim.get_number_of_particles())
plot_density()
print('CIC density at x = y = 0: ', Ne[ny - 1 - int(ny/2), int(nx/2)])
fig.savefig(outputFolder + '/im0.png')