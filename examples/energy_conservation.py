#Demonstration of energy conservation (figs. 5 and 6 in arXiv:2302.01893)
import pipic
from pipic.tools import *
import matplotlib.pyplot as plt
import numpy as np
from numba import cfunc, carray


#===========================SIMULATION INITIALIZATION===========================
temperature = (2/3) * 1e-6 * electron_mass * light_velocity**2
density = 1e+18
debye_length = np.sqrt(temperature / (4*np.pi * density * electron_charge**2))
plasma_period = np.sqrt(np.pi * electron_mass / (density * electron_charge**2))
xmin, xmax = -612.4*debye_length, 612.4*debye_length
field_amplitude = 0.001 * 4*np.pi * (xmax - xmin) * electron_charge * density
nx = 32 # number of cells
iterationsPerPlasmaPeriod = 8
time_step = plasma_period/iterationsPerPlasmaPeriod
fig_label = 'ec2'

#---------------------setting solver and simulation region----------------------
if fig_label == 'boris':
    sim = pipic.init(solver='fourier_boris', nx=nx, xmin=xmin, xmax=xmax)
if fig_label == 'boris_no_div_cleaning':
    sim = pipic.init(solver='fourier_boris', nx=nx, xmin=xmin, xmax=xmax)
    sim.fourier_solver_settings(divergence_cleaning=False)
if fig_label == 'ec':
    sim = pipic.init(solver='ec', nx=nx, xmin=xmin, xmax=xmax)
if fig_label == 'ec2':
    sim = pipic.init(solver='ec2', nx=nx, xmin=xmin, xmax=xmax)

#------------------------------adding electrons---------------------------------
@cfunc(add_particles_callback)
def density_callback(r, data_double, data_int):# callback function
    return density # can be any function of coordinate r[0] 
sim.add_particles(name='electron', number=nx*100,
                 charge=-electron_charge, mass=electron_mass,
                 temperature=temperature, density=density_callback.address)

#---------------------------setting initial field-------------------------------
shift = np.pi/nx # a shift to mitigate aliasing

@cfunc(field_loop_callback)
def setField_callback(ind, r, E, B, data_double, data_int):
    E[0] = field_amplitude*np.sin(2*np.pi*r[0]/(xmax - xmin) + shift)

sim.field_loop(handler=setField_callback.address)


#=================================OUTPUT========================================
fig, axs = plt.subplots(2, constrained_layout=True)

#-------------preparing output for electron distrribution f(x, px)--------------
xpx_dist = np.zeros((64, 128), dtype=np.double)
pxLim = 5*np.sqrt(temperature * electron_mass)
inv_dx_dpx = (xpx_dist.shape[1]/(xmax - xmin))*(xpx_dist.shape[0]/(2*pxLim))

@cfunc(particle_loop_callback)
def xpx_callback(r, p, w, id, data_double, data_int):   
    ix = int(xpx_dist.shape[1]*(r[0] - xmin)/(xmax - xmin))
    iy = int(xpx_dist.shape[0]*0.5*(1 + p[0]/pxLim))
    data = carray(data_double, xpx_dist.shape, dtype=np.double)
    if iy >= 0 and iy < xpx_dist.shape[0]:
        data[iy, ix] += w[0]*inv_dx_dpx

axs[0].set_title('$\partial N / \partial x \partial p_x$ (s g$^{-1}$cm$^{-2}$)')
axs[0].set(ylabel='$p_x$ (cm g/s)')
axs[0].xaxis.set_ticklabels([])
plot0 = axs[0].imshow(xpx_dist, vmax=3*density/pxLim,
                      extent=[xmin,xmax,-pxLim,pxLim], interpolation='none',
                      aspect='auto', cmap='YlOrBr')
fig.colorbar(plot0, ax=axs[0], location='right')

def plot_xpx():
    xpx_dist.fill(0)
    sim.particle_loop(name='electron', handler=xpx_callback.address,
                     data_double=address_of(xpx_dist))
    plot0.set_data(xpx_dist)

#--------------------preparing output of kinetic energy-------------------------
mc2 = electron_mass * light_velocity**2
m2c2 = (electron_mass * light_velocity)**2
kineticEnergy = np.zeros((1, ), dtype=np.double)

@cfunc(particle_loop_callback)
def kineticEn_cb(r, p, w, id, data_double, data_int):   
    data_double[0] += w[0]*mc2*(np.sqrt(1+(p[0]**2+p[1]**2+p[2]**2)/m2c2)-1)

def getKineticEn(sim):
    kineticEnergy[0] = 0
    sim.particle_loop(name='electron', handler=kineticEn_cb.address,
                     data_double=address_of(kineticEnergy))
    return kineticEnergy[0]

#--------------------preparing output of field energy-------------------------
factor = ((sim.xmax-sim.xmin)/sim.nx)*((sim.ymax-sim.ymin)/sim.ny)* \
        ((sim.zmax-sim.zmin)/sim.nz)/(8*np.pi)
fieldEnergy = np.zeros((1, ), dtype=np.double)

@cfunc(field_loop_callback)
def fieldEn_cb(ind, r, E, B, data_double, data_int):
    data_double[0] += factor*(E[0]**2+E[1]**2+E[2]**2+B[0]**2+B[1]**2+B[2]**2)

def getFieldEnergy(sim):
    fieldEnergy[0] = 0
    sim.field_loop(handler=fieldEn_cb.address, data_double=address_of(fieldEnergy))
    return fieldEnergy[0]


#===============================SIMULATION======================================
nPlasmaOscillations = 100 # number of plasma oscillations to be simulated
experimentDuration = nPlasmaOscillations*plasma_period
nIterations = nPlasmaOscillations*iterationsPerPlasmaPeriod
timeP = np.zeros((0, ), dtype=np.double)
Efield = np.zeros((0, ), dtype=np.double)
Etot = np.zeros((0, ), dtype=np.double)
E0 = 0
for i in range(nPlasmaOscillations*iterationsPerPlasmaPeriod):
    print(i, "/", nPlasmaOscillations*iterationsPerPlasmaPeriod)
    Ek = getKineticEn(sim)
    Ef = getFieldEnergy(sim)
    if i == 0:
        E0 = Ek+Ef  
    timeP = np.append(timeP, i*time_step/plasma_period)
    Efield = np.append(Efield, Ef/E0)
    Etot = np.append(Etot, (Ek+Ef)/E0)
    print('relative energy deviation = ', ((Ek+Ef)-E0)/E0)
    if Ef > 2*E0:
        break
    sim.advance(time_step=time_step)


#=====================SAVING SIMULATION RESULT TO FILE==========================
Name = fig_label + '_' + str(iterationsPerPlasmaPeriod)
with open(Name + '.npy', 'wb') as f:
    np.save(f, Etot)
    np.save(f, Efield)
    np.save(f, timeP)
