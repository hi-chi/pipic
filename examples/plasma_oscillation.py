# This is the file used for producing fig. 2 in arXiv:2302.01893.
import pipic
from pipic import types
from pipic.consts import *
from pipic.tools import *
import matplotlib.pyplot as plt
import numpy as np
from numba import cfunc, carray
import os


#===========================SIMULATION INITIALIZATION===========================
temperature = 1e-6 * electron_mass * light_velocity**2
density = 1e+18
debye_length = np.sqrt(temperature / (4*np.pi * density * electron_charge**2))
plasma_period = np.sqrt(np.pi * electron_mass / (density * electron_charge**2))
xmin, xmax = -64*debye_length, 64*debye_length
field_amplitude = 0.01 * 4*np.pi * (xmax - xmin) * electron_charge * density
nx = 128 # number of cells
time_step = plasma_period/64

#---------------------setting solver and simulation region----------------------
sim = pipic.init(solver='ec', nx=nx, xmin=xmin, xmax=xmax) 
# by default ny=nz=1, YMax=ZMax=0.5, YMin=ZMin=-0.5

#------------------------------adding electrons---------------------------------
@cfunc(types.add_particles)
def density_callback(r, data_double, data_int):# callback function
    return density # can be any function of coordinate r[0] 
sim.add_particles(name='electron', number=nx*1000,
                 charge=-electron_charge, mass=electron_mass,
                 temperature=temperature, density=density_callback.address)

#---------------------------setting initial field-------------------------------
@cfunc(types.field_loop)
def setField_callback(ind, r, E, B, data_double, data_int):
    E[0] = field_amplitude*np.sin(2*np.pi*r[0]/(xmax - xmin))
sim.field_loop(handler=setField_callback.address)


#=================================OUTPUT========================================
fig, axs = plt.subplots(2, constrained_layout=True)

#-------------preparing output for electron distrribution f(x, px)--------------
xpx_dist = np.zeros((64, 128), dtype=np.double)
pxLim = 5*np.sqrt(temperature * electron_mass)
inv_dx_dpx = (xpx_dist.shape[1]/(xmax - xmin))*(xpx_dist.shape[0]/(2*pxLim))

@cfunc(types.particle_loop)
def xpx_callback(r, p, w, id, data_double, data_int):   
    ix = int(xpx_dist.shape[1]*(r[0] - xmin)/(xmax - xmin))
    iy = int(xpx_dist.shape[0]*0.5*(1 + p[0]/pxLim))
    data = carray(data_double, xpx_dist.shape, dtype=np.double)
    if 0 <= iy < xpx_dist.shape[0]:
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

#-------------------------preparing output of Ex(x)-----------------------------
Ex = np.zeros((32, ), dtype=np.double)

@cfunc(types.it2r)
def Ex_it2r(it, r, data_double, data_int):
    r[0] = xmin + (it[0] + 0.5)*(xmax - xmin)/Ex.shape[0]

@cfunc(types.field2data)
def get_Ex(it, r, E, B, data_double, data_int):
    data_double[it[0]] = E[0]

axs[1].set_xlim([xmin, xmax])
axs[1].set_ylim([-field_amplitude, field_amplitude])
axs[1].set(xlabel='$x$ (cm)', ylabel='$E_x$ (cgs units)')
x_axis = np.linspace(xmin, xmax, Ex.shape[0])
sim.custom_field_loop(Ex.shape[0], Ex_it2r.address, get_Ex.address, address_of(Ex))
plot_Ex_, = axs[1].plot(x_axis, Ex)

def plot_Ex():
    sim.custom_field_loop(number_of_iterations=Ex.shape[0], it2r=Ex_it2r.address,
                        field2data=get_Ex.address, data_double=address_of(Ex))
    plot_Ex_.set_ydata(Ex)


#===============================SIMULATION======================================
outputFolder = 'plasma_oscillation_output'
if not os.path.exists(outputFolder):
   os.makedirs(outputFolder)
for i in range(32):
    sim.advance(time_step=time_step, number_of_iterations=2)
    plot_xpx()
    plot_Ex()
    fig.savefig(outputFolder + '/im' + str(i) + '.png')
    if i == 25:
        fig.savefig(outputFolder + '/fig2.pdf')
    print(i, '/', 32)