#Basic setup for a laser pulse interation with a solid-density plasma layer 
#for results see fig. 6 in arXiv:2302.01893
import sys
import pipic
from pipic import consts,types
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from numba import cfunc, carray, types as nbt
import moving_window

#===========================SIMULATION INITIALIZATION===========================
wavelength = 1e-4 #units? cm?
pulseWidth = 10e-15*consts.light_velocity
nx, xmin, xmax = 2**7, -20*wavelength, 15*wavelength
ny, ymin, ymax = 2**6, -8*pulseWidth, 8*pulseWidth
nz, zmin, zmax = 2**6, -8*pulseWidth, 8*pulseWidth
dx, dy, dz = (xmax - xmin)/nx, (ymax - ymin)/ny, (zmax - zmin)/nz
timestep = 0.5*dx/consts.light_velocity
thickness = 20 # thickness (in dx) of the area where the density and field is restored/removed 

#---------------------setting solver and simulation region----------------------
sim=pipic.init(solver='ec',nx=nx,ny=ny,nz=nz,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,zmin=zmin,zmax=zmax)
#---------------------------setting field of the pulse--------------------------

amplitude_a0 = 100.
omega = 2 * np.pi * consts.light_velocity / wavelength
a0 = consts.electron_mass * consts.light_velocity * omega / (-consts.electron_charge)
fieldAmplitude = amplitude_a0 * a0
print(a0)
PulseDuration = 10e-15 #s
k = 2*np.pi / wavelength
spotsize = PulseDuration*consts.light_velocity #1e-2*1e-4 #cm
R = 1e-2 # wavefront curvature #0./(PulseDuration*consts.light_velocity)

N_cr = 1e-1*consts.electron_mass * omega ** 2 / (4 * np.pi * consts.electron_charge ** 2)
density = N_cr
debye_length = .1*wavelength/64.0 #>>3/(4*pi*density), <dx ???
temperature = 4 * np.pi * density * (consts.electron_charge ** 2) * debye_length ** 2
particles_per_cell = 5

pwl = np.sqrt(4*np.pi*density*consts.electron_charge**2/consts.electron_mass)
print('lambda:',2*np.pi*consts.light_velocity/pwl)
print('pulse_length:',PulseDuration*consts.light_velocity)

fp = f'data.txt'


@cfunc(types.field_loop_callback)
def initiate_field_callback(ind, r, E, B, data_double, data_int):
    if data_int[0] == 0:        

        rx_curv = r[0]-(r[1]**2+r[2]**2)/(2*R)
        # y-polarized 
        E[1] = fieldAmplitude*np.exp(-(rx_curv**2/(PulseDuration*consts.light_velocity)**2 + 
                                      (r[1]**2 + r[2]**2)/spotsize**2))*np.cos(k*rx_curv)
        B[2] = fieldAmplitude*np.exp(-(rx_curv**2/(PulseDuration*consts.light_velocity)**2 + 
                                       (r[1]**2 + r[2]**2)/spotsize**2))*np.cos(k*rx_curv)        
        # z-polarized
        E[2] = fieldAmplitude*np.exp(-(rx_curv**2/(PulseDuration*consts.light_velocity)**2 +
                                      (r[1]**2 + r[2]**2)/spotsize**2))*np.cos(k*rx_curv)
        B[1] = -fieldAmplitude*np.exp(-(rx_curv**2/(PulseDuration*consts.light_velocity)**2 +
                                       (r[1]**2 + r[2]**2)/spotsize**2))*np.cos(k*rx_curv)
 

@cfunc(types.add_particles_callback)
def initiate_density(r, data_double, data_int):# Just to add electrons to particle list

    rollback = np.floor(data_int[0]*timestep*consts.light_velocity/dx)   

    r_rel = xmin + dx*(rollback%nx)

    r_min = r_rel 
    r_max = r_rel + dx
    if (r[0] > r_min and r[0] < r_max) or (r[0] > xmax - (xmin - r_min)) or (r[0] < xmin + (r_max - xmax)):
        return density
    else:
        return 0 


@cfunc(types.field_loop_callback)
def remove_field(ind, r, E, B, data_double, data_int):
    rollback = np.floor(data_int[0]*timestep*consts.light_velocity/dx)
    r_rel = xmin + dx*(rollback%nx) 
    
    r_min = r_rel - thickness*dx
    r_max = r_rel #+ thickness*dx
    if (r[0] > r_min and r[0] < r_max) or (r[0] > xmax - (xmin - r_min)) or (r[0] < xmin + (r_max - xmax)): 
        E[1] = 0
        B[2] = 0
        E[2] = 0 
        B[1] = 0
        E[0] = 0
        B[0] = 0 
#=================================OUTPUT========================================

#-------------------------preparing output of fields (x)-----------------------------
Ey = np.zeros((nx,ny,nz), dtype=np.double) 
rho_rad = np.zeros((nx,ny//2), dtype=np.double) 
ps = np.zeros((nx,nx//2), dtype=np.double) 
pmin = -1
pmax = 1
r = np.arange(dy/2,ymax+dy/2,dy)
norm_rad = (np.pi*((r+dy/2)**2-(r-dy/2)**2))[np.newaxis,:]

#------------------get functions-----------------------------------------------

@cfunc(types.particle_loop_callback)
def get_radial_density(r, p, w, id, data_double, data_int):   
    ix = int(rho_rad.shape[0]*(r[0] - xmin)/(xmax - xmin))
    ir = int(rho_rad.shape[1]*np.sqrt(r[1]**2+r[2]**2)/ymax)
    data = carray(data_double, rho_rad.shape, dtype=np.double)
    
    if ir < rho_rad.shape[1] and ix < rho_rad.shape[0]:
        data[ix, ir] += w[0]/dx

@cfunc(types.particle_loop_callback)
def get_phase_space(r, p, w, id, data_double, data_int):   
    ix = int(ps.shape[0]*(r[0] - xmin)/(xmax - xmin))
    ip = int(ps.shape[1]*(p[0] - pmin)/(pmax - pmin))
    data = carray(data_double, ps.shape, dtype=np.double)
    
    if ip>=0 and ip < ps.shape[1] and ix < ps.shape[0]:
        data[ix, ip] = p[0]#w[0]/(dx*dy*dz)


@cfunc(types.field_loop_callback)
def get_field_Ey(ind, r, E, B, data_double, data_int):
    _E = carray(data_double, Ey.shape, dtype=np.double)
    _E[ind[0], ind[1], ind[2]] = E[1]

def load_fields():
    sim.field_loop(handler=get_field_Ey.address, data_double=pipic.addressof(Ey))

#------------------- plot fields -------------------------------------------
fig, axs = plt.subplots(3,1, constrained_layout=False)

E_plot = axs[0].imshow(Ey[:,:,nz//2].T,vmin=-fieldAmplitude, vmax=fieldAmplitude, 
        extent=[xmin,xmax,ymin,ymax], interpolation='none', \
        aspect='auto', cmap='inferno_r', origin='lower')
fig.colorbar(E_plot,ax=axs)
rho_plot = axs[1].imshow(rho_rad[:,:].T,vmin=0, vmax=density*2, 
        extent=[xmin,xmax,0,ymax/2], interpolation='none', \
        aspect='auto', cmap='inferno_r', origin='lower')

phase_space_plot = axs[2].imshow(ps.T,vmin=0, vmax=density*2*1e2, 
        extent=[xmin,xmax,pmin,pmax], interpolation='none', \
        aspect='auto', cmap='inferno_r', origin='lower')

#===============================SIMULATION======================================

data_int = np.zeros((1, ), dtype=np.intc) # data for passing the iteration number
s = 8000 #nb of iterations
window_speed = consts.light_velocity #speed of moving window

#-----------------------adding the handler of extension-------------------------

#for passing variables to the handler  
data_double = np.zeros((4, ), dtype=np.double)
data_double[0] = density

#data_double=np.zeros((1,), dtype=np.double)



density_handler_adress = moving_window.handler(#density=density,
                                               thickness=thickness,
                                               number_of_particles=particles_per_cell,
                                               temperature=temperature,)

sim.add_handler(name=moving_window.name, 
                subject='electron,cells',
                handler=density_handler_adress,
                data_int=pipic.addressof(data_int),
                data_double=pipic.addressof(data_double))

#-----------------------initiate field and plasma-------------------------
sim.field_loop(handler=initiate_field_callback.address, data_int=pipic.addressof(data_int),
                   use_omp=True)

sim.add_particles(name='electron', number=int(ny*nz*nx),#particles_per_cell),
                 charge=consts.electron_charge, mass=consts.electron_mass,
                 temperature=temperature, density=initiate_density.address,
                 data_int=pipic.addressof(data_int))


#-----------------------run simulation-------------------------
for i in range(1,1000):
    
    data_int[0] = i 

    sim.advance(time_step=timestep, number_of_iterations=1,use_omp=True)
    
    sim.field_loop(handler=remove_field.address, data_int=pipic.addressof(data_int),
                   use_omp=True)

    if i==200:
        data_double[0] = density*5e-1

    if i%10==0:
        rho_rad.fill(0)
        ps.fill(0)

        sim.particle_loop(name='electron', handler=get_radial_density.address,
                      data_double=pipic.addressof(rho_rad))
        sim.particle_loop(name='electron', handler=get_phase_space.address,
                      data_double=pipic.addressof(ps))
        load_fields()

        E = Ey[:,:,nz//2]
        roll_back =  int(np.floor(i*timestep*window_speed/dx))  
        E = np.roll(E,-roll_back,axis=0)
        E_plot.set_data((E.T))
        
        rho_rad = np.roll(rho_rad,-roll_back,axis=0)/norm_rad
        rho_plot.set_data(rho_rad.T)

        ps = np.roll(ps,-roll_back,axis=0)
        phase_space_plot.set_data(ps.T)
        print(ps.max(),ps.min())

        rho_ = rho_rad[-1,-10:].mean()
        print(i,f':{rho_,density}')
        fig.savefig('im' + str(i) + '.png')

#np.savetxt(fp,u)
