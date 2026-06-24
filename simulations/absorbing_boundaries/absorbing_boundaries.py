#Basic setup for a laser pulse interation with a solid-density plasma layer 
#for results see fig. 6 in arXiv:2302.01893
import sys
import pipic
import matplotlib.pyplot as plt
import matplotlib.colors as plt_col
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from pipic.extensions import absorbing_boundaries
from pipic import consts,types
from numba import cfunc, carray, types as nbt

# arguments from command line

arg = sys.argv
fall = float(arg[1])
boundarySize_mult = float(arg[2])
angle = float(arg[3])

#===========================SIMULATION INITIALIZATION===========================
wavelength = 1e-4 #units? cm?
boundarySize = boundarySize_mult*wavelength
nx, XMin, XMax = 2**8, -16*wavelength, 16*wavelength
ny, YMin, YMax = 2**8 + int(2**8/(32*wavelength)*2*boundarySize), -16*wavelength-boundarySize, 16*wavelength+boundarySize
dx, dy = (XMax - XMin)/nx, (YMax - YMin)/ny
ts = 1/4 
timeStep = ts*dx/consts.light_velocity

#---------------------setting solver and simulation region----------------------
#sim = pipic.ec(nx=nx, ny=ny, XMin=XMin, XMax=XMax, YMin=YMin, YMax=YMax) 
sim = pipic.init(solver='fourier_boris',nx=nx,ny=ny,xmin=XMin,xmax=XMax,ymin=YMin,ymax=YMax) 
#---------------------------setting field of the pulse--------------------------
amplitude_a0 = 10.
omega = 2*np.pi*consts.light_velocity/wavelength
a0 = consts.electron_mass*consts.light_velocity*omega/(-consts.electron_charge)
fieldAmplitude = amplitude_a0*a0
PulseDuration = 10e-15 


k = 2*np.pi/wavelength
#fp = f'data_ts_0,25_fall_{str(fall)}_bs_4.txt'
fp = f'./data/data_ts_{ts}_fall_{str(fall)}_bs_{str(boundarySize_mult)}_angle_{str(angle)}.txt'
rot = np.pi*angle/180


@cfunc(types.field_loop_callback)
def setField_callback(ind, r, E, B, dataDouble, dataInt):
    if dataInt[0] == 0:

        r0_rot = r[0]*np.cos(rot) - r[1]*np.sin(rot)
        r1_rot = r[0]*np.sin(rot) + r[1]*np.cos(rot)
        r2_rot = r[2]

        x = r0_rot #- init_laser_pos
        rho2 = r1_rot**2 + r0_rot**2
        
        k = 2*np.pi/wavelength
        # Rayleigh length
        waist = PulseDuration*consts.light_velocity
        pulseWidth_x = PulseDuration*consts.light_velocity
        focusPosition = np.sqrt(YMax**2 + XMax**2) # distance from the focus to the center of the simulation box
        Zr = np.pi*waist**2/wavelength 
    
        # curvature
        R = focusPosition*(1+(Zr/focusPosition)**2)        
        spotsize_init = waist*np.sqrt(1+(focusPosition/Zr)**2)
        phase = np.arctan(focusPosition/Zr)    
        amp = fieldAmplitude*(waist/spotsize_init)*np.exp(-rho2/spotsize_init**2)
        curvature = np.exp(1j*k*rho2/(2*R))
        gp = np.real(amp*curvature*np.exp(-1j*(k*x + phase))*np.exp(-x**2/(2*pulseWidth_x**2)))

        # z-polarized 
        E[2] = gp
        B1 = -gp
        B0 = 0
        B[0] = B0*np.cos(-rot) - B1*np.sin(-rot)
        B[1] = B0*np.sin(-rot) + B1*np.cos(-rot)

        # y-polarized
        B[2] = gp
        E1 = gp
        E0 = 0
        E[0] = E0*np.cos(-rot) - E1*np.sin(-rot)
        E[1] = E0*np.sin(-rot) + E1*np.cos(-rot)

#=================================OUTPUT========================================

#-------------------------preparing output of fields (x)-----------------------------
Ez = np.zeros((nx,ny), dtype=np.double) 
By = np.zeros((nx,ny), dtype=np.double) 
Ey = np.zeros((nx,ny), dtype=np.double) 
Bz = np.zeros((nx,ny), dtype=np.double) 
Ex = np.zeros((nx,ny), dtype=np.double) 
Bx = np.zeros((nx,ny), dtype=np.double) 


#------------------get functions-----------------------------------------------

@cfunc(types.field_loop_callback)
def get_field_Bz(ind, r, E, B, dataDouble, dataInt):
    _Bz = carray(dataDouble, Bz.shape, dtype=np.double)
    _Bz[ind[0], ind[1]] = B[2]

@cfunc(types.field_loop_callback)
def get_field_Ez(ind, r, E, B, dataDouble, dataInt):
    _Ez = carray(dataDouble, Ez.shape, dtype=np.double)
    _Ez[ind[0], ind[1]] = E[2]

@cfunc(types.field_loop_callback)
def get_field_By(ind, r, E, B, dataDouble, dataInt):
    _By = carray(dataDouble, By.shape, dtype=np.double)
    _By[ind[0], ind[1]] = B[1]

@cfunc(types.field_loop_callback)
def get_field_Ey(ind, r, E, B, dataDouble, dataInt):
    _Ey = carray(dataDouble, Ey.shape, dtype=np.double)
    _Ey[ind[0], ind[1]] = E[1]

@cfunc(types.field_loop_callback)
def get_field_Ex(ind, r, E, B, dataDouble, dataInt):
    _Ex = carray(dataDouble, Ex.shape, dtype=np.double)
    _Ex[ind[0], ind[1]] = E[0]

@cfunc(types.field_loop_callback)
def get_field_Bx(ind, r, E, B, dataDouble, dataInt):
    _Bx = carray(dataDouble, Bx.shape, dtype=np.double)
    _Bx[ind[0], ind[1]] = B[0]


def load_fields():
    sim.field_loop(handler=get_field_Ey.address, data_double=pipic.addressof(Ey), use_omp=True)
    sim.field_loop(handler=get_field_Ez.address, data_double=pipic.addressof(Ez), use_omp=True)
    sim.field_loop(handler=get_field_By.address, data_double=pipic.addressof(By), use_omp=True)
    sim.field_loop(handler=get_field_Bz.address, data_double=pipic.addressof(Bz), use_omp=True)
    sim.field_loop(handler=get_field_Bx.address, data_double=pipic.addressof(Bx), use_omp=True)
    sim.field_loop(handler=get_field_Ex.address, data_double=pipic.addressof(Ex), use_omp=True)



def get_energy():
    return (Ex**2 + Ey**2 + Ez**2 + Bx**2 + By**2 + Bz**2)/2

#------------------set absorbing boundaries-----------------------------------
field_handler_adress = absorbing_boundaries.field_handler(simulation_box =sim.simulation_box(),
                                                          characteristic_wavelength = wavelength,
                                                          boundary_size = boundarySize,
                                                          axis = 'y',
                                                          fall = fall)

sim.add_handler(name=absorbing_boundaries.name,
                subject='fields', 
                field_handler=field_handler_adress,)

#------------------- plot fields -------------------------------------------
fig, axs = plt.subplots(1, constrained_layout=True)

E = get_energy()
E_plot = axs.imshow(E.T,vmin=0, vmax=fieldAmplitude**2, 
        extent=[XMin,XMax,YMin,YMax], interpolation='none', \
        aspect='equal', cmap='inferno_r', origin='lower')
fig.colorbar(E_plot,ax=axs)

#===============================SIMULATION======================================
dataInt = np.zeros((1, ), dtype=np.intc) # data for passing the iteration number
#s = round((dx/consts.light_velocity/4)/timeStep*3000) 

prop_length = (YMax-YMin)/np.sin(rot) # propagation length of the laser pulse in the simulation box
s = round(prop_length/consts.light_velocity/timeStep+10)
#s = round(2**5*np.sqrt(2)*wavelength/consts.light_velocity/timeStep + 10) # total number of time steps
u = np.empty((s//20+1,2))
count = 0


sim.fourier_solver_settings(divergence_cleaning=1, sin2_kfilter=-1)
sim.advance(time_step=0, number_of_iterations=1,use_omp=True)
sim.field_loop(handler=setField_callback.address, data_int=pipic.addressof(dataInt), \
                  use_omp=True)
sim.advance(time_step=0, number_of_iterations=1,use_omp=True)
sim.fourier_solver_settings(divergence_cleaning=0, sin2_kfilter=-1)

for i in range(s):
    dataInt[0] = i
    sim.advance(time_step=timeStep, number_of_iterations=1,use_omp=True)
    if i%20==0:
        load_fields()
        E = get_energy()
        E_plot.set_data(E.T)
        u_ = np.sum(E)*dx*dy
        u[count,:] = [i*timeStep,u_]
        print(f'{count,i}/{s}: ',f':{u_}')
        count += 1
        #fig.savefig('im' + str(i) + '.png')

np.savetxt(fp, u, header='time [s] energy')

