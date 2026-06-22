#Basic setup for a laser pulse interation with a solid-density plasma layer 
#for results see fig. 6 in arXiv:2302.01893
import sys
import pipic
import numpy as np
from pipic.extensions import absorbing_boundaries
from pipic import consts,types
from numba import cfunc, carray


def get_pic_steps(default_steps):
    if "--steps" in sys.argv:
        idx = sys.argv.index("--steps")
        return max(1, int(sys.argv[idx + 1]))
    return default_steps




#===========================SIMULATION INITIALIZATION===========================
wavelength = 1e-4 #units? cm?
nx, XMin, XMax = 2**8, -16*wavelength, 16*wavelength
ny, YMin, YMax = 2**8, -16*wavelength, 16*wavelength
dx, dy = (XMax - XMin)/nx, (YMax - YMin)/ny
timeStep = dx/consts.light_velocity/4

#---------------------setting solver and simulation region----------------------
#sim = pipic.ec(nx=nx, ny=ny, XMin=XMin, XMax=XMax, YMin=YMin, YMax=YMax) 
sim = pipic.init(solver='ec2',nx=nx,ny=ny,xmin=XMin,xmax=XMax,ymin=YMin,ymax=YMax) 
#---------------------------setting field of the pulse--------------------------
amplitude_a0 = 10.
omega = 2*np.pi*consts.light_velocity/wavelength
a0 = consts.electron_mass*consts.light_velocity*omega/(-consts.electron_charge)
fieldAmplitude = amplitude_a0*a0
PulseDuration = 5e-15
k = 2*np.pi/wavelength

boundarySize = (4)*wavelength
fall = 1
rot = -np.pi/4

@cfunc(types.field_loop_callback)
def setField_callback(ind, r, E, B, dataDouble, dataInt):
    if dataInt[0] == 0:

        r_rot = (r[0]*np.cos(rot)+r[1]*np.sin(rot), -r[0]*np.sin(rot)+r[1]*np.cos(rot))
        

        # P-polarized
        Ep_rot = fieldAmplitude*np.exp(-(r_rot[0]**2+r_rot[1]**2)/(PulseDuration*consts.light_velocity)**2)*np.cos(k*r_rot[0])
        Bp_rot = fieldAmplitude*np.exp(-(r_rot[0]**2+r_rot[1]**2)/(PulseDuration*consts.light_velocity)**2)*np.cos(k*r_rot[0])

        # S-polarized
        Es_rot = fieldAmplitude*np.exp(-(r_rot[0]**2+r_rot[1]**2)/(PulseDuration*consts.light_velocity)**2)*np.cos(k*r_rot[0])
        Bs_rot = -fieldAmplitude*np.exp(-(r_rot[0]**2+r_rot[1]**2)/(PulseDuration*consts.light_velocity)**2)*np.cos(k*r_rot[0])


        E[0] = Ep_rot*np.sin(-rot)
        E[1] = Ep_rot*np.cos(rot) 
        B[2] = Bp_rot

        B[0] = Bs_rot*np.sin(-rot)
        B[1] = Bs_rot*np.cos(rot) 
        E[2] = Es_rot

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
    _Bz[ind[1], ind[0]] = B[2]

@cfunc(types.field_loop_callback)
def get_field_Ez(ind, r, E, B, dataDouble, dataInt):
    _Ez = carray(dataDouble, Ez.shape, dtype=np.double)
    _Ez[ind[1], ind[0]] = E[2]

@cfunc(types.field_loop_callback)
def get_field_By(ind, r, E, B, dataDouble, dataInt):
    _By = carray(dataDouble, By.shape, dtype=np.double)
    _By[ind[1], ind[0]] = B[1]

@cfunc(types.field_loop_callback)
def get_field_Ey(ind, r, E, B, dataDouble, dataInt):
    _Ey = carray(dataDouble, Ey.shape, dtype=np.double)
    _Ey[ind[1], ind[0]] = E[1]

@cfunc(types.field_loop_callback)
def get_field_Ex(ind, r, E, B, dataDouble, dataInt):
    _Ex = carray(dataDouble, Ex.shape, dtype=np.double)
    _Ex[ind[1], ind[0]] = E[0]

@cfunc(types.field_loop_callback)
def get_field_Bx(ind, r, E, B, dataDouble, dataInt):
    _Bx = carray(dataDouble, Bx.shape, dtype=np.double)
    _Bx[ind[1], ind[0]] = B[0]


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
field_handler_adress = absorbing_boundaries.field_handler(sim.simulation_box(), wavelength, boundarySize,'y',fall)

sim.add_handler(name=absorbing_boundaries.name,
                subject='fields', 
                field_handler=field_handler_adress,)

#===============================SIMULATION======================================
dataInt = np.zeros((1, ), dtype=np.intc) # data for passing the iteration number
s = 1000
s = get_pic_steps(s)
for i in range(s):
    dataInt[0] = i
    sim.field_loop(handler=setField_callback.address, data_int=pipic.addressof(dataInt), \
                  use_omp=True)
    sim.advance(time_step=timeStep, number_of_iterations=1,use_omp=True)
    if i%20==0:
        load_fields()
    energy = get_energy()
    print(i, f': energy_sum={np.sum(energy)*dx*dy}')

