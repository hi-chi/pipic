#Basic setup for a laser pulse interation with a solid-density plasma layer 
#for results see fig. 6 in arXiv:2302.01893
import sys
import pipic
import numpy as np
from pipic.extensions import moving_window
from pipic import consts,types
from numba import cfunc, carray


def get_pic_steps(default_steps):
    if "--steps" in sys.argv:
        idx = sys.argv.index("--steps")
        return max(1, int(sys.argv[idx + 1]))
    return default_steps


#===========================SIMULATION INITIALIZATION===========================
wavelength = 1e-4 #units? cm?
nx, XMin, XMax = 2**8, -24*wavelength, 24*wavelength
ny, YMin, YMax = 2**8, -24*wavelength, 24*wavelength
dx, dy = (XMax - XMin)/nx, (YMax - YMin)/ny
timeStep = dx/consts.light_velocity/4

#---------------------setting solver and simulation region----------------------
#sim = pipic.ec(nx=nx, ny=ny, XMin=XMin, XMax=XMax, YMin=YMin, YMax=YMax) 
sim = pipic.init(solver='ec',nx=nx,ny=ny,xmin=XMin,xmax=XMax,ymin=YMin,ymax=YMax) 
#---------------------------setting field of the pulse--------------------------
density = 1e18 # in cm^-3
pulseWidth_x = wavelength*2#*0.001 # [cm] (radial size of the laser focus)
a0 = 10.
omega = 2*np.pi*consts.light_velocity/wavelength
E0 = -a0 * consts.electron_mass * consts.light_velocity * omega / consts.electron_charge #[statV/cm] 
initial_pos = XMax - XMax/2
focusPosition = XMax

@cfunc(types.field_loop_callback)
def initiate_field_callback(ind, r, E, B, data_double, data_int):
    if data_int[0] == 0:  
        x = r[0] - initial_pos     
        rho = x**2+r[1]**2 #+ initial_pos

        # Rayleigh length
        Zr = np.pi*pulseWidth_x**2/wavelength 
    
        # curvature
        fp_rel = focusPosition-initial_pos
        R = fp_rel*(1+(Zr/fp_rel)**2)        
        spotsize_init = pulseWidth_x*np.sqrt(1+(fp_rel/Zr)**2)
        phase = np.arctan(fp_rel/Zr)    
        amp = E0*(pulseWidth_x/spotsize_init)*np.exp(-rho/spotsize_init**2)
        k = 2*np.pi/wavelength
        curvature = np.exp(1j*k*rho/(2*R))
        gp = np.real(amp*curvature*np.exp(-1j*(k*x + phase))*np.exp(-x**2/(2*pulseWidth_x**2)))
        # z-polarized 
        E[2] = gp       
        B[1] = -gp   


start_plasma = XMax
@cfunc(types.add_particles_callback)
def density_callback_electron(r, data_double, data_int):  # callback function
    if r[0] > start_plasma:
        return density
    else:
        return 0

#=================================OUTPUT========================================

#-------------------------preparing output of fields (x)-----------------------------
Ez = np.zeros((nx,ny), dtype=np.double) 
By = np.zeros((nx,ny), dtype=np.double) 
Ey = np.zeros((nx,ny), dtype=np.double) 
Bz = np.zeros((nx,ny), dtype=np.double) 
Ex = np.zeros((nx,ny), dtype=np.double) 
Bx = np.zeros((nx,ny), dtype=np.double) 
rho = np.zeros((nx,ny), dtype=np.double) # density

#------------------get functions-----------------------------------------------
@cfunc(types.particle_loop_callback)
def get_density(r, p, w, id, data_double, data_int):   
    ix = int(nx*(r[0] - XMin)/(XMax - XMin))
    iy = int(ny*(r[1] - YMin)/(YMax - YMin))
    data = carray(data_double, rho.shape, dtype=np.double)
    
    if ix < rho.shape[0] and iy < rho.shape[1]:
        data[ix,iy] += w[0]/(dx*dy)


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
    sim.particle_loop(name='electron', handler=get_density.address, data_double=pipic.addressof(rho))


def get_energy():
    return (Ex**2 + Ey**2 + Ez**2 + Bx**2 + By**2 + Bz**2)/2

#---------------------add particles---------------------------------------

sim.add_particles(
    name="electron",
    number=nx*ny*5,
    charge= -consts.electron_charge,
    mass= consts.electron_mass,
    temperature=1e-12,
    density=density_callback_electron.address,
)

#-----------------------initiate field-------------------------
sim.field_loop(handler=initiate_field_callback.address, data_int=pipic.addressof(np.zeros((1,), dtype=np.intc)),
use_omp=True)

#------------------set absorbing boundaries-----------------------------------

dataInt = np.zeros((1, ), dtype=np.intc) # data for passing the iteration number

#-------------------moving_window--------------------------------------


field_handler_adress = moving_window.field_handler(sim.simulation_box(),
                                                   thickness = 16)
particle_handler_adress = moving_window.handler(sim.simulation_box(),
                                                thickness=16,
                                                particles_per_cell=5,
                                                temperature=1e-12,
                                                density_profile=density_callback_electron.address,)
sim.add_handler(name=moving_window.name,
                subject='cells,electron',
                field_handler=field_handler_adress,
                handler=particle_handler_adress,
                data_int=pipic.addressof(dataInt),)



#===============================SIMULATION======================================
default_steps = 2000
s = get_pic_steps(default_steps)
for i in range(s):
    dataInt[0] = i
    sim.advance(time_step=timeStep, number_of_iterations=1,use_omp=True)
    if i%20==0:
        rho.fill(0)
        load_fields()
        energy = get_energy()
        print(i, f': energy_sum={np.sum(energy)*dx*dy}, rho_sum={np.sum(rho)*dx*dy}')

