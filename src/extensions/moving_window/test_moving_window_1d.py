#Basic setup for a laser pulse interation with a solid-density plasma layer 
#for results see fig. 6 in arXiv:2302.01893
import sys
import pipic
from pipic.extensions import moving_window
from pipic import consts,types
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from numba import cfunc, carray, types as nbt


if __name__ == '__main__':


    #===========================SIMULATION INITIALIZATION===========================
    # Electron number density
    n0 = 8e18 #1e18 # [1/cm^3]
    # plasma frequency, [e] = statC, [me] = g, [n0] = 1/cm^3, [c] = cm/s
    omega_p = np.sqrt(np.pi*4*consts.electron_charge**2*n0/consts.electron_mass) # [1/s]
    # laser wavelength
    wl = 1e-4 # [cm]

    pulseWidth_x = 2*wl # [cm] (radial size of the laser focus)
    nx, xmin, xmax = 2**8, -10*wl, 10*wl 
    dx = (xmax - xmin)/nx
    # 10 timesteps per laser cycle
    timestep = dx/consts.light_velocity/2
    thickness = 2**5 # thickness (in dz) of the area where the density and field is restored/removed 


    #---------------------setting solver and simulation region----------------------
    sim=pipic.init(solver='ec',nx=nx,xmin=xmin,xmax=xmax)

    #=================================FIELD========================================
    # laser (radial) frequency
    omega = 2 * np.pi * consts.light_velocity / wl # [1/s]
    a0 = 5 # [unitless]
    # Field amplitude
    E0 = -a0 * consts.electron_mass * consts.light_velocity * omega / consts.electron_charge #[statV/cm] 


    @cfunc(types.field_loop_callback)
    def initiate_field_callback(ind, r, E, B, data_double, data_int):

        if data_int[0] == 0:       
            x = r[0]
            
            k = 2*np.pi/wl
            gp = np.real(E0*np.exp(-1j*(k*x))*np.exp(-x**2/(2*pulseWidth_x**2)))

            # z-polarized 
            E[2] = -gp       
            B[1] = gp        

    #=================================PLASMA PROFILE========================================
    # plasma profile
    debye_length = 1e-2 # [cm] 
    temperature = 0#1e-12 * 4 * np.pi * n0 * consts.electron_charge ** 2 * debye_length ** 2 # [erg/kB] (?)
    particles_per_cell = 5
    start_upramp = xmax
    end_upramp = start_upramp + xmax
    density_drop = end_upramp
    end_of_plasma = density_drop + 2*xmax

    @cfunc(types.add_particles_callback)
    def density_profile(r, data_double, data_int):
        # r is the position in the 'lab frame'  
        R = r[0]  # moving window position
            
        if R>start_upramp and R < end_upramp: 
            return n0*((R-start_upramp)/(end_upramp-start_upramp))
        elif R >= end_upramp and R < density_drop:
            return n0
        elif R >= density_drop and R < end_of_plasma:
            return n0*0.5
        else:
            return 0 

    #=================================OUTPUT========================================
    #-------------------------preparing output of fields (x)-----------------------------
    Ex = np.zeros((nx,), dtype=np.double) 
    Ey = np.zeros((nx,), dtype=np.double) 
    Ez = np.zeros((nx,), dtype=np.double) 
    rho = np.zeros((nx,), dtype=np.double) 


    #------------------get functions-----------------------------------------------
        
    @cfunc(types.particle_loop_callback)
    def get_density(r, p, w, id, data_double, data_int):  
        ix = int(rho.shape[0]*(r[0] - xmin)/(xmax - xmin))


        data = carray(data_double, rho.shape, dtype=np.double)
        
        if ix < rho.shape[0]:
            data[ix] += w[0]/(dx)
    

    @cfunc(types.field_loop_callback)
    def get_field_Ey(ind, r, E, B, data_double, data_int):
        _E = carray(data_double, Ey.shape, dtype=np.double)
        _E[ind[0],] = E[1]
 
    @cfunc(types.field_loop_callback)
    def get_field_Ex(ind, r, E, B, data_double, data_int):       
        _E = carray(data_double, Ex.shape, dtype=np.double)
        _E[ind[0],] = E[0]
  
    @cfunc(types.field_loop_callback)
    def get_field_Ez(ind, r, E, B, data_double, data_int):      
        _E = carray(data_double, Ez.shape, dtype=np.double)
        _E[ind[0],] = E[2]


    def load_fields():
        sim.field_loop(handler=get_field_Ey.address, data_double=pipic.addressof(Ey))
        sim.field_loop(handler=get_field_Ex.address, data_double=pipic.addressof(Ex))
        sim.field_loop(handler=get_field_Ez.address, data_double=pipic.addressof(Ez))

    #===============================SIMULATION======================================
    data_int = np.zeros((1, ), dtype=np.intc) # data for passing the iteration number

    #-----------------------initiate field and plasma-------------------------
    sim.field_loop(handler=initiate_field_callback.address, data_int=pipic.addressof(data_int),
                    use_omp=True)

    # This part is just for initiating the electron species, 
    # so that the algorithm knows that there is a species called electron
    # it is therefore not important where the electrons are or what density they have
    sim.add_particles(name='electron', number=1,#particles_per_cell),
                    charge=consts.electron_charge, mass=consts.electron_mass,
                    temperature=temperature, density=density_profile.address,)

    #-----------------------adding the handler of extension-------------------------
    window_speed = consts.light_velocity #speed of moving window
    density_handler_adress = moving_window.handler(sim.ensemble_data(),
                                                   sim.simulation_box(),
                                                   thickness=thickness,
                                                   particles_per_cell=particles_per_cell,
                                                   temperature=temperature,
                                                   density_profile=density_profile.address,
                                                   velocity=window_speed,)
    field_handler_adress = moving_window.field_handler(sim.simulation_box(),timestep)

    sim.add_handler(name=moving_window.name, 
                    subject='electron,cells',
                    handler=density_handler_adress,
                    field_handler=field_handler_adress,
                    data_int=pipic.addressof(data_int),)
    
    #-----------------------run simulation-------------------------
    s = int((end_of_plasma+xmax)/consts.light_velocity/timestep) # number of steps in the simulation 
    checkpoint = 10

    fig, ax = plt.subplots(2, 1, figsize=(15, 5))
    for i in range(s):
        sim.advance(time_step=timestep, number_of_iterations=1,use_omp=False)

        if i%checkpoint==0:
            print(f'Iteration {i}/{s}')
            rho.fill(0)

            # load fields and densities            
            sim.particle_loop(name='electron', handler=get_density.address,
                        data_double=pipic.addressof(rho))
            load_fields()

            # plot fields and densities
            ax[0].cla()
            ax[1].cla()
            im = ax[0].plot(rho, label='rho')
            ax[0].set_ylim(0,n0)
            ax[1].plot(Ez, label='Ex')
            plt.savefig(f'./rho_{i:04d}.png')

