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
    s = 1500 
    checkpoint = 100

    #===========================SIMULATION INITIALIZATION===========================
    # Electron number density
    n0 = 8e18 #1e18 # [1/cm^3]
    # plasma frequency, [e] = statC, [me] = g, [n0] = 1/cm^3, [c] = cm/s
    omega_p = np.sqrt(np.pi*4*consts.electron_charge**2*n0/consts.electron_mass) # [1/s]
    # plasma wavelength
    wp = 2*np.pi*consts.light_velocity/omega_p # [cm]
    # laser wavelength
    wl = 1e-4 # [cm]

    pulseWidth_x = 2*wl # [cm] (radial size of the laser focus)
    spotsize = 100*wl # [cm] (radial size of the laser focus)
    nx, xmin, xmax = 2**4, -5*wl, 5*wl 
    ny, ymin, ymax = 2**4, -5*wl, 5*wl
    nz, zmin, zmax = 2**8, -10*wl, 10*wl
    dx, dy, dz = (xmax - xmin)/nx, (ymax - ymin)/ny, (zmax - zmin)/nz
    # 10 timesteps per laser cycle
    timestep = dz/consts.light_velocity/2
    thickness = 2**5 # thickness (in dz) of the area where the density and field is restored/removed 


    #---------------------setting solver and simulation region----------------------
    sim=pipic.init(solver='ec',nx=nx,ny=ny,nz=nz,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,zmin=zmin,zmax=zmax)

    #---------------------------setting field of the pulse--------------------------
    # laser (radial) frequency
    omega = 2 * np.pi * consts.light_velocity / wl # [1/s]
    a0 = 5 # [unitless]
    # Field amplitude
    E0 = -a0 * consts.electron_mass * consts.light_velocity * omega / consts.electron_charge #[statV/cm] 

    focusPosition = wp*2

    k = 2*np.pi/wl
    Zr = np.pi*omega**1/wl 
    R = focusPosition*(1+(Zr/focusPosition)**2)
    spotsize_init = spotsize*np.sqrt(1+(focusPosition/Zr)**2)

    # laser wave number
    k = 2*np.pi / wl # [1/cm]
    # Critical density
    N_cr = consts.electron_mass * omega ** 2 / (4 * np.pi * consts.electron_charge ** 2)

    debye_length = 1e-2 # [cm] 
    temperature = 0 #4 * np.pi * n0 * consts.electron_charge ** 2 * debye_length ** 2 # [erg/kB] (?)

    particles_per_cell = 5

    @cfunc(types.field_loop_callback)
    def initiate_field_callback(ind, r, E, B, data_double, data_int):

        if data_int[0] == 0:       
            x = r[2]
            rho2 = r[1]**2 + r[0]**2
            
            k = 2*np.pi/wl
            # Rayleigh length
            Zr = np.pi*spotsize**2/wl 
            # curvature
            R = focusPosition*(1+(Zr/focusPosition)**2)
            
            spotsize_init = spotsize*np.sqrt(1+(focusPosition/Zr)**2)
            phase = 0 #np.arctan(focusPosition/Zr)    
            amp = E0*(spotsize/spotsize_init)*np.exp(-rho2/spotsize_init**2)
            curvature = 1 #np.exp(1j*k*rho2/(2*R))

            gp = np.real(amp*curvature*np.exp(-1j*(k*x + phase))*np.exp(-x**2/(2*pulseWidth_x**2)))

            # x-polarized 
            E[0] = gp       
            B[1] = gp        

    # plasma profile
    simulation_length = s*timestep*consts.light_velocity
    upramp = 0.5*simulation_length
    density_drop = 0.5*simulation_length
    end_of_plasma = 0.9*simulation_length 

    @cfunc(types.add_particles_callback)
    def density_profile(r, data_double, data_int):
        # r is the position in the 'lab frame'  
        R = r[2] - (zmax - thickness*dz)
            
        if R > 0 and R < upramp: 
            return n0*(R/upramp)
        elif R > upramp and R < density_drop:
            return n0
        elif R >= density_drop and R < end_of_plasma:
            return n0*0.5
        else:
            return 0 


    @cfunc(types.field_loop_callback)
    def remove_field(ind, r, E, B, data_double, data_int):

        rollback = np.floor(data_int[0]*timestep*consts.light_velocity/dz)
        if rollback%(thickness//2)==0:
            r_rel = zmin + dz*(rollback%nz)  
            r_min = r_rel - thickness*dz
            r_max = r_rel 
            if (r[2] > r_min and r[2] < r_max) or (r[2] > zmax - (zmin - r_min)) or (r[2] < zmin + (r_max - zmax)): 
                E[1] = 0
                B[2] = 0
                E[2] = 0 
                B[1] = 0
                E[0] = 0
                B[0] = 0 

    #=================================OUTPUT========================================
    #-------------------------preparing output of fields (x)-----------------------------
    Ex = np.zeros((nx,ny,nz), dtype=np.double) 
    Ey = np.zeros((nx,ny,nz), dtype=np.double) 
    Ez = np.zeros((nx,ny,nz), dtype=np.double) 

    rho = np.zeros((nx,ny,nz), dtype=np.double) 
    pmin = 0
    pmax = 1e-15
    nps = 2**8
    ps = np.zeros((nps,nz), dtype=np.double) 

    #------------------get functions-----------------------------------------------
        
    @cfunc(types.particle_loop_callback)
    def get_density(r, p, w, id, data_double, data_int):  
        ix = int(rho.shape[0]*(r[0] - xmin)/(xmax - xmin))
        iy = int(rho.shape[1]*(r[1] - ymin)/(ymax - ymin))
        iz = int(rho.shape[2]*(r[2] - zmin)/(zmax - zmin))

        data = carray(data_double, rho.shape, dtype=np.double)
        
        if (iy < rho.shape[1] and 
            ix < rho.shape[0] and
            iz < rho.shape[2]):
            data[ix, iy, iz] += w[0]/(dx*dy*dz)
    
    @cfunc(types.particle_loop_callback)
    def get_phase_space(r, p, w, id, data_double, data_int):   
        iz = int(ps.shape[0]*(r[2] - zmin)/(zmax - zmin))
        ip = int(ps.shape[1]*(p[2] - pmin)/(pmax - pmin))
        data = carray(data_double, ps.shape, dtype=np.double)
        
        if ip>=0 and ip < ps.shape[1] and iz < ps.shape[0]:
            data[ip,iz] += w[0]/(dx*dy*dz) 


    @cfunc(types.field_loop_callback)
    def get_field_Ey(ind, r, E, B, data_double, data_int):
        _E = carray(data_double, Ey.shape, dtype=np.double)
        _E[ind[0], ind[1], ind[2]] = E[1]
 
    @cfunc(types.field_loop_callback)
    def get_field_Ex(ind, r, E, B, data_double, data_int):       
        _E = carray(data_double, Ex.shape, dtype=np.double)
        _E[ind[0], ind[1], ind[2]] = E[0]
  
    @cfunc(types.field_loop_callback)
    def get_field_Ez(ind, r, E, B, data_double, data_int):      
        _E = carray(data_double, Ez.shape, dtype=np.double)
        _E[ind[0], ind[1], ind[2]] = E[2]


    def load_fields():
        sim.field_loop(handler=get_field_Ey.address, data_double=pipic.addressof(Ey))
        sim.field_loop(handler=get_field_Ex.address, data_double=pipic.addressof(Ex))
        sim.field_loop(handler=get_field_Ez.address, data_double=pipic.addressof(Ez))

    #===============================SIMULATION======================================
    data_int = np.zeros((1, ), dtype=np.intc) # data for passing the iteration number

    #-----------------------initiate field and plasma-------------------------
    
    sim.fourier_solver_settings(divergence_cleaning=1, sin2_kfilter=-1)
    sim.advance(time_step=0, number_of_iterations=1,use_omp=True)
    sim.field_loop(handler=initiate_field_callback.address, data_int=pipic.addressof(data_int),
                    use_omp=True)
    sim.advance(time_step=0, number_of_iterations=1,use_omp=True)
    sim.fourier_solver_settings(divergence_cleaning=0, sin2_kfilter=-1)
    

    x_axis = np.linspace(xmin, xmax, nx)
    y_axis = np.linspace(ymin, ymax, ny)
    z_axis = np.linspace(zmin, zmax, nz)
    p_axis = np.linspace(pmin, pmax, nps)
    # This part is just for initiating the electron species, 
    # so that the algorithm knows that there is a species called electron
    # it is therefore not important where the electrons are or what density they have
    sim.add_particles(name='electron', number=1,#particles_per_cell),
                    charge=consts.electron_charge, mass=consts.electron_mass,
                    temperature=temperature, density=density_profile.address,)
    #-----------------------adding the handler of extension-------------------------
    data_double = np.zeros((1, ), dtype=np.double) # data for passing some other data
    window_speed = consts.light_velocity #speed of moving window
    density_handler_adress = moving_window.handler(sim.ensemble_data(),
                                                   thickness=thickness,
                                                   particles_per_cell=particles_per_cell,
                                                   temperature=temperature,
                                                   density=density_profile.address,
                                                   velocity=window_speed,)
    field_handler_adress = moving_window.field_handler(sim.simulation_box(),timestep)
    sim.add_handler(name=moving_window.name, 
                    subject='electron,cells',
                    handler=density_handler_adress,
                    field_handler=field_handler_adress,
                    data_int=pipic.addressof(data_int),)
    
    #-----------------------run simulation-------------------------
    fig, ax = plt.subplots(2, 1, figsize=(15, 5))
    for i in range(s):
        data_int[0] = i 
        sim.advance(time_step=timestep, number_of_iterations=1,use_omp=False)

        if i%checkpoint==0:
            print(f'Iteration {i}')
            rho.fill(0)
            ps.fill(0)
            # load fields and densities            
            
            sim.particle_loop(name='electron', handler=get_density.address,
                        data_double=pipic.addressof(rho))
            sim.particle_loop(name='electron', handler=get_phase_space.address,
                        data_double=pipic.addressof(ps))
            load_fields()
            
            roll_back = np.round((i*timestep*window_speed)/dz,decimals=0).astype(int)

            im = ax[0].imshow(rho[nx//2, :, :], origin='lower', aspect='auto',
                        extent=[ymin, ymax, zmin, zmax], cmap='Reds',vmin=0,vmax=n0)
            ax[1].imshow(Ex[nx//2, :, :], origin='lower', aspect='auto',
                        extent=[ymin, ymax, zmin, zmax], cmap='seismic')
            plt.savefig(f'./rho_{i:04d}.png')

