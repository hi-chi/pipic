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
import h5py


def create_hdf5(fp, shape, dsets=['Ex','Ez','Ey','Bx','Bz','By','rho'],mode="w"):

    # Create a new HDF5 file
    file = h5py.File(fp, mode)

    # Create a dataset
    for d in dsets:
        file.create_dataset(d, shape=shape, dtype=np.double)
    
    file.close()


if __name__ == '__main__':
    s = 100 
    checkpoint = 10   

    #===========================SIMULATION INITIALIZATION===========================
    # Electron number density
    n0 = 8e18 #1e18 # [1/cm^3]
    # plasma frequency, [e] = statC, [me] = g, [n0] = 1/cm^3, [c] = cm/s
    omega_p = np.sqrt(np.pi*4*consts.electron_charge**2*n0/consts.electron_mass) # [1/s]
    # plasma wavelength
    wp = 2*np.pi*consts.light_velocity/omega_p # [cm]
    # laser wavelength
    wl = 1e-4 # [cm]
    # laser pulse width (spatial in x and radially in focus)
    pulseDuration_FWHM = 16e-15

    pulseWidth_x = (pulseDuration_FWHM/2.355)*consts.light_velocity # [cm]
    spotsize = 2*pulseWidth_x
    nx, xmin, xmax = 2**7, -5*spotsize, 5*spotsize 
    ny, ymin, ymax = 2**7, -5*spotsize, 5*spotsize
    nz, zmin, zmax = 2**7, -5*pulseWidth_x, 5*pulseWidth_x
    dx, dy, dz = (xmax - xmin)/nx, (ymax - ymin)/ny, (zmax - zmin)/nz
    # 10 timesteps per laser cycle
    timestep = 0.1/(2*np.pi*consts.light_velocity/wl) #0.5*dx_/consts.light_velocity
    thickness = 10 # thickness (in dz) of the area where the density and field is restored/removed 


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
    end_of_plasma = 0.8*simulation_length 

    @cfunc(types.add_particles_callback)
    def density_profile(r, data_double, data_int):
        # r is the position in the 'lab frame'  
        R = r[2] - (zmax - thickness*dz)
            
        if R > 0 and R < upramp: 
            return n0*(R/upramp)
        elif R > upramp and R < density_drop:
            return n0
        elif R > density_drop and R < end_of_plasma:
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
    window_speed = consts.light_velocity #speed of moving window

    #-----------------------adding the handler of extension-------------------------

    # variable for passing density to the handler  
    data_double = np.zeros((1, ), dtype=np.double)
    data_double[0] = n0*1e-10 # initiate with something low

    #-----------------------initiate field and plasma-------------------------

    
    sim.fourier_solver_settings(divergence_cleaning=1, sin2_kfilter=-1)
    sim.advance(time_step=0, number_of_iterations=1,use_omp=True)
    sim.field_loop(handler=initiate_field_callback.address, data_int=pipic.addressof(data_int),
                    use_omp=True)
    sim.advance(time_step=0, number_of_iterations=1,use_omp=True)
    sim.fourier_solver_settings(divergence_cleaning=0, sin2_kfilter=-1)

    # This part is just for initiating the electron species, 
    # so that the algorithm knows that there is a species called electron
    # it is therefore not important where the electrons are or what density they have
    sim.add_particles(name='electron', number=1,#particles_per_cell),
                    charge=consts.electron_charge, mass=consts.electron_mass,
                    temperature=temperature, density=density_profile.address,
                    data_int=pipic.addressof(data_int))
    
    density_handler_adress = moving_window.handler(sim.ensemble_data(),
                                                   thickness=thickness,
                                                   particles_per_cell=particles_per_cell,
                                                   temperature=temperature,
                                                   density=density_profile.address,)
    sim.add_handler(name=moving_window.name, 
                    subject='electron,cells',
                    handler=density_handler_adress,
                    data_int=pipic.addressof(data_int),
                    data_double=pipic.addressof(data_double))
    
    #-----------------------run simulation-------------------------
    
    dsets = ['Ex','Ey','Ez','rho','ps']
    fields = [Ex,Ey,Ez,rho,ps]
    ncp = s//checkpoint

    hdf5_fp = 'junk.h5'
    create_hdf5(hdf5_fp, shape=(ncp,nx,ny,nz),dsets=dsets[:-1])
    create_hdf5(hdf5_fp, shape=(ncp,nps,nz),dsets=['ps'],mode="r+")

    create_hdf5(hdf5_fp,shape=(nx,),dsets=['x_axis',],mode="r+")
    create_hdf5(hdf5_fp,shape=(ny,),dsets=['y_axis',],mode="r+")
    create_hdf5(hdf5_fp,shape=(nz,),dsets=['z_axis',],mode="r+")
    create_hdf5(hdf5_fp,shape=(nps,),dsets=['p_axis',],mode="r+")

    x_axis = np.linspace(xmin, xmax, nx)
    y_axis = np.linspace(ymin, ymax, ny)
    z_axis = np.linspace(zmin, zmax, nz)
    p_axis = np.linspace(pmin, pmax, nps)

    with h5py.File(hdf5_fp,"r+") as file:
        file['x_axis'][:] = x_axis
        file['y_axis'][:] = y_axis
        file['z_axis'][:] = z_axis
        file['p_axis'][:] = p_axis

    
    for i in range(s):
        
        data_int[0] = i 

        sim.advance(time_step=timestep, number_of_iterations=1,use_omp=True)
        
        sim.field_loop(handler=remove_field.address, data_int=pipic.addressof(data_int),
                    use_omp=True)
        
        

        if i%checkpoint==0:
            print(i)
            rho.fill(0)
            ps.fill(0)
            # load fields and densities            
            
            sim.particle_loop(name='electron', handler=get_density.address,
                        data_double=pipic.addressof(rho))
            sim.particle_loop(name='electron', handler=get_phase_space.address,
                        data_double=pipic.addressof(ps))
            load_fields()
            
            roll_back = np.round((i*timestep*window_speed)/dz,decimals=0).astype(int)
            
            try:
                with h5py.File(hdf5_fp,"r+") as file:
                    for j,f in enumerate(dsets):
                        file[f][i//checkpoint] = np.roll(fields[j],-roll_back,axis=-1) 
            except IOError:
                input('Another process is accessing the hdf5 file. Close the file and press enter.')
                with h5py.File(hdf5_fp,"r+") as file:
                    for j,f in enumerate(dsets):
                        file[f][i//checkpoint] = np.roll(fields[j],-roll_back,axis=-1)    

