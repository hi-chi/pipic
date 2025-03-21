import pipic
from pipic.tools import *
from pipic.extensions import x_converter_c
from numba import cfunc, carray
import matplotlib.pyplot as plt
import math, os, time
import numpy as np
from pipic import consts, types
import numba

#============ setting the simulation box ============
wavelength = 1e-4
L = 16*wavelength
sim = pipic.init(solver='ec',
                 nx=128, xmin=-L/2, xmax=L/2,
                 ny=128, ymin=-L/2, ymax=L/2,
                 nz=128, zmin=-L/2, zmax=L/2)

#====== setting the field of a focused pulse ========
from pipic.extensions import focused_pulse
focused_pulse.set_box(nx=sim.nx, ny=sim.ny, nz=sim.nz,
                      xmin=sim.xmin, ymin=sim.ymin, zmin=sim.zmin,
                      xmax=sim.xmax, ymax=sim.ymax, zmax=sim.xmax)
dist = 30*wavelength # distance to the initial location of the pulse (can be beyond the box)
starting_dist = 6*wavelength # the distance, which the pulse is brought to
focused_pulse.set_path(-dist/np.sqrt(2.0), -dist/np.sqrt(2.0), 0) # the path to go (defines the axis)
focused_pulse.set_e_axis(-1/np.sqrt(2.0), 1/np.sqrt(2.0), 0) # polarization axis (electric field)
focused_pulse.set_l_size(2*wavelength) # longitudinal size of the pulse
@cfunc(numba.types.double(numba.types.CPointer(numba.types.double)))
def pulse_shape(par): # initial field amplitude (at point r) as a function of the following parameters
    # par[0] = dist - abs(r), longitudinal coordinate relative to dist (\in [-l_size/2, l_size/2])
    # par[1] is the angle between r and -path, theta angle in polar coordinates (\in [0, pi))
    # par[2] is the angle between e_axis and r in the plane perpendicular to path (\in [0, 2*pi))
    # par[3] is the angle between r and e_axis (\in [0, pi))
    return math.sin(2*np.pi*par[0]/wavelength)*math.sin(par[3])*10*wavelength/dist
focused_pulse.set_shape(pulse_shape.address)
# now we do divergence cleaning and advance the field state to bring the pulse to the initial state of a simulation
sim.fourier_solver_settings(divergence_cleaning=True)
sim.advance(time_step=0)
time_start = time.time() # this is to set rho to zero for divergence cleaning
sim.field_loop(handler=focused_pulse.field_loop_cb()) # setting the field
print('focused_pulse took', time.time() - time_start, 's.')
sim.advance(time_step=(dist - starting_dist)/consts.light_velocity) # bringing the pulse to the initial location of a further simulation
sim.fourier_solver_settings(divergence_cleaning=False) # unless using fourier_boris, for which divergence cleaning is recommended

#============ show field in xy plane ==================
fig, ax = plt.subplots()
field_output = np.zeros((sim.ny, sim.nx), dtype=np.double)
xmin, ymin = sim.xmin, sim.ymin
nx, ny = sim.nx, sim.ny
stepx, stepy = L/nx, L/ny
@cfunc(types.it2r_callback)
def fieldPlot_it2r(it, r, data_double, data_int):
    r[0] = xmin + stepx*(it[0] % nx)
    r[1] = ymin + stepy*(it[0] - (it[0] % nx))/nx
    r[2] = 0
@cfunc(types.field_loop_callback)
def fieldPlot_cb(it, r, E, B, data_double, data_int):
    Field = carray(data_double, field_output.shape, dtype=np.double)
    Field[int((it[0] - (it[0] % nx))/nx), (it[0] % nx)] = np.sqrt(E[0]*E[0] + E[1]*E[1] + E[2]*E[2] + B[0]*B[0] + B[1]*B[1] + B[2]*B[2])
sim.custom_field_loop(number_of_iterations=nx*ny, it2r=fieldPlot_it2r.address,
                      field2data=fieldPlot_cb.address, data_double=pipic.addressof(field_output))
plot0 = ax.imshow(field_output, vmin=0, vmax=3,
                  extent=[sim.xmin,sim.xmax,sim.ymin,sim.ymax], interpolation='none',
                  aspect='equal', cmap='YlOrBr', origin='lower')
ax.set(xlabel='$x$ (cm)', ylabel='$y$ (cm)')
ax.ticklabel_format(axis='both', scilimits=(0,0), useMathText=True)

# =============== simulation =====================
time_step = 0.125*wavelength/consts.light_velocity
frames = 100
output_folder = 'test_focused_pulse_output'
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
for i in range(frames):
    sim.advance(time_step=time_step)
    sim.custom_field_loop(number_of_iterations=nx*ny, it2r=fieldPlot_it2r.address,
                      field2data=fieldPlot_cb.address, data_double=pipic.addressof(field_output))
    plot0.set_data(field_output)
    fig.savefig(output_folder + '/im' + str(i) + '.png')
    print(i)


# An example of setting field directly in Python
# @cfunc(types.field_loop_callback)
# def setField_cb(ind, r, E, B, data_double, data_int):
#     R = np.sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2])
#     if R == 0:
#        return
#     r0 = types.Double3(r[0], r[1], r[2])
#     r0.normalize()
#     e0 = types.Double3(-np.sqrt(1/2.0), np.sqrt(1/2.0), 0)
#     b0 = r0.cross(e0)
#     e0 = r0.cross(b0)
#     dist = 6*wavelength
#     lf = np.exp(-0.5*((R - dist)/wavelength)**2)*np.sin(2*np.pi*(R - dist)/wavelength)
#     E[0], E[1], E[2] = lf*e0.x, lf*e0.y, lf*e0.z
#     B[0], B[1], B[2] = lf*b0.x, lf*b0.y, lf*b0.z
# time_start = time.time()
# sim.field_loop(handler=setField_cb.address)
# print('focused_pulse py took', time.time() - time_start, 's.')
