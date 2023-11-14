import pipic
from pipic.tools import *
import matplotlib.pyplot as plt
import numpy as np


# ===========================SIMULATION INITIALIZATION===========================
temperature = 1e-6 * electron_mass * light_velocity**2
density = 1e+18
debye_length = sqrt(temperature / (4*pi * density * electron_charge**2))
plasma_period = sqrt(pi * electron_mass / (density * electron_charge**2))
l = 128*debye_length
xmin, xmax = -l/2, l/2
field_amplitude = 0.01*4*pi * (xmax-xmin) * electron_charge * density
nx = 128
time_step = plasma_period/64

# ---------------------setting solver and simulation region----------------------
sim = pipic.init(solver='ec', nx=nx, xmin=xmin, xmax=xmax)

# ------------------------------adding electrons---------------------------------
@cfunc(type_addParticles)
def density_callback(r, dataDouble, dataInt):
    return density * (abs(r[0]) < l/4)

sim.add_particles(name='electron', number=500*nx,
                  charge=-electron_charge, mass=electron_mass,
                  temperature=temperature, density=density_callback.address)

# ---------------------------setting initial field-------------------------------
@cfunc(type_fieldLoop)
def setField_callback(ind, r, E, B, dataDouble, dataInt):
    E[0] = field_amplitude * sin(4*pi * r[0] / (xmax-xmin)) * (abs(r[0]) < l/4)

sim.field_loop(handler=setField_callback.address)


# =================================OUTPUT========================================
fig, axs = plt.subplots(2, constrained_layout=True)

# -------------preparing output for electron distribution f(x, px)--------------
xpx_dist = numpy.zeros((64, 128), dtype=numpy.double)
pxLim = 5 * sqrt(temperature * electron_mass)
inv_dx_dpx = (xpx_dist.shape[1] / (xmax-xmin)) * (xpx_dist.shape[0] / (2 * pxLim))

@cfunc(type_particleLoop)
def xpx_callback(r, p, w, id, dataDouble, dataInt):
    ix = int(xpx_dist.shape[1] * (r[0] - xmin) / (xmax-xmin))
    iy = int(xpx_dist.shape[0] * 0.5 * (1 + p[0] / pxLim))
    data = carray(dataDouble, xpx_dist.shape, dtype=numpy.double)
    if iy >= 0 and iy < xpx_dist.shape[0]:
        data[iy, ix] += w[0] * inv_dx_dpx / (3*density/pxLim)

axs[0].set_title('$\partial N / \partial x \partial p_x$ (s g$^{-1}$cm$^{-2}$)')
axs[0].set(ylabel='$p_x$ (cm g/s)')
axs[0].xaxis.set_ticklabels([])
plot0 = axs[0].imshow(xpx_dist, vmin=0, vmax=1,
                      extent=[xmin, xmax, -pxLim, pxLim], interpolation='none',
                      aspect='auto', cmap='YlOrBr')
fig.colorbar(plot0, ax=axs[0], location='right')

def plot_xpx():
    xpx_dist.fill(0)
    sim.particle_loop(name='electron', handler=xpx_callback.address,
                      dataDouble=address_of(xpx_dist))
    plot0.set_data(xpx_dist)

# -------------------------preparing output of Ex(x)-----------------------------
Ex = numpy.zeros((32,), dtype=numpy.double)

@cfunc(type_it2r)
def Ex_it2r(it, r, dataDouble, dataInt):
    r[0] = xmin + (it[0] + 0.5) * (xmax-xmin) / Ex.shape[0]

@cfunc(type_field2data)
def get_Ex(it, r, E, B, dataDouble, dataInt):
    dataDouble[it[0]] = E[0]

axs[1].set_xlim([xmin, xmax])
axs[1].set_ylim([-field_amplitude, field_amplitude])
axs[1].set(xlabel='$x$ (cm)', ylabel='$E_x$ (cgs units)')
x_axis = np.linspace(xmin, xmax, Ex.shape[0])
plot_Ex_, = axs[1].plot(x_axis, Ex)

def plot_Ex():
    sim.custom_field_loop(number_of_iterations=Ex.shape[0], it2r=Ex_it2r.address,
                          field2data=get_Ex.address, dataDouble=address_of(Ex))
    plot_Ex_.set_ydata(Ex)


# ===============================SIMULATION======================================
output_folder = 'basic_example_output'
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
time_start = time.time()
for i in range(32):
    sim.advance(time_step=time_step, number_of_iterations=8)
    plot_xpx()
    plot_Ex()
    fig.savefig(output_folder + '/im' + str(i) + '.png')
    print(i, '/', 32)
