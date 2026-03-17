# This is the file used for producing fig. 2 in arXiv:2302.01893.
import pipic
from pipic import consts, types
import numpy as np
from numba import cfunc, carray
import sys


def get_pic_steps(default_steps):
    if "--steps" in sys.argv:
        idx = sys.argv.index("--steps")
        return max(1, int(sys.argv[idx + 1]))
    return default_steps


# ===========================SIMULATION INITIALIZATION===========================
temperature = 1e-6 * consts.electron_mass * consts.light_velocity**2
density = 1e18
debye_length = np.sqrt(temperature / (4 * np.pi * density * consts.electron_charge**2))
plasma_period = np.sqrt(np.pi * consts.electron_mass / (density * consts.electron_charge**2))
xmin, xmax = -64 * debye_length, 64 * debye_length
field_amplitude = 0.01 * 4 * np.pi * (xmax - xmin) * consts.electron_charge * density
nx = 128  # number of cells
time_step = plasma_period / 64

# ---------------------setting solver and simulation region----------------------
sim = pipic.init(solver="ec", nx=nx, xmin=xmin, xmax=xmax)
# by default ny=nz=1, YMax=ZMax=0.5, YMin=ZMin=-0.5


# ------------------------------adding electrons---------------------------------
@cfunc(types.add_particles_callback)
def density_callback(r, data_double, data_int):  # callback function
    return density  # can be any function of coordinate r[0]


sim.add_particles(
    name="electron",
    number=nx * 1000,
    charge=-consts.electron_charge,
    mass=consts.electron_mass,
    temperature=temperature,
    density=density_callback.address,
)


# ---------------------------setting initial field-------------------------------
@cfunc(types.field_loop_callback)
def setField_callback(ind, r, E, B, data_double, data_int):
    E[0] = field_amplitude * np.sin(2 * np.pi * r[0] / (xmax - xmin))


sim.field_loop(handler=setField_callback.address)


# =================================OUTPUT========================================
# -------------preparing output for electron distrribution f(x, px)--------------
xpx_dist = np.zeros((64, 128), dtype=np.double)
pxLim = 5 * np.sqrt(temperature * consts.electron_mass)
inv_dx_dpx = (xpx_dist.shape[1] / (xmax - xmin)) * (xpx_dist.shape[0] / (2 * pxLim))


@cfunc(types.particle_loop_callback)
def xpx_callback(r, p, w, id, data_double, data_int):
    ix = int(xpx_dist.shape[1] * (r[0] - xmin) / (xmax - xmin))
    iy = int(xpx_dist.shape[0] * 0.5 * (1 + p[0] / pxLim))
    data = carray(data_double, xpx_dist.shape, dtype=np.double)
    if 0 <= iy < xpx_dist.shape[0]:
        data[iy, ix] += w[0] * inv_dx_dpx


def plot_xpx():
    xpx_dist.fill(0)
    sim.particle_loop(
        name="electron", handler=xpx_callback.address, data_double=pipic.addressof(xpx_dist)
    )
    return np.sum(xpx_dist)


# -------------------------preparing output of Ex(x)-----------------------------
Ex = np.zeros((32,), dtype=np.double)


@cfunc(types.it2r_callback)
def Ex_it2r(it, r, data_double, data_int):
    r[0] = xmin + (it[0] + 0.5) * (xmax - xmin) / Ex.shape[0]


@cfunc(types.field2data_callback)
def get_Ex(it, r, E, B, data_double, data_int):
    data_double[it[0]] = E[0]


def plot_Ex():
    sim.custom_field_loop(
        number_of_iterations=Ex.shape[0],
        it2r=Ex_it2r.address,
        field2data=get_Ex.address,
        data_double=pipic.addressof(Ex),
    )
    return np.sum(Ex)


# ===============================SIMULATION======================================
steps = get_pic_steps(64)
for i in range(steps):
    sim.advance(time_step=time_step, number_of_iterations=1)
    xpx_sum = plot_xpx()
    ex_sum = plot_Ex()
    print(i, "/", steps, "sum(xpx)=", xpx_sum, "sum(Ex)=", ex_sum)
