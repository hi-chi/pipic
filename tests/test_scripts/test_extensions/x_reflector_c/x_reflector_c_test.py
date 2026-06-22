import pipic
from pipic.tools import *
from pipic.extensions import x_reflector_c
from numba import cfunc, carray
import math, numpy, os, time
import sys


def get_pic_steps(default_steps):
    if "--steps" in sys.argv:
        idx = sys.argv.index("--steps")
        return max(1, int(sys.argv[idx + 1]))
    return default_steps


# ===========================SIMULATION INITIALIZATION===========================
temperature = 1e-6 * electron_mass * light_velocity**2
density = 1e18
debye_length = math.sqrt(temperature / (4 * math.pi * density * electron_charge**2))
plasma_period = math.sqrt(math.pi * electron_mass / (density * electron_charge**2))
L = 128 * debye_length
xmin, xmax = -L / 2, L / 2
field_amplitude = -0.01 * 4 * math.pi * (xmax - xmin) * electron_charge * density
nx = 128  # number of cells
time_step = plasma_period / 64

# ---------------------setting solver and simulation region----------------------
sim = pipic.init(solver="ec", nx=nx, xmin=xmin, xmax=xmax)

# -----------------------adding the handler of extension-------------------------
sim.add_handler(
    name=x_reflector_c.name,
    subject="electron",
    handler=x_reflector_c.handler(location=L / 4 + L / 32, thickness=L / 16),
)


# ------------------------------adding electrons---------------------------------
@cfunc(add_particles_callback)
def density_callback(r, data_double, data_int):
    return density * (abs(r[0]) < L / 4)


sim.add_particles(
    name="electron",
    number=nx * 500,
    charge=electron_charge,
    mass=electron_mass,
    temperature=temperature,
    density=density_callback.address,
)


# ---------------------------setting initial field-------------------------------
@cfunc(field_loop_callback)
def setField_callback(ind, r, E, B, data_double, data_int):
    E[0] = field_amplitude * math.sin(4 * math.pi * r[0] / (xmax - xmin)) * (abs(r[0]) < L / 4)


sim.field_loop(handler=setField_callback.address)


# =================================OUTPUT========================================

# -------------preparing output for electron distrribution f(x, px)--------------
xpx_dist = numpy.zeros((64, 128), dtype=numpy.double)
pxLim = 5 * math.sqrt(temperature * electron_mass)
inv_dx_dpx = (xpx_dist.shape[1] / (xmax - xmin)) * (xpx_dist.shape[0] / (2 * pxLim))


@cfunc(particle_loop_callback)
def xpx_callback(r, p, w, id, data_double, data_int):
    ix = int(xpx_dist.shape[1] * (r[0] - xmin) / (xmax - xmin))
    iy = int(xpx_dist.shape[0] * 0.5 * (1 + p[0] / pxLim))
    data = carray(data_double, xpx_dist.shape, dtype=numpy.double)
    if iy >= 0 and iy < xpx_dist.shape[0]:
        data[iy, ix] += w[0] * inv_dx_dpx / (3 * density / pxLim)


def plot_xpx():
    xpx_dist.fill(0)
    sim.particle_loop(
        name="electron", handler=xpx_callback.address, data_double=addressof(xpx_dist)
    )


# -------------------------preparing output of Ex(x)-----------------------------
Ex = numpy.zeros((32,), dtype=numpy.double)


@cfunc(it2r_callback)
def Ex_it2r(it, r, data_double, data_int):
    r[0] = xmin + (it[0] + 0.5) * (xmax - xmin) / Ex.shape[0]


@cfunc(field2data_callback)
def get_Ex(it, r, E, B, data_double, data_int):
    data_double[it[0]] = E[0]


def plot_Ex():
    sim.custom_field_loop(
        number_of_iterations=Ex.shape[0],
        it2r=Ex_it2r.address,
        field2data=get_Ex.address,
        data_double=addressof(Ex),
    )


# ===============================SIMULATION======================================

i_max = get_pic_steps(32)
time_start = time.time()
for i in range(i_max):
    sim.advance(time_step=time_step, number_of_iterations=8)
    plot_xpx()
    plot_Ex()
    print(Ex.sum(), xpx_dist.sum())
print("Total time of simulation and output is", time.time() - time_start, "s.")
print("For the time taken by components of pi-PIC see pipic_performance.txt.")
