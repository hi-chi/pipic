import pipic
from pipic import consts, types
import numpy as np
from numba import cfunc, carray
import os, time
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
l = 128 * debye_length
xmin, xmax = -l / 2, l / 2
field_amplitude = 0.01 * 4 * np.pi * (xmax - xmin) * consts.electron_charge * density
nx = 128
time_step = plasma_period / 64

# ---------------------setting solver and simulation region----------------------
sim = pipic.init(solver="electrostatic_1d", nx=nx, xmin=xmin, xmax=xmax)

# ------------------------------adding electrons---------------------------------
@cfunc(types.add_particles_callback)
def density_callback(r, data_double, data_int):
    return density #* (abs(r[0]) < l / 4)


sim.add_particles(
    name="electron",
    number=500 * nx,
    charge=consts.electron_charge,
    mass=consts.electron_mass,
    temperature=temperature,
    density=density_callback.address,
)


# ---------------------------setting initial field-------------------------------
@cfunc(types.field_loop_callback)
def setField_callback(ind, r, E, B, data_double, data_int):
    E[0] = field_amplitude * np.sin(4 * np.pi * r[0] / (xmax - xmin)) #* (abs(r[0]) < l / 4)

sim.field_loop(handler=setField_callback.address)


# ===============================SIMULATION======================================
s = get_pic_steps(64)
for i in range(s):
    sim.advance(time_step=time_step, number_of_iterations=1)
    print(i, "/", s)
