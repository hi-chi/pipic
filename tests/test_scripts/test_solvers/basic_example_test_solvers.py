import pipic
from pipic import consts, types
import numpy as np
from numba import cfunc, carray
import argparse


parser = argparse.ArgumentParser(description="Run basic pi-PIC example with a selected solver.")
parser.add_argument(
    "--solver",
    default="electrostatic_1d",
    help="Solver to use.",
)
parser.add_argument(
    "--steps",
    type=int,
    default=64,
    help="Number of PIC advancements.",
)
args = parser.parse_args()
solver_name = args.solver
s = max(1, args.steps)


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
sim = pipic.init(solver=solver_name, nx=nx, xmin=xmin, xmax=xmax)

# ------------------------------adding electrons---------------------------------
@cfunc(types.add_particles_callback)
def density_callback(r, data_double, data_int):
    return density #* (abs(r[0]) < l / 4)

print(consts.electron_charge)

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


# -------------------------preparing output of Ex(x)-----------------------------
Ex = np.zeros((nx,1), dtype=np.double)


@cfunc(types.field_loop_callback)
def get_Ex(ind, r, E, B, data_double, data_int):
    data = carray(data_double, Ex.shape, dtype=np.double)
    if ind[0] < nx:
        data[ind[0],0] = E[0]



# ===============================SIMULATION======================================
for i in range(s):

    sim.advance(time_step=time_step, number_of_iterations=1)    
    sim.field_loop(
        handler=get_Ex.address,
        data_double=pipic.addressof(Ex),
    )
    print(Ex.sum())
