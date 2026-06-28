# As a test of extension 'qed_volokitin2023' this code reproduces data of fig. 10 in [A. Gonoskov et al. PRE (2015), arXiv:1412.6426]

import pipic
from pipic import consts, types
import numpy as np
from numba import cfunc
from pipic.extensions import qed_volokitin2023
import math
import sys


def get_pic_steps(default_steps):
    if "--steps" in sys.argv:
        idx = sys.argv.index("--steps")
        return max(1, int(sys.argv[idx + 1]))
    return default_steps

# ===========================SIMULATION INITIALIZATION===========================
gamma = 2e5
E_cr = (
    consts.electron_mass**2 * consts.light_velocity**3 / (abs(consts.electron_charge) * consts.hbar)
)
B_y = 0.2 * E_cr

t_rad = (
    3.85
    * gamma ** (1 / 3.0)
    * (E_cr / B_y) ** (2 / 3)
    * consts.hbar**2
    / (consts.electron_mass * consts.light_velocity * consts.electron_charge**2)
)
t_sim = 8 * t_rad
N_steps = 32
time_step = t_sim / N_steps

density = 1e18
nx = 64
xmin, xmax = -nx * consts.light_velocity * time_step, nx * consts.light_velocity * time_step

# ---------------------setting solver and simulation region----------------------
sim = pipic.init(solver="fourier_boris", nx=nx, xmin=xmin, xmax=xmax)


# ------------------------------adding electrons---------------------------------
@cfunc(types.add_particles)
def density(r, data_double, data_int):
    return density


sim.add_particles(
    name="electron",
    number=1 * nx,
    charge=consts.electron_charge,
    mass=consts.electron_mass,
    temperature=0,
    density=density.address,
)


# ----------------------adding photons and postirons----------------------------
@cfunc(types.add_particles)
def null_callback(r, data_double, data_int):
    return 0


sim.add_particles(
    name="photon", number=0, charge=0, mass=0, temperature=0, density=null_callback.address
)

sim.add_particles(
    name="positron",
    number=0,
    charge=-consts.electron_charge,
    mass=consts.electron_mass,
    temperature=0,
    density=null_callback.address,
)

sim.set_rng_seed(14)
# ------------------------------adding extension---------------------------------
qed_handler = qed_volokitin2023.handler(
    electron_type=sim.get_type_index("electron"),
    positron_type=sim.get_type_index("positron"),
    photon_type=sim.get_type_index("photon"),
)
sim.add_handler(
    name=qed_volokitin2023.name, subject="electron, positron, photon", handler=qed_handler
)


# --------------------------setting electron momentum----------------------------
@cfunc(types.particle_loop)
def set_p(r, p, w, id, data_double, data_int):
    p[0] = gamma * consts.electron_mass * consts.light_velocity
    p[1] = 0
    p[2] = 0


sim.particle_loop(name="electron", handler=set_p.address)


# ---------------------------setting initial field-------------------------------
@cfunc(types.field_loop)
def setField(ind, r, E, B, data_double, data_int):
    E[0] = 0
    E[1] = 0
    E[2] = 0
    B[0] = 0
    B[1] = B_y
    B[2] = 0


sim.field_loop(handler=setField.address)

# -------------preparing output for electron distribution f(x, px)--------------
N_e_ep = np.zeros((1,), dtype=np.double)


@cfunc(types.particle_loop)
def N_e_ep_cb(r, p, w, id, data_double, data_int):
    if (
        math.sqrt(p[0] ** 2 + p[1] ** 2 + p[2] ** 2)
        > (1e-3) * gamma * consts.electron_mass * consts.light_velocity
    ):
        data_double[0] += w[0]
    else:
        w[0] = 0


def get_N_e_ep():
    N_e_ep[0] = 0
    sim.particle_loop(
        name="electron", handler=N_e_ep_cb.address, data_double=pipic.addressof(N_e_ep)
    )
    sim.particle_loop(
        name="positron", handler=N_e_ep_cb.address, data_double=pipic.addressof(N_e_ep)
    )
    return N_e_ep[0]


@cfunc(types.particle_loop)
def ph_cb(r, p, w, id, data_double, data_int):
    if (
        math.sqrt(p[0] ** 2 + p[1] ** 2 + p[2] ** 2)
        < (1e-3) * gamma * consts.electron_mass * consts.light_velocity
    ):
        w[0] = 0


def remove_low_en_ph():
    sim.particle_loop(name="photon", handler=ph_cb.address)


# ===============================SIMULATION======================================
N0 = get_N_e_ep()
Nsteps = get_pic_steps(round(t_sim / time_step))
for i in range(Nsteps):
    sim.advance(time_step=time_step, number_of_iterations=1)
    remove_low_en_ph()
    if i % (N_steps / 32) == 0:
        print(
            i / (N_steps / 32),
            "/",
            Nsteps,
            ": N_{ee+} =",
            get_N_e_ep() / N0,
            ", N_{macro-particles} =",
            sim.get_number_of_particles(),
        )
