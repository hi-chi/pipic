import numpy as np
from numba import types, int32
from numba.experimental import jitclass

# __all__ is used to not pollute the namespace when using 'from xyz import *'
__all__ = ['CellInterface']

_struct_cellInterface = [
    ('I', types.CPointer(types.int32)),
    ('D', types.CPointer(types.double)),
    ('F', types.CPointer(types.double)),
    ('P', types.CPointer(types.double)),
    ('NP', types.CPointer(types.double)),
    ('c_', types.double[:]),
    ('c_init', int32),
]

_double_uint64_1 = 4.9406564584124654417657e-324  # double with binary code of uint64(1), workaround for particle type setting


@jitclass(_struct_cellInterface)
class CellInterface:
    def __init__(self, I, D, F, P, NP):
        self.I = I
        self.D = D
        self.F = F
        self.P = P
        self.NP = NP
        self.c_init = 0

    @property
    def ix(self):  # x-index of the cell being processed
        return self.I[0]

    @property
    def iy(self):  # y-index of the cell being processed
        return self.I[1]

    @property
    def iz(self):  # z-index of the cell being processed
        return self.I[2]

    @property
    def nx(self):  # x-size of the grid
        return self.I[3]

    @property
    def ny(self):  # y-size of the grid
        return self.I[4]

    @property
    def nz(self):  # z-size of the grid
        return self.I[5]

    @property
    def dim(self):  # dimensionality of the simulation at hand
        return self.I[6]

    @property
    def a_num(self):  # number of additional attributes of particles
        return self.I[7]

    @property
    def p_type(self):  # type index of the particles being processed, can be fetched via getTypeIndex(typeName)
        return self.I[8]

    @property
    def PSize(self):  # size of the subset of particles to be processed
        return self.I[9]

    @property
    def NPcapacity(self):
        return self.I[10]  # capacity of the buffer for new particles, cannot be exceeded

    @property
    def NPSize(self):
        return self.I[11]  # number of current thread

    @property
    def grid_type(self):
        return self.I[12]  # number of current thread

    @property
    def thread_num(self):
        return self.I[13]  # number of current thread

    @property
    def global_min_x(self):  # x-minimum of the entire simulation region
        return self.D[0]

    @property
    def global_min_y(self):  # y-minimum of the entire simulation region
        return self.D[1]

    @property
    def global_min_z(self):  # z-minimum of the entire simulation region
        return self.D[2]

    @property
    def global_max_x(self):  # x-maximum of the entire simulation region
        return self.D[3]

    @property
    def global_max_y(self):  # y-maximum of the entire simulation region
        return self.D[4]

    @property
    def global_max_z(self):  # z-maximum of the entire simulation region
        return self.D[5]

    @property
    def step_x(self):  # x-step of the grid
        return self.D[6]

    @property
    def step_y(self):  # y-step of the grid
        return self.D[7]

    @property
    def step_z(self):  # z-step of the grid
        return self.D[8]

    @property
    def inv_step_x(self):  # inverse step along x
        return self.D[9]

    @property
    def inv_step_y(self):  # inverse step along y
        return self.D[10]

    @property
    def inv_step_z(self):  # inverse step along z
        return self.D[11]

    @property
    def time_step(self):  # time step
        return self.D[12]

    @property
    def p_charge(self):  # charge of particles being processed
        return self.D[13]

    @property
    def p_mass(self):  # charge of particles being processed
        return self.D[14]

    @property
    def cell_min(self):  # minimum of the cell being processed
        return Double3(self.D[0] + self.I[0] * self.D[6], self.D[1] + self.I[1] * self.D[7],
                       self.D[2] + self.I[2] * self.D[8])

    @property
    def cell_max(self):  # minimum of the cell being processed
        return Double3(self.D[0] + (self.I[0] + 1) * self.D[6], self.D[1] + (self.I[1] + 1) * self.D[7],
                       self.D[2] + (self.I[2] + 1) * self.D[8])

    def get_r(self, ip, r):
        size = 8 + self.I[7]
        r.x = self.P[ip * size + 0]
        r.y = self.P[ip * size + 1]
        r.z = self.P[ip * size + 2]

    def get_p(self, ip, p):
        size = 8 + self.I[7]
        p.x = self.P[ip * size + 3]
        p.y = self.P[ip * size + 4]
        p.z = self.P[ip * size + 5]

    def get_w(self, ip, w):
        size = 8 + self.I[7]
        w = self.P[ip * size + 6]

    def get_id_d(self, ip, id_d):
        size = 8 + self.I[7]
        id_d = self.P[ip * size + 7]

    def get_a(self, ip, ia, a):
        size = 8 + self.I[7]
        a = self.P[ip * size + 7 + ia]

    def set_p(self, ip, p):
        size = 8 + self.I[7]
        self.P[ip * size + 3] = p.x
        self.P[ip * size + 4] = p.y
        self.P[ip * size + 5] = p.z

    def set_w(self, ip, w):
        size = 8 + self.I[7]
        self.P[ip * size + 6] = w

    def set_a(self, ip, ia, a):
        size = 8 + self.I[7]
        self.P[ip * size + 7 + ia] = a

    def add_particle(self, r, p, w, typeIndex):  # additional attributes (if any) myst be set by NP_set_a()
        if self.I[11] < self.I[10]:
            self.I[11] += 1
            ip = self.I[11]
            size = 8 + self.I[7]
            self.NP[ip * size + 0] = r.x
            self.NP[ip * size + 1] = r.y
            self.NP[ip * size + 2] = r.z
            self.NP[ip * size + 3] = p.x
            self.NP[ip * size + 4] = p.y
            self.NP[ip * size + 5] = p.z
            self.NP[ip * size + 6] = w
            self.NP[ip * size + 7] = typeIndex * _double_uint64_1
            return True
        else:
            return False  # indicates that the buffer is overloaded (it will be extended on the next call of handler)

    def NP_set_a(self, ip, ia, a):
        size = 8 + self.I[7]
        self.NP[ip * size + 7 + ia] = a

    def interpolate_field(self, r, E, B):
        if self.I[12] == 0:
            if self.c_init == 0:
                self.c_ = np.zeros(8, dtype=np.double)
                self.c_init = 1

            self.c_[1] = (r.x - (self.D[0] + self.I[0] * self.D[6])) * self.D[9]
            self.c_[0] = 1 - self.c_[1]
            dim = self.I[6]
            cn = 2
            if dim > 1:
                wy = (r.y - (self.D[1] + self.I[1] * self.D[7])) * self.D[10]
                wy_ = 1 - wy
                for i in range(2):
                    self.c_[i + 2] = self.c_[i] * wy
                    self.c_[i] *= wy_
                cn = 4
            if dim > 2:
                wz = (r.z - (self.D[2] + self.I[2] * self.D[8])) * self.D[11]
                wz_ = 1 - wz
                for i in range(4):
                    self.c_[i + 4] = self.c_[i] * wz
                    self.c_[i] *= wz_
                cn = 8
            E.x, E.y, E.z, B.x, B.y, B.z = 0, 0, 0, 0, 0, 0
            for i in range(cn):
                E.x += self.c_[i] * self.F[6 * i]
                E.y += self.c_[i] * self.F[6 * i + 1]
                E.z += self.c_[i] * self.F[6 * i + 2]
                B.x += self.c_[i] * self.F[6 * i + 3]
                B.y += self.c_[i] * self.F[6 * i + 4]
                B.z += self.c_[i] * self.F[6 * i + 5]
            return True
        else:
            return False  # grid type is unknown
