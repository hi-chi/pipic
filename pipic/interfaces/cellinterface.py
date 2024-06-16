"""
-------------------------------------------------------------------------------------------------------
This file is part of pi-PIC.
pi-PIC, Copyright 2023 Arkady Gonoskov
---------------------------------------------------------------------------------------------------------
pi-PIC is free software: you can redistribute it and/or modify it under the terms of the GNU General
Public License as published by the Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

pi-PIC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with pi-PIC. If not, se
<https://www.gnu.org/licenses/>.
---------------------------------------------------------------------------------------------------------
Website: https://github.com/hi-chi/pipic
Contact: arkady.gonoskov@gu.se.
-------------------------------------------------------------------------------------------------------*/
"""

import numpy as np
from numba import int32, types
from numba.experimental import jitclass

from pipic.types import Double3

# __all__ is used to not pollute the namespace when using 'from xyz import *'
__all__ = ["CellInterface"]


_struct_cellInterface = [
    ("I", types.CPointer(types.int32)),
    ("D", types.CPointer(types.double)),
    ("F", types.CPointer(types.double)),
    ("P", types.CPointer(types.double)),
    ("NP", types.CPointer(types.double)),
    ("c_", types.double[:]),
    ("c_init", int32),
]

_double_uint64_1 = 4.9406564584124654417657e-324  # double with binary code of uint64(1), workaround for particle type setting


@jitclass(_struct_cellInterface)
class CellInterface:
    def __init__(self, I, D, F, P, NP):  # noqa: E741
        self.I = I
        self.D = D
        self.F = F
        self.P = P
        self.NP = NP
        self.c_init = 0

    @property
    def ix(self):
        """x-index of the cell being processed"""
        return self.I[0]

    @property
    def iy(self):
        """y-index of the cell being processed"""
        return self.I[1]

    @property
    def iz(self):
        """z-index of the cell being processed"""
        return self.I[2]

    @property
    def nx(self):
        """x-size of the grid"""
        return self.I[3]

    @property
    def ny(self):
        """y-size of the grid"""
        return self.I[4]

    @property
    def nz(self):
        """z-size of the grid"""
        return self.I[5]

    @property
    def dim(self):
        """dimensionality of the simulation at hand"""
        return self.I[6]

    @property
    def number_of_attributes(
        self,
    ):
        """number of additional attributes of particles"""
        return self.I[7]

    @property
    def particle_type_index(
        self,
    ):
        """type index of the particles being processed, can be fetched via getTypeIndex(typeName)"""
        return self.I[8]

    @property
    def particle_subset_size(
        self,
    ):
        """size of the subset of particles to be processed"""
        return self.I[9]

    @property
    def particle_buffer_capacity(self):
        """capacity of the buffer for new particles, cannot be exceeded"""
        return self.I[10]

    @property
    def particle_buffer_size(self):
        return self.I[11]

    @property
    def grid_type(self):
        return self.I[12]

    @property
    def thread_num(self):
        """number of current thread"""
        return self.I[13]

    @property
    def rng_seed(self):
        """random integer (from -2147483648 to 2147483647, generated for each cell/iteration) to be used as a seed for keeping deterministism"""
        return self.I[14]

    @property
    def global_xmin(self):
        """x-minimum of the entire simulation region"""
        return self.D[0]

    @property
    def global_ymin(self):
        """y-minimum of the entire simulation region"""
        return self.D[1]

    @property
    def global_zmin(self):
        """z-minimum of the entire simulation region"""
        return self.D[2]

    @property
    def global_xmax(self):
        """x-maximum of the entire simulation region"""
        return self.D[3]

    @property
    def global_ymax(self):
        """y-maximum of the entire simulation region"""
        return self.D[4]

    @property
    def global_zmax(self):
        """z-maximum of the entire simulation region"""
        return self.D[5]

    @property
    def step_x(self):
        """x-step of the grid"""
        return self.D[6]

    @property
    def step_y(self):
        """y-step of the grid"""
        return self.D[7]

    @property
    def step_z(self):
        """z-step of the grid"""
        return self.D[8]

    @property
    def inv_step_x(self):
        """inverse step along x"""
        return self.D[9]

    @property
    def inv_step_y(self):
        """inverse step along y"""
        return self.D[10]

    @property
    def inv_step_z(self):
        """inverse step along z"""
        return self.D[11]

    @property
    def time_step(self):
        """time step"""
        return self.D[12]

    @property
    def particle_charge(self):
        """charge of particles being processed"""
        return self.D[13]

    @property
    def particle_mass(self):
        """mass of particles being processed"""
        return self.D[14]

    @property
    def cell_min(self):
        """minimum of the cell being processed"""
        return Double3(
            self.D[0] + self.I[0] * self.D[6],
            self.D[1] + self.I[1] * self.D[7],
            self.D[2] + self.I[2] * self.D[8],
        )

    @property
    def cell_max(self):
        """minimum of the cell being processed"""
        return Double3(
            self.D[0] + (self.I[0] + 1) * self.D[6],
            self.D[1] + (self.I[1] + 1) * self.D[7],
            self.D[2] + (self.I[2] + 1) * self.D[8],
        )

    def get_r(self, idx, r):
        """r is 3d position coordinate."""
        size = 8 + self.I[7]
        r.x = self.P[idx * size + 0]
        r.y = self.P[idx * size + 1]
        r.z = self.P[idx * size + 2]

    def get_p(self, idx, p):
        """p is 3d momentum coordinate."""
        size = 8 + self.I[7]
        p.x = self.P[idx * size + 3]
        p.y = self.P[idx * size + 4]
        p.z = self.P[idx * size + 5]

    def get_w(self, idx, w):
        """w is particle weight."""
        size = 8 + self.I[7]
        w = self.P[idx * size + 6]  # noqa: F841

    def get_id_d(self, idx, id_d):
        """id given as a 'double'"""
        size = 8 + self.I[7]
        id_d = self.P[idx * size + 7]  # noqa: F841

    def get_a(self, idx, attribute_idx, a):
        """a is particle attribute."""
        size = 8 + self.I[7]
        a = self.P[idx * size + 7 + attribute_idx]  # noqa: F841

    def set_p(self, idx, p):
        size = 8 + self.I[7]
        self.P[idx * size + 3] = p.x
        self.P[idx * size + 4] = p.y
        self.P[idx * size + 5] = p.z

    def set_w(self, idx, w):
        size = 8 + self.I[7]
        self.P[idx * size + 6] = w

    def set_a(self, idx, attribute_idx, a):
        size = 8 + self.I[7]
        self.P[idx * size + 7 + attribute_idx] = a

    def add_particle(self, r, p, w, type_index):
        """additional attributes (if any) myst be set by particle_buffer_set_a()"""
        if self.I[11] < self.I[10]:
            self.I[11] += 1
            ip = self.I[11] - 1
            size = 8 + self.I[7]
            self.NP[ip * size + 0] = r.x
            self.NP[ip * size + 1] = r.y
            self.NP[ip * size + 2] = r.z
            self.NP[ip * size + 3] = p.x
            self.NP[ip * size + 4] = p.y
            self.NP[ip * size + 5] = p.z
            self.NP[ip * size + 6] = w
            self.NP[ip * size + 7] = type_index * _double_uint64_1
            return True
        else:
            return False  # indicates that the buffer is overloaded (it will be extended on the next call of handler)

    def particle_buffer_set_a(self, idx, attribute_idx, a):
        size = 8 + self.I[7]
        self.NP[idx * size + 7 + attribute_idx] = a

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
