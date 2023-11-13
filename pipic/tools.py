'''
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
// Description: Here we define necessary cfunc decorators, as well as interfaces for code extension.
'''
import _pipic as pipic
import ctypes
from numba import cfunc, types, carray, double, int32
import numpy, os, time
from numba.experimental import jitclass
from math import *

def addressOf(data):
    if data.dtype == numpy.double:
        return ctypes.addressof(data.ctypes.data_as(ctypes.POINTER(ctypes.c_double)).contents)
    if data.dtype == numpy.intc:
        return ctypes.addressof(data.ctypes.data_as(ctypes.POINTER(ctypes.c_int)).contents)
        
type_addParticles = types.double(types.CPointer(types.double), types.CPointer(types.double), types.CPointer(types.int32))
type_particleLoop = types.void(types.CPointer(types.double), types.CPointer(types.double), types.CPointer(types.double), types.CPointer(types.uint64), types.CPointer(types.double), types.CPointer(types.int32))
type_fieldLoop = types.void(types.CPointer(types.int32), types.CPointer(types.double), types.CPointer(types.double), types.CPointer(types.double), types.CPointer(types.double), types.CPointer(types.int32))
type_it2r = types.void(types.CPointer(types.int32), types.CPointer(types.double), types.CPointer(types.double), types.CPointer(types.int32))
type_field2data = types.void(types.CPointer(types.int32), types.CPointer(types.double), types.CPointer(types.double), types.CPointer(types.double), types.CPointer(types.double), types.CPointer(types.int32))
type_handler = types.void(types.CPointer(types.int32), types.CPointer(types.double), types.CPointer(types.double), types.CPointer(types.double), types.CPointer(types.double), types.CPointer(types.double), types.CPointer(types.int32))

electronCharge = pipic.electronCharge
electronMass = pipic.electronMass
lightVelocity = pipic.lightVelocity

@jitclass([('x', int32), ('y', int32), ('z', int32),])
class int3:
    def __init__(self, x = 0, y = 0, z = 0):
        self.x = x
        self.y = y
        self.z = z

@jitclass([('x', double), ('y', double), ('z', double),])
class double3:
    def __init__(self, x = 0, y = 0, z = 0):
        self.x = x
        self.y = y
        self.z = z
    def normalize(self):
        r = sqrt(self.x*self.x + self.y*self.y + self.z*self.z)
        if r > 0:
            inv_r = 1/r
            self.x *= inv_r
            self.y *= inv_r
            self.z *= inv_r
    def norm(self) -> float:
        return sqrt(self.x*self.x + self.y*self.y + self.z*self.z)
    def norm2(self) -> float:
        return self.x*self.x + self.y*self.y + self.z*self.z
    def __add__(self, other):
        return double3(self.x + other.x, self.y + other.y, self.z + other.z)
    def __sub__(self, other):
        return double3(self.x - other.x, self.y - other.y, self.z - other.z)
    def __iadd__(self, other):
        self.x += other.x
        self.y += other.y
        self.z += other.z
        return self
    def __mul__(self, other):
        return double3(self.x * other, self.y * other, self.z * other)
    def __imul__(self, other):
        self.x *= other
        self.y *= other
        self.z *= other
        return self
    
    def cross(self, b):
        return double3(self.y*b.z-self.z*b.y, self.z*b.x-self.x*b.z, self.x*b.y-self.y*b.x)
    
    def dot(self, b) -> float:
        return self.x*b.x + self.y*b.y + self.z*b.z


struct_cellInterface = [
    ('I', types.CPointer(types.int32)),
    ('D', types.CPointer(types.double)),
    ('F', types.CPointer(types.double)),
    ('P', types.CPointer(types.double)),
    ('NP', types.CPointer(types.double)),
    ('c_', types.double[:]),
    ('c_init', int32),
]

double_uint64_1 = 4.9406564584124654417657e-324 # double with binary code of uint64(1), workaround for particle type setting

@jitclass(struct_cellInterface)
class cellInterface:
    def __init__(self, I, D, F, P, NP):
        self.I = I
        self.D = D
        self.F = F
        self.P = P
        self.NP = NP
        self.c_init = 0
    @property
    def ix(self): # x-index of the cell being processed
        return self.I[0]
    @property
    def iy(self): # y-index of the cell being processed
        return self.I[1]
    @property
    def iz(self): # z-index of the cell being processed
        return self.I[2]
    @property
    def nx(self): # x-size of the grid
        return self.I[3]
    @property
    def ny(self): # y-size of the grid
        return self.I[4]
    @property
    def nz(self): # z-size of the grid
        return self.I[5]
    @property
    def dim(self): # dimensionality of the simulation at hand
        return self.I[6]
    @property
    def aNum(self): # number of additional attributes of particles
        return self.I[7]
    @property
    def PType(self): # type index of the particles being processed, can be fetched via getTypeIndex(typeName)
        return self.I[8]
    @property
    def PSize(self): # size of the subset of particles to be processed
        return self.I[9]
    @property
    def NPcapacity(self):
        return self.I[10] # capacity of the buffer for new particles, cannot be exceeded
    @property
    def NPSize(self):
        return self.I[11] # number of current thread
    @property
    def gridType(self):
        return self.I[12] # number of current thread
    @property
    def threadNum(self):
        return self.I[13] # number of current thread
    @property
    def gloablMinX(self): # x-minimum of the entire simulation region
        return self.D[0]
    @property
    def gloablMinY(self): # y-minimum of the entire simulation region
        return self.D[1]
    @property
    def gloablMinZ(self): # z-minimum of the entire simulation region
        return self.D[2]
    @property
    def gloablMaxX(self): # x-maximum of the entire simulation region
        return self.D[3]
    @property
    def gloablMaxY(self): # y-maximum of the entire simulation region
        return self.D[4]
    @property
    def gloablMaxZ(self): # z-maximum of the entire simulation region
        return self.D[5]
    @property
    def stepX(self): # x-step of the grid
        return self.D[6]
    @property
    def stepY(self): # y-step of the grid
        return self.D[7]
    @property
    def stepZ(self): # z-step of the grid
        return self.D[8]
    @property
    def invStepX(self): # inverse step along x
        return self.D[9]
    @property
    def invStepY(self): # inverse step along y
        return self.D[10]
    @property
    def invStepZ(self): # inverse step along z
        return self.D[11]
    @property
    def timeStep(self): # time step
        return self.D[12] 
    @property
    def PCharge(self): # charge of particles being processed
        return self.D[13] 
    @property
    def PMass(self): # charge of particles being processed
        return self.D[14] 
    @property
    def cellMin(self): # minimum of the cell being processed
        return double3(self.D[0] + self.I[0]*self.D[6], self.D[1] + self.I[1]*self.D[7], self.D[2] + self.I[2]*self.D[8])
    @property
    def cellMax(self): # minimum of the cell being processed
        return double3(self.D[0] + (self.I[0]+1)*self.D[6], self.D[1] + (self.I[1]+1)*self.D[7], self.D[2] + (self.I[2]+1)*self.D[8])
    
    def get_r(self, ip, r): 
        size = 8 + self.I[7]
        r.x = self.P[ip*size + 0]
        r.y = self.P[ip*size + 1]
        r.z = self.P[ip*size + 2]
    
    def get_p(self, ip, p):
        size = 8 + self.I[7]
        p.x = self.P[ip*size + 3]
        p.y = self.P[ip*size + 4]
        p.z = self.P[ip*size + 5]

    def get_w(self, ip, w):
        size = 8 + self.I[7]
        w = self.P[ip*size + 6]

    def get_id_d(self, ip, id_d):
        size = 8 + self.I[7]
        id_d = self.P[ip*size + 7]
    
    def get_a(self, ip, ia, a):
        size = 8 + self.I[7]
        a = self.P[ip*size + 7 + ia]
        
    def set_p(self, ip, p):
        size = 8 + self.I[7]
        self.P[ip*size + 3] = p.x
        self.P[ip*size + 4] = p.y
        self.P[ip*size + 5] = p.z
    
    def set_w(self, ip, w):
        size = 8 + self.I[7]
        self.P[ip*size + 6] = w

    def set_a(self, ip, ia, a):
        size = 8 + self.I[7]
        self.P[ip*size + 7 + ia] = a
    
    def addParticle(self, r, p, w, typeIndex): # additional attributes (if any) myst be set by NP_set_a()
        if self.I[11] < self.I[10]:
            self.I[11] += 1
            ip = self.I[11]
            size = 8 + self.I[7]
            self.NP[ip*size + 0] = r.x
            self.NP[ip*size + 1] = r.y
            self.NP[ip*size + 2] = r.z
            self.NP[ip*size + 3] = p.x
            self.NP[ip*size + 4] = p.y
            self.NP[ip*size + 5] = p.z
            self.NP[ip*size + 6] = w   
            self.NP[ip*size + 7] = typeIndex*double_uint64_1
            return True
        else:
            return False # indicates that the buffer is overloaded (it will be extended on the next call of handler)
    
    def NP_set_a(self, ip, ia, a):
        size = 8 + self.I[7]
        self.NP[ip*size + 7 + ia] = a

    def interpolateField(self, r, E, B):
        if self.I[12] == 0:
            if self.c_init == 0:
                self.c_ = numpy.zeros(8, dtype=numpy.double)
                self.c_init = 1
            
            self.c_[1] = (r.x - (self.D[0] + self.I[0]*self.D[6]))*self.D[9]
            self.c_[0] = 1 - self.c_[1]
            dim = self.I[6]
            cn = 2
            if dim > 1:
                wy = (r.y - (self.D[1] + self.I[1]*self.D[7]))*self.D[10]
                wy_ = 1 - wy
                for i in range(2):
                    self.c_[i+2] = self.c_[i]*wy
                    self.c_[i] *= wy_
                cn = 4
            if dim > 2:
                wz = (r.z - (self.D[2] + self.I[2]*self.D[8]))*self.D[11]
                wz_ = 1 - wz
                for i in range(4):
                    self.c_[i+4] = self.c_[i]*wz
                    self.c_[i] *= wz_
                cn = 8
            E.x, E.y, E.z, B.x, B.y, B.z = 0, 0, 0, 0, 0, 0
            for i in range(cn):
                E.x += self.c_[i]*self.F[6*i]
                E.y += self.c_[i]*self.F[6*i+1]
                E.z += self.c_[i]*self.F[6*i+2]
                B.x += self.c_[i]*self.F[6*i+3]
                B.y += self.c_[i]*self.F[6*i+4]
                B.z += self.c_[i]*self.F[6*i+5]  
            return True
        else:
            return False # grid type is unknown
