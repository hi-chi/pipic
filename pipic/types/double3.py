import numba as nb
from numba.experimental import jitclass

# __all__ is used to not pollute the namespace when using 'from xyz import *'
__all__ = ['Double3']


@jitclass([('x', nb.types.double), ('y', nb.types.double), ('z', nb.types.double),])
class Double3:
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
        return Double3(self.x + other.x, self.y + other.y, self.z + other.z)
    def __sub__(self, other):
        return Double3(self.x - other.x, self.y - other.y, self.z - other.z)
    def __iadd__(self, other):
        self.x += other.x
        self.y += other.y
        self.z += other.z
        return self
    def __mul__(self, other):
        return Double3(self.x * other, self.y * other, self.z * other)
    def __imul__(self, other):
        self.x *= other
        self.y *= other
        self.z *= other
        return self

    def cross(self, b):
        return Double3(self.y * b.z - self.z * b.y, self.z * b.x - self.x * b.z, self.x * b.y - self.y * b.x)

    def dot(self, b) -> float:
        return self.x*b.x + self.y*b.y + self.z*b.z