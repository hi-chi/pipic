import numba as nb
from numba.experimental import jitclass

# __all__ is used to not pollute the namespace when using 'from xyz import *'
__all__ = ['Int3']


@jitclass([('x', nb.types.int32), ('y', nb.types.int32), ('z', nb.types.int32),])
class Int3:
    def __init__(self, x = 0, y = 0, z = 0):
        self.x = x
        self.y = y
        self.z = z
