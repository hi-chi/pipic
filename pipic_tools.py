import pipic
import ctypes
from numba import cfunc, types, carray
import numpy

def addressOf(data):
    if data.dtype == numpy.double:
        return ctypes.addressof(data.ctypes.data_as(ctypes.POINTER(ctypes.c_double)).contents)
    if data.dtype == int:
        return ctypes.addressof(data.ctypes.data_as(ctypes.POINTER(ctypes.c_int)).contents)
        
type_addParticles = types.double(types.CPointer(types.double), types.CPointer(types.double), types.CPointer(types.int32))
type_particleLoop = types.void(types.CPointer(types.double), types.CPointer(types.double), types.CPointer(types.double), types.CPointer(types.uint64), types.CPointer(types.double), types.CPointer(types.int32))
type_fieldLoop = types.void(types.CPointer(types.int32), types.CPointer(types.double), types.CPointer(types.double), types.CPointer(types.double), types.CPointer(types.double), types.CPointer(types.int32))
type_it2r = types.void(types.CPointer(types.int32), types.CPointer(types.double), types.CPointer(types.double), types.CPointer(types.int32))
type_field2data = types.void(types.CPointer(types.int32), types.CPointer(types.double), types.CPointer(types.double), types.CPointer(types.double), types.CPointer(types.double), types.CPointer(types.int32))

electronCharge = pipic.electronCharge
electronMass = pipic.electronMass
lightVelocity = pipic.lightVelocity
