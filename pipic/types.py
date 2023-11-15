import numba as nb

# __all__ is used to not pollute the namespace when using 'from xyz import *'
__all__ = ['add_particles',
           'particle_loop',
           'field_loop',
           'it2r',
           'field2data',
           'handler']

add_particles = nb.types.double(nb.types.CPointer(nb.types.double),
                                nb.types.CPointer(nb.types.double),
                                nb.types.CPointer(nb.types.int32))

particle_loop = nb.types.void(nb.types.CPointer(nb.types.double),
                              nb.types.CPointer(nb.types.double),
                              nb.types.CPointer(nb.types.double),
                              nb.types.CPointer(nb.types.uint64),
                              nb.types.CPointer(nb.types.double),
                              nb.types.CPointer(nb.types.int32))

field_loop = nb.types.void(nb.types.CPointer(nb.types.int32),
                           nb.types.CPointer(nb.types.double),
                           nb.types.CPointer(nb.types.double),
                           nb.types.CPointer(nb.types.double),
                           nb.types.CPointer(nb.types.double),
                           nb.types.CPointer(nb.types.int32))

it2r = nb.types.void(nb.types.CPointer(nb.types.int32),
                     nb.types.CPointer(nb.types.double),
                     nb.types.CPointer(nb.types.double),
                     nb.types.CPointer(nb.types.int32))

field2data = nb.types.void(nb.types.CPointer(nb.types.int32),
                           nb.types.CPointer(nb.types.double),
                           nb.types.CPointer(nb.types.double),
                           nb.types.CPointer(nb.types.double),
                           nb.types.CPointer(nb.types.double),
                           nb.types.CPointer(nb.types.int32))

handler = nb.types.void(nb.types.CPointer(nb.types.int32),
                        nb.types.CPointer(nb.types.double),
                        nb.types.CPointer(nb.types.double),
                        nb.types.CPointer(nb.types.double),
                        nb.types.CPointer(nb.types.double),
                        nb.types.CPointer(nb.types.double),
                        nb.types.CPointer(nb.types.int32))
