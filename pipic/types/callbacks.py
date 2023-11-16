import numba as nb

# __all__ is used to not pollute the namespace when using 'from xyz import *'
__all__ = ['add_particles_callback',
           'particle_loop_callback',
           'field_loop_callback',
           'it2r_callback',
           'field2data_callback',
           'handler_callback']


add_particles_callback = nb.types.double(nb.types.CPointer(nb.types.double),
                                         nb.types.CPointer(nb.types.double),
                                         nb.types.CPointer(nb.types.int32))

particle_loop_callback = nb.types.void(nb.types.CPointer(nb.types.double),
                                       nb.types.CPointer(nb.types.double),
                                       nb.types.CPointer(nb.types.double),
                                       nb.types.CPointer(nb.types.uint64),
                                       nb.types.CPointer(nb.types.double),
                                       nb.types.CPointer(nb.types.int32))

field_loop_callback = nb.types.void(nb.types.CPointer(nb.types.int32),
                                    nb.types.CPointer(nb.types.double),
                                    nb.types.CPointer(nb.types.double),
                                    nb.types.CPointer(nb.types.double),
                                    nb.types.CPointer(nb.types.double),
                                    nb.types.CPointer(nb.types.int32))

it2r_callback = nb.types.void(nb.types.CPointer(nb.types.int32),
                              nb.types.CPointer(nb.types.double),
                              nb.types.CPointer(nb.types.double),
                              nb.types.CPointer(nb.types.int32))

field2data_callback = nb.types.void(nb.types.CPointer(nb.types.int32),
                                    nb.types.CPointer(nb.types.double),
                                    nb.types.CPointer(nb.types.double),
                                    nb.types.CPointer(nb.types.double),
                                    nb.types.CPointer(nb.types.double),
                                    nb.types.CPointer(nb.types.int32))

handler_callback = nb.types.void(nb.types.CPointer(nb.types.int32),
                                 nb.types.CPointer(nb.types.double),
                                 nb.types.CPointer(nb.types.double),
                                 nb.types.CPointer(nb.types.double),
                                 nb.types.CPointer(nb.types.double),
                                 nb.types.CPointer(nb.types.double),
                                 nb.types.CPointer(nb.types.int32))