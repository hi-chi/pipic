import numpy
from numba import cfunc

from pipic.tools import CellInterface, Double3, addressof, handler_callback

name = "x_reflector_py"  # extension's name
DataDouble = numpy.zeros((2,), dtype=numpy.double)  # extension's data


@cfunc(handler_callback)
def Handler(CI_I, CI_D, CI_F, CI_P, CI_NP, data_double, data_int):
    """
    function called for each cell and for each requested type (or cells)
    that can be used to modify states, add or remove particles;
    'data_double' refers to 'DataDouble' (usage must be thread-safe)
    """
    C = CellInterface(CI_I, CI_D, CI_F, CI_P, CI_NP)  # unpacking data
    r = Double3(0, 0, 0)  # memory allocation (can be expensive in numba)
    p = Double3(0, 0, 0)
    for ip in range(C.particle_subset_size):
        C.get_r(ip, r)  # setting r to the coordinate of ip-th particle
        if (
            r.x > data_double[0] - 0.5 * data_double[1]
            and r.x < data_double[0] + 0.5 * data_double[1]
        ):
            C.get_p(ip, p)  # setting p to the momentum of ip-th particle
            if p.x > 0:
                p.x *= -1
                C.set_p(ip, p)  # setting the momentum of ip-th particle to p


def handler(location, thickness):
    """
    sets the data of the extension based on parameters and
    returns the address of the Handler
    """
    global DataDouble
    DataDouble[0] = location
    DataDouble[1] = thickness
    return Handler.address


def data_double():
    """
    returns the address of extension's data
    """
    return addressof(DataDouble)
