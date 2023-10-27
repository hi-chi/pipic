import pipic
from pipic_tools import *
name = 'x_reflector_py'
DataDouble = numpy.zeros((2, ), dtype=numpy.double) 
@cfunc(type_handler)
def Handler(CI_I, CI_D, CI_F, CI_P, CI_NP, dataDouble, dataInt):
    C = cellInterface(CI_I, CI_D, CI_F, CI_P, CI_NP) # unpacking data   
    r = double3(0, 0, 0) # memory allocation (can be expensive in numba)
    p = double3(0, 0, 0) 
    for ip in range(C.PSize):
        C.get_r(ip, r)
        if  r.x > dataDouble[0] - 0.5*dataDouble[1]\
        and r.x < dataDouble[0] + 0.5*dataDouble[1]:
            C.get_p(ip, p)
            if p.x > 0:
                p.x *= -1
                C.set_p(ip, p)
def handler(location, thickness):
    global DataDouble
    DataDouble[0] = location
    DataDouble[1] = thickness
    return Handler.address
def dataDouble():
    return addressOf(DataDouble) 