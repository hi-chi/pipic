# Making extensions

$\pi$-PIC offers a possibility to develop extensions in Python, C/C++, Fortran and other programming languages. Apart from modifying the field state via `fieldLoop()`, the extensions can modify, add and remove particles depending on the local field state. For example, this can be used to account for ionization, radiation reaction, QED particle and photon generation, etc.

Python extensions
--
1.	Copy the following $\pi$-PIC files to a local folder:
    - so-file of $\pi$-PIC
    - `pipic_tools.py`
2.	Create an `<extension>.py` file and develop the code of extension following the layout of examples in `pipic/extensions` (for more controls see methods of struct cellInterface in `pipic_tools.py`).
3.	To use the extension, place `<extension>.py`, so-file of $\pi$-PIC and `pipic_tools.py` in a folder with Python file of your simulation. Add and configure your extension (replace `<*>` with respective entries):
    ```
    - import <extension> 
    - <simulation>. addCellHandler(name=<extension>.name, subject=’<particleType>’, \
    handler=<extension>.handler(<handler parameters>)
    ```

As an example see below an extension that can be used to reflect particles from an $x$-limited region.
```
import pipic
from pipic_tools import *
name = 'x_reflector_py'
DataDouble = numpy.zeros((2, ), dtype=numpy.double) 
@cfunc(type_handler)
def Handler(CI_I, CI_D, CI_F, CI_P, CI_NP, dataDouble, dataInt):
    C = cellInterface(CI_I, CI_D, CI_F, CI_P, CI_NP) # unpacking data   
    r = double3(0, 0, 0) # memory allocation (can be expensive in numba)
    p = double3(0, 0, 0) 
    for ip in range(C.PSize): # making a loop over particles in the cell being processed
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
```

C/C++ extensions
--
1.	From `src` folder copy the following interface files in a local folder:
    - `primitives.h`
    - `interfaces.h`
    - `CMakeLists.txt`
2.	In the file `CMakeLists.txt` replace “pipic” with the extension name.
3.	Download `pybind11` to the local folder
4.	Create an `<extension>.cpp` file and develop the code of extension following the layout of examples in `pipic/extensions` (for more controls see methods of struct cellInterface in `interfaces.h`).
5.	Compile your extension by sequentially running:
    ```
    - cmake .
    - make
    ```
6.	To use the extension, place the generated so-file, so-file of pi-PIC and pipic_tools.py in a folder with Python file of your simulation. Add and configure your extension (replace `<*>` with respective entries):
    ```
    - import <extension> 
    - <simulation>. addCellHandler(name=<extension>.name, subject=’<particleType>’, \
    handler=<extension>.handler(<handler parameters>)
    ```

As an example see below an extension that can be used to convert type of particles in an $x$-limited region.

```
#include "interfaces.h"
#include <pybind11/pybind11.h>
#include "pybind11/stl.h"
#include <pybind11/operators.h>

const string name = "x_converter_c";
static double xMin, xMax; // limits of the action
static int typeTo; // the typeIndex to set for affected particles

void Handler(int *I, double *D, double *F, double *P, double *NP, double *dataDouble, int *dataInt){
    cellInterface CI(I, D, F, P, NP);
    // changing type requires replacing existing particle with a new particle with altered typeIndex
    for(int ip = 0; ip < CI.PSize; ip++)
        if((CI.P(ip)->r.x > xMin)&&(CI.P(ip)->r.x < xMax))
            if(CI.NPSize < CI.NPcapacity){ // checking if the buffer permits adding a particle 
                *CI.NP(CI.NPSize) = *CI.P(ip); // copy particle to a new particle (buffer)
                CI.NP(CI.NPSize)->id = typeTo; // put the new particle type in id (convention)
                CI.P(ip)->w = 0; // zero weight indicates that the paricle has to be removed
                CI.NPSize++;
            }
    
};

int64_t handler(double location, double thickness, int TypeTo){ // initialization
    xMin = location - 0.5*thickness;
    xMax = location + 0.5*thickness;
    typeTo = TypeTo;
    return (int64_t)Handler;
};

namespace py = pybind11;
PYBIND11_MODULE(x_converter_c, object) {
    object.attr("name") = name;
    object.def("handler", &handler, py::arg("location"), py::arg("thickness"), py::arg("typeTo"));
}
```