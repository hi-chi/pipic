# Making extensions

$\pi$-PIC offers a possibility to develop extensions in Python, C/C++, Fortran and other programming languages that can be used to produce a callable function for Python. Extensions can be developed, compiled and used independently, as well as contributed to the toolbox delivered with $\pi$-PIC (see [Python extensions](https://github.com/hi-chi/pipic/blob/main/pipic/extensions) and [C++ extensions](https://github.com/hi-chi/pipic/blob/main/src/extensions)). Apart from modifying the field state via `fieldLoop()`, the extensions can modify, add and remove particles based on the local field state. For example, this can be used to account for ionization, radiation reaction, QED particle and photon generation, etc.

Python extensions
--
Creating an extension within Python is based on developing a handler in the form of a [Numba C callback](https://numba.readthedocs.io/en/stable/user/cfunc.html) that is called for each cell and offers a possibility to add new particles and/or modify states of all particles in the cell sequentially for all requested types. When called the handler receives data pointers to be passed to a `CellInterface` structure that offers basic information and options for modifying existing and/or adding new particles within the processed cell (see complete list of interfaces in [`cellinterface.py`](https://github.com/hi-chi/pipic/blob/main/pipic/interfaces/cellinterface.py)). The state of the field can be modified whenever needed using `field_loop()`, see [User interfaces](https://github.com/hi-chi/pipic/blob/main/docs/guides/INTERFACES.md). Despite some limitations on what can be used within Numba C callbacks they can provide sufficiently high performance and flexibility for many cases of interest.

The procedure for developing an extension within Python can look as follows:
1.	Create an `<extension>.py` file and develop the code of extension following the layout of an [example](https://github.com/hi-chi/pipic/blob/main/pipic/extensions/x_reflector_py.py)
2.	To use the extension, place `<extension>.py` in a folder with Python file of your simulation.
3.	In the file of your simulation add and configure your extension (replace `<*>` with respective entries):
    ```
    - import <extension> 
    - <simulation>. add_cell_handler(name=<extension>.name, subject=’<particleType>’, \
    handler=<extension>.handler(<handler parameters>)
    ```
    `<particleType>` is a list of particle types to be processed by the extension; `all_types` can be used to process all types, `cells` to call the handler for each cell independently of particles there (can be used to add new particles).

As an example see an extension [`x_reflector_py.py`](https://github.com/hi-chi/pipic/blob/main/pipic/extensions/x_reflector_py.py) that can be used to reflect particles from an $x$-limited region and its use in [`x_reflector_py_test.py`](https://github.com/hi-chi/pipic/blob/main/examples/x_reflector_py_test.py).

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
