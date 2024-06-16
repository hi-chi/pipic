/*-------------------------------------------------------------------------------------------------------
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
// Description: A C/C++ extension that can be used to reflect particles from an x-limited region.

#include "interfaces.h"
#include <pybind11/pybind11.h>
#include "pybind11/stl.h"
#include <pybind11/operators.h>

const string name = "x_reflector_c";
static double xMin, xMax;

// function that is called to process particles in a given cell
void Handler(int *I, double *D, double *F, double *P, double *NP, double *dataDouble, int *dataInt){
    cellInterface CI(I, D, F, P, NP); // interface for manipulating with the content of a cell
    for(int ip = 0; ip < CI.particleSubsetSize; ip++)
        if((CI.Particle(ip)->r.x > xMin)&&(CI.Particle(ip)->r.x < xMax)&&(CI.Particle(ip)->p.x > 0))
            CI.Particle(ip)->p.x *= -1;
};

// extension initialization
int64_t handler(double location, double thickness){
    xMin = location - 0.5*thickness;
    xMax = location + 0.5*thickness;
    return (int64_t)Handler;
};

namespace py = pybind11;
PYBIND11_MODULE(_x_reflector_c, object) {
    object.attr("name") = name;
    object.def("handler", &handler, py::arg("location"), py::arg("thickness"));
}
