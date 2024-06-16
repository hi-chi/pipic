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
// Description: A C/C++ extension that can be used to convert type of particles in an x-limited region.

#include "interfaces.h"
#include <pybind11/pybind11.h>
#include "pybind11/stl.h"
#include <pybind11/operators.h>

const string name = "x_converter_c";
static double xMin, xMax; // limits of the action
static int typeTo; // the typeIndex to set for affected particles

// function that is called to process particles in a given cell
void Handler(int *I, double *D, double *F, double *P, double *NP, double *dataDouble, int *dataInt){
    cellInterface CI(I, D, F, P, NP);
    // changing type requires replacing existing particle with a new particle with altered typeIndex
    for(int ip = 0; ip < CI.particleSubsetSize; ip++)
        if((CI.Particle(ip)->r.x > xMin)&&(CI.Particle(ip)->r.x < xMax))
            if(CI.particleBufferSize < CI.particleBufferCapacity){ // checking if the buffer permits adding a particle
                *CI.newParticle(CI.particleBufferSize) = *CI.Particle(ip); // copy particle to a new particle (buffer)
                CI.newParticle(CI.particleBufferSize)->id = typeTo; // put the new particle type in id (convention)
                CI.Particle(ip)->w = 0; // zero weight indicates that the particle has to be removed
                CI.particleBufferSize++;
            }

};

// extension initialization
int64_t handler(double location, double thickness, int TypeTo){
    xMin = location - 0.5*thickness;
    xMax = location + 0.5*thickness;
    typeTo = TypeTo;
    return (int64_t)Handler;
};

namespace py = pybind11;
PYBIND11_MODULE(_x_converter_c, object) {
    object.attr("name") = name;
    object.def("handler", &handler, py::arg("location"), py::arg("thickness"), py::arg("typeTo"));
}
