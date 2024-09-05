/*-------------------------------------------------------------------------------------------------------
This file is an implementation of an extension landau_lifshitz, which is compartible with pi-PIC.
landau_lifshitz, Copyright 2023 Joel Magnusson
---------------------------------------------------------------------------------------------------------
landau_lifshitz is free software: you can redistribute it and/or modify it under the terms of the GNU General
Public License as published by the Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

landau_lifshitz is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with pi-PIC. If not, se
<https://www.gnu.org/licenses/>.
---------------------------------------------------------------------------------------------------------
Website: https://github.com/hi-chi/pipic
Contact: joel.magnusson@physics.gu.se.
-------------------------------------------------------------------------------------------------------*/
// Description: A C/C++ extension that can be used to account for Landau-Lifshitz radiation reaction.

#include "interfaces.h"
#include <pybind11/pybind11.h>
#include "pybind11/stl.h"
#include <pybind11/operators.h>

const string name = "landau_lifshitz";

// function that is called to process particles in a given cell
void Handler(int *I, double *D, double *F, double *P, double *NP, double *dataDouble, int *dataInt){
    cellInterface CI(I, D, F, P, NP); // interface for manipulating with the content of a cell
    for(int ip = 0; ip < CI.particleSubsetSize; ip++) {
        double c = lightVelocity;
        double m = CI.particleMass;
        double q = CI.particleCharge;
        double dt = CI.timeStep;
        double3 E, B;

        particle *P = CI.Particle(ip);
        CI.interpolateField(P->r, E, B);

        double gamma = sqrt(1 + sqr(P->p)/sqr(m*c));
        double3 v = P->p / (gamma*m);
        double c_inv = 1/c;

        double3 dp = dt * (2.0 / 3.0) * sqr(sqr(q) / (m * sqr(c))) *
            (cross(E, B) + c_inv * (cross(B, cross(B, v)) + dot(v, E) * E) -
                c_inv * sqr(gamma) * (sqr(E + c_inv * cross(v, B)) - sqr(dot(E, v) * c_inv)) * v);

        P->p += dp;
    }
};

// extension initialization
int64_t handler(){
    return (int64_t)Handler;
};

namespace py = pybind11;
PYBIND11_MODULE(_landau_lifshitz, object) {
    object.attr("name") = name;
    object.def("handler", &handler);
}
