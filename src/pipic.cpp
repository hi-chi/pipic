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
// Description of the file: File introduces Python interfaces.

#include <algorithm>
#include <pybind11/pybind11.h>
#include "pybind11/stl.h"
#include <pybind11/operators.h>
#include "pipic.h"
#include <iostream>
#include <pthread.h>
#include <omp.h>
#include <iostream>
#include <vector>
#include <string>
#include <cstdint>
#include <cmath>
#include <pthread.h>
#include <omp.h>
#include <fftw3.h>
#include <complex>

namespace py = pybind11;

typedef complex<double> Complex;

void pipic_info()
{
    cout << "---------------------------------------------------------------------------------------------------------" << endl;
    cout << "pi-PIC, Copyright 2023 Arkady Gonoskov" << endl;
    cout << "---------------------------------------------------------------------------------------------------------" << endl;
    cout << "pi-PIC is free software: you can redistribute it and/or modify it under the terms of the GNU General" << endl;
    cout << "Public License as published by the Free Software Foundation, either version 3 of the License, or (at your" << endl;
    cout << "option) any later version." << endl;
    cout << endl;
    cout << "pi-PIC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even" << endl;
    cout << "the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the" << endl;
    cout << "GNU General Public License for more details." << endl;
    cout << endl;
    cout << "You should have received a copy of the GNU General Public License along with pi-PIC. If not, see"<< endl;
    cout << "<https://www.gnu.org/licenses/>." << endl;
    cout << "---------------------------------------------------------------------------------------------------------" << endl;
    cout << "You are currently using version 0.1." << endl;
    cout << "Website: https://github.com/hi-chi/pipic" << endl;
    cout << "To notify about bugs and/or request functionality extensions contact arkady.gonoskov@gu.se." << endl;
    cout << "---------------------------------------------------------------------------------------------------------" << endl;
}
 
PYBIND11_MODULE(pipic, object) {
    object.def("info", &pipic_info, "Info about piPIC");
    object.attr("lightVelocity") = lightVelocity;
    object.attr("electronCharge") = electronCharge;
    object.attr("electronMass") = electronMass;
    object.attr("protonMass") = protonMass;

    py::class_<pipic_boris>(object, "boris")
        .def(py::init<int, double, double, int, double, double, int, double, double>(), 
        py::arg("nx"), py::arg("XMin"), py::arg("XMax"), py::arg("ny")=1, py::arg("YMin")=-0.5, py::arg("YMax")=0.5, py::arg("nz")=1, py::arg("ZMin")=-0.5, py::arg("ZMax")=0.5)
        .def_readonly("nx", &pipic_boris::nx)
        .def_readonly("XMin", &pipic_boris::XMin)
        .def_readonly("XMax", &pipic_boris::XMax)
        .def_readonly("ny", &pipic_boris::ny)
        .def_readonly("YMin", &pipic_boris::YMin)
        .def_readonly("YMax", &pipic_boris::YMax)
        .def_readonly("nz", &pipic_boris::nz)
        .def_readonly("ZMin", &pipic_boris::ZMin)
        .def_readonly("ZMax", &pipic_boris::ZMax)
        .def("getNumberOfParticles", &pipic_boris::getNumberOfParticles)
        .def("addParticles", &pipic_boris::pyAddParticles, py::arg("name"), py::arg("number"), py::arg("charge"), py::arg("mass"), py::arg("temperature"), py::arg("density"), py::arg("dataDouble") = 0, py::arg("dataInt") = 0)
        .def("particleLoop", &pipic_boris::pyParticleLoop, py::arg("name"), py::arg("handler"), py::arg("dataDouble") = 0, py::arg("dataInt") = 0)
        .def("fieldLoop", &pipic_boris::pyFieldLoop, py::arg("handler"), py::arg("dataDouble") = 0, py::arg("dataInt") = 0, py::arg("useOmp") = false)
        .def("customFieldLoop", &pipic_boris::pyCustomFieldLoop, py::arg("numberOfIterations"), py::arg("it2r"), py::arg("field2data"), py::arg("dataDouble") = 0, py::arg("dataInt") = 0)
        .def("advance", &pipic_boris::pyAdvance, py::arg("timeStep"), py::arg("numberOfIterations") = 1)
        .def("setShufflingPolicy", &pipic_boris::pySetShufflingPolicy, py::arg("cellSuffleProbability") = 1, py::arg("doubleLoop") = true)
    ;

    py::class_<pipic_ec>(object, "ec")
        .def(py::init<int, double, double, int, double, double, int, double, double>(), 
        py::arg("nx"), py::arg("XMin"), py::arg("XMax"), py::arg("ny")=1, py::arg("YMin")=-0.5, py::arg("YMax")=0.5, py::arg("nz")=1, py::arg("ZMin")=-0.5, py::arg("ZMax")=0.5)
        .def_readonly("nx", &pipic_ec::nx)
        .def_readonly("XMin", &pipic_ec::XMin)
        .def_readonly("XMax", &pipic_ec::XMax)
        .def_readonly("ny", &pipic_ec::ny)
        .def_readonly("YMin", &pipic_ec::YMin)
        .def_readonly("YMax", &pipic_ec::YMax)
        .def_readonly("nz", &pipic_ec::nz)
        .def_readonly("ZMin", &pipic_ec::ZMin)
        .def_readonly("ZMax", &pipic_ec::ZMax)
        .def("getNumberOfParticles", &pipic_ec::getNumberOfParticles)
        .def("addParticles", &pipic_ec::pyAddParticles, py::arg("name"), py::arg("number"), py::arg("charge"), py::arg("mass"), py::arg("temperature"), py::arg("density"), py::arg("dataDouble") = 0, py::arg("dataInt") = 0)
        .def("particleLoop", &pipic_ec::pyParticleLoop, py::arg("name"), py::arg("handler"), py::arg("dataDouble") = 0, py::arg("dataInt") = 0)
        .def("fieldLoop", &pipic_ec::pyFieldLoop, py::arg("handler"), py::arg("dataDouble") = 0, py::arg("dataInt") = 0, py::arg("useOmp") = false)
        .def("customFieldLoop", &pipic_ec::pyCustomFieldLoop, py::arg("numberOfIterations"), py::arg("it2r"), py::arg("field2data"), py::arg("dataDouble") = 0, py::arg("dataInt") = 0)
        .def("advance", &pipic_ec::pyAdvance, py::arg("timeStep"), py::arg("numberOfIterations") = 1)
        .def("setShufflingPolicy", &pipic_boris::pySetShufflingPolicy, py::arg("cellSuffleProbability") = 1, py::arg("doubleLoop") = true)
        .def("logPolicy", &pipic_ec::pyLogPolicy, py::arg("logToFile") = true, py::arg("logToScreen") = false)
    ;
}
