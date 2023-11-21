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
    cout << "You are currently using version 1.0" << endl;
    cout << "Website: https://github.com/hi-chi/pipic" << endl;
    cout << "To notify about bugs and/or request functionality extensions contact arkady.gonoskov@gu.se." << endl;
    cout << "---------------------------------------------------------------------------------------------------------" << endl;
}

PYBIND11_MODULE(pipic, object) {
    object.def("info", &pipic_info, "Info about pi-PIC");
    object.attr("lightVelocity") = lightVelocity;
    object.attr("electronCharge") = electronCharge;
    object.attr("electronMass") = electronMass;
    object.attr("protonMass") = protonMass;

    py::class_<pipic>(object, "init")
        .def(py::init<string, int, double, double, int, double, double, int, double, double>(), 
        py::arg("solver"), py::arg("nx"), py::arg("XMin"), py::arg("XMax"), py::arg("ny")=1, py::arg("YMin")=-0.5, py::arg("YMax")=0.5, py::arg("nz")=1, py::arg("ZMin")=-0.5, py::arg("ZMax")=0.5)
        .def_readonly("nx", &pipic::nx)
        .def_readonly("XMin", &pipic::XMin)
        .def_readonly("XMax", &pipic::XMax)
        .def_readonly("ny", &pipic::ny)
        .def_readonly("YMin", &pipic::YMin)
        .def_readonly("YMax", &pipic::YMax)
        .def_readonly("nz", &pipic::nz)
        .def_readonly("ZMin", &pipic::ZMin)
        .def_readonly("ZMax", &pipic::ZMax)
        .def("getNumberOfParticles", &pipic::getNumberOfParticles)
        .def("addParticles", &pipic::pyAddParticles, py::arg("name"), py::arg("number"), py::arg("charge"), py::arg("mass"), py::arg("temperature"), py::arg("density"), py::arg("dataDouble") = 0, py::arg("dataInt") = 0)
        .def("particleLoop", &pipic::pyParticleLoop, py::arg("name"), py::arg("handler"), py::arg("dataDouble") = 0, py::arg("dataInt") = 0)
        .def("fieldLoop", &pipic::pyFieldLoop, py::arg("handler"), py::arg("dataDouble") = 0, py::arg("dataInt") = 0, py::arg("useOmp") = false)
        .def("customFieldLoop", &pipic::pyCustomFieldLoop, py::arg("numberOfIterations"), py::arg("it2r"), py::arg("field2data"), py::arg("dataDouble") = 0, py::arg("dataInt") = 0)
        .def("advance", &pipic::pyAdvance, py::arg("timeStep"), py::arg("numberOfIterations") = 1)
        .def("fourierSolverSettings", &pipic::pyFourierSolverSettings, py::arg("divergenceCleaning") = -1, py::arg("sin2_kFilter") = -1)
        .def("logPolicy", &pipic::pyLogPolicy, py::arg("logToFile") = true, py::arg("logToScreen") = false)
        .def("setRngGenSeed", &pipic::setRngGenSeed, py::arg("seed"))
        .def("getTypeIndex", &pipic::getTypeIndex, py::arg("typeName"))
        .def("addHandler", &pipic::addHandler, py::arg("name"), py::arg("subject"), py::arg("handler"), py::arg("dataDouble") = 0, py::arg("dataInt") = 0)
    ;
}