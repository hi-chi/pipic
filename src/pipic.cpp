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

PYBIND11_MODULE(_pipic, object) {
    object.attr("light_velocity") = lightVelocity;
    object.attr("electron_charge") = electronCharge;
    object.attr("electron_mass") = electronMass;
    object.attr("proton_mass") = protonMass;

    py::class_<pipic>(object, "init")
        .def(py::init<string, int, double, double, int, double, double, int, double, double>(),
        py::arg("solver"), py::arg("nx"), py::arg("xmin"), py::arg("xmax"), py::arg("ny")=1, py::arg("ymin")=-0.5, py::arg("ymax")=0.5, py::arg("nz")=1, py::arg("zmin")=-0.5, py::arg("zmax")=0.5)
        .def_readonly("nx", &pipic::nx)
        .def_readonly("xmin", &pipic::XMin)
        .def_readonly("xmax", &pipic::XMax)
        .def_readonly("ny", &pipic::ny)
        .def_readonly("ymin", &pipic::YMin)
        .def_readonly("ymax", &pipic::YMax)
        .def_readonly("nz", &pipic::nz)
        .def_readonly("zmin", &pipic::ZMin)
        .def_readonly("zmax", &pipic::ZMax)
        .def("get_number_of_particles", &pipic::getNumberOfParticles)
        .def("add_particles", &pipic::pyAddParticles, py::arg("name"), py::arg("number"), py::arg("charge"), py::arg("mass"), py::arg("temperature"), py::arg("density"), py::arg("data_double") = 0, py::arg("data_int") = 0)
        .def("particle_loop", &pipic::pyParticleLoop, py::arg("name"), py::arg("handler"), py::arg("data_double") = 0, py::arg("data_int") = 0)
        .def("field_loop", &pipic::pyFieldLoop, py::arg("handler"), py::arg("data_double") = 0, py::arg("data_int") = 0, py::arg("use_omp") = false)
        .def("custom_field_loop", &pipic::pyCustomFieldLoop, py::arg("number_of_iterations"), py::arg("it2r"), py::arg("field2data"), py::arg("data_double") = 0, py::arg("data_int") = 0)
        .def("advance", &pipic::pyAdvance, py::arg("time_step"), py::arg("number_of_iterations") = 1, py::arg("use_omp") = true)
        .def("fourier_solver_settings", &pipic::pyFourierSolverSettings, py::arg("divergence_cleaning") = -1, py::arg("sin2_kfilter") = -1)
        .def("log_policy", &pipic::pyLogPolicy, py::arg("log_to_file") = true, py::arg("log_to_screen") = false)
        .def("set_rng_seed", &pipic::setRngGenSeed, py::arg("seed"))
        .def("get_type_index", &pipic::getTypeIndex, py::arg("type_name"))
        .def("add_handler", &pipic::addHandler, py::arg("name"), py::arg("subject"), py::arg("handler"), py::arg("data_double") = 0, py::arg("data_int") = 0)
        .def("ensemble_data", &pipic::ensembleData)
        .def("en_corr_type", &pipic::en_corr_type, py::arg("correction_type") = 2)
    ;
}
