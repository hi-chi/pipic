/*-------------------------------------------------------------------------------------------------------
This file is an implementation of an extension "focused_pulse" for pi-PIC.
focused_pulse, Copyright 2024 Christoffer Olofsson and Arkady Gonoskov
---------------------------------------------------------------------------------------------------------
focused_pulse is free software: you can redistribute it and/or modify it under the terms of the GNU General
Public License as published by the Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

focused_pulse is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with focused_pulse. If not, se
<https://www.gnu.org/licenses/>.
---------------------------------------------------------------------------------------------------------
Website: https://github.com/hi-chi/pipic
Contact: arkady.gonoskov@gu.se.
-------------------------------------------------------------------------------------------------------*/
// Description: This extension is for setting initial tightly focused pulses with arbitrary angle of
// focusing, transverse and longitudinal shapes, as well as phase front imperfections. The method 
// is a development of an approach described in Ref. 
// [E. Panova et al. Appl. Sci., 11(3), 956 (2021); https://doi.org/10.3390/app11030956].

#include "interfaces.h"
#include <pybind11/pybind11.h>
#include "pybind11/stl.h"
#include <pybind11/operators.h>

const string name = "focused_pulse";

struct pulse{
    int3 n;
    double3 min, max, center, path, e_axis;
    double theta_max, l_size;
    double(*shape)(double*);
    bool checklist[9];
    bool readyToGo;
    double dist;
    double3 step, npath, size_; // step, normalized path, inverse size of the box
    pulse(): readyToGo(false) {
        for(int i = 0; i < 9; i++) checklist[i] = false;
    }
    void singlePointSet(int3 i, double *E, double *B)
    {
        double r[3];
        r[0] = min.x + i.x*step.x;
        r[1] = min.y + i.y*step.y;
        r[2] = min.z + i.z*step.z;
        double R = sqrt(sqr(r[0] - center.x) + sqr(r[1] - center.y) + sqr(r[2] - center.z));
        if(R == 0) return;
        double rd = dist - R; // relative distance
        if(abs(rd) > 0.5*l_size) return;
        double3 r0(r[0] - center.x, r[1] - center.y, r[2] - center.z);
        r0 = (1/R)*r0;
        double3 b0 = cross(r0, e_axis);
        b0.normalize();
        double3 e0 = cross(r0, b0);
    
        double par[4]; // relative distance, theta, alpha, beta = acos(dot(e_axis, r0))
        par[0] = rd;
        double r0npath_dot = dot(r0, npath);
        if(1 - abs(r0npath_dot) > 1e-10) par[1] = acos(-r0npath_dot); else par[1] = pi*(r0npath_dot > 0);
        if(par[1] > theta_max) return;
        double3 r1 = r0 - dot(r0, npath)*npath;
        r1.normalize();
        double r1e0_dot = dot(r1, e0);
        if(1 - abs(r1e0_dot) > 1e-10) par[2] = acos(r1e0_dot); else par[1] = pi*(r1e0_dot < 0);
        if(dot(cross(e_axis, npath), r0) < 0) par[2] = 2*pi - par[2];
        
        double e_axis_r0_dot = dot(e_axis, r0);
        if(1 - abs(e_axis_r0_dot) > 1e-10) par[3] = acos(e_axis_r0_dot); else par[3] = pi*(e_axis_r0_dot < 0);
        
        double f = shape(&par[0]) * (dist/R);

        E[0] += f*e0.x;
        E[1] += f*e0.y;
        E[2] += f*e0.z;
        
        B[0] += f*b0.x;
        B[1] += f*b0.y;
        B[2] += f*b0.z;
    }
    void setField(int *i, double *E, double *B){
        if(!readyToGo){ // check if all the required parameters are set
            for(int i = 0; i < 9; i++) if(!checklist[i]){ 
                cout << "pipic::focused_pulse: Error: Not all required parameters are provided." << endl;
                cout << "check list: "; for(int j = 0; j < 9; j++) cout << checklist[j] << ":";
                exit(0);
            }
            readyToGo = true;
        }
        //int3 ind(i[0], i[1], i[2]);
        //singlePointSet(ind, E, B);
        
        int3 hi_min, hi_max; // limits for hyper index
        hi_min.x = floor((center.x - dist - 0.5*l_size - min.x)*size_.x);
        hi_min.y = floor((center.y - dist - 0.5*l_size - min.y)*size_.y);
        hi_min.z = floor((center.z - dist - 0.5*l_size - min.z)*size_.z);
        hi_max.x = floor((center.x + dist + 0.5*l_size - min.x)*size_.x) + 1;
        hi_max.y = floor((center.y + dist + 0.5*l_size - min.y)*size_.y) + 1;
        hi_max.z = floor((center.z + dist + 0.5*l_size - min.z)*size_.z) + 1;
        
        //full loop over all hypercells (~(dist/l_size)^3):
        //for(int hix = hi_min.x; hix <= hi_max.x; hix += 1)
        //for(int hiy = hi_min.y; hiy <= hi_max.y; hiy += 1)
        //for(int hiz = hi_min.z; hiz <= hi_max.z; hiz += 1)
        //    singlePointSet({hix*n.x + i[0], hiy*n.y + i[1], hiz*n.z + i[2]}, E, B);
        //return;

        //optimized loop (~(dist/l_size)^2):
        for(int hiy = hi_min.y; hiy <= hi_max.y; hiy += 1)
        for(int hiz = hi_min.z; hiz <= hi_max.z; hiz += 1)
        {
            double rsx = min.x + i[0]*step.x - center.x;
            double rsy = min.y + i[1]*step.y - center.y;
            double rsz = min.z + i[2]*step.z - center.z;

            double rt2 = sqr(hiy*(max.y - min.y) + rsy) + sqr(hiz*(max.z - min.z) + rsz);
            double rmin2 = sqr(dist - 0.5*l_size);
            double rmax2 = sqr(dist + 0.5*l_size);
            int hiMin, hiMax;

            double rxmin = -1, rxmax = -1; 

            if(rt2 >= rmin2) hiMin = (rsx <= 0); 
            else {
                rxmin = sqrt(rmin2 - rt2);
                hiMin = floor((rxmin - rsx)*size_.x) + 1;
            }
            if(rt2 > rmax2) hiMax = -1;
            else {
                rxmax = sqrt(rmax2 - rt2);
                hiMax = floor((rxmax - rsx)*size_.x);
            }
            for(int hix = hiMin; hix <= hiMax; hix += 1)
                singlePointSet({hix*n.x + i[0], hiy*n.y + i[1], hiz*n.z + i[2]}, E, B);

            if(rt2 >= rmin2) hiMin = (rsx > 0);
            else hiMin = floor((rxmin + rsx)*size_.x) + 1;
            if(rt2 > rmax2) hiMax = -1;
            else hiMax = floor((rxmax + rsx)*size_.x);

            for(int hix = hiMin; hix <= hiMax; hix += 1)
                singlePointSet({-hix*n.x + i[0], hiy*n.y + i[1], hiz*n.z + i[2]}, E, B);   
        }
    }
};
static vector<pulse> Pulse(1);

void add_pulse(){ // a function to add a pulse and start configuring it
    Pulse.push_back(pulse());
};

size_t get_pulse_count() { //Added method to retrieve the number of pulses placed in the domain. Solely for debug and flexibility.
    return Pulse.size();
};

void clear_pulse() { //Added method to clear entire pulse vector, for reusability purposes.
    while (Pulse.size() > 0) Pulse.pop_back();
};

void set_box(int nx, int ny, int nz, double xmin, double ymin, double zmin, double xmax, double ymax, double zmax){
    int cp = Pulse.size() - 1; // current pulse to set
    Pulse[cp].n = int3(nx, ny, nz); Pulse[cp].checklist[0] = true;
    Pulse[cp].min = double3(xmin, ymin, zmin); Pulse[cp].checklist[1] = true;
    Pulse[cp].max = double3(xmax, ymax, zmax); Pulse[cp].checklist[2] = true;
    Pulse[cp].size_ = double3(1/(xmax - xmin), 1/(ymax - ymin), 1/(zmax - zmin));
    Pulse[cp].step = double3((xmax - xmin)/nx, (ymax - ymin)/ny, (zmax - zmin)/nz);
    Pulse[cp].center = double3(0.5*(xmin + xmax), 0.5*(ymin + ymax), 0.5*(zmin + zmax)); Pulse[cp].checklist[3] = true;
    Pulse[cp].theta_max = pi; Pulse[cp].checklist[6] = true; // default value
};

void set_center(double xcenter, double ycenter, double zcenter){
    Pulse[Pulse.size() - 1].center = double3(xcenter, ycenter, zcenter); Pulse[Pulse.size() - 1].checklist[3] = true;
};

void set_path(double xpath, double ypath, double zpath){
    Pulse[Pulse.size() - 1].path = double3(xpath, ypath, zpath); Pulse[Pulse.size() - 1].checklist[4] = true;
    Pulse[Pulse.size() - 1].dist = Pulse[Pulse.size() - 1].path.norm();
    Pulse[Pulse.size() - 1].npath = Pulse[Pulse.size() - 1].path;
    Pulse[Pulse.size() - 1].npath.normalize();
};

void set_e_axis(double x_e_axis, double y_e_axis, double z_e_axis){
    Pulse[Pulse.size() - 1].e_axis = double3(x_e_axis, y_e_axis, z_e_axis); Pulse[Pulse.size() - 1].checklist[5] = true;
    Pulse[Pulse.size() - 1].e_axis.normalize();
};

void set_theta_max(double theta_max){
    Pulse[Pulse.size() - 1].theta_max = theta_max; Pulse[Pulse.size() - 1].checklist[6] = true;
};

void set_l_size(double l_size){
    Pulse[Pulse.size() - 1].l_size = l_size; Pulse[Pulse.size() - 1].checklist[7] = true;
};

void set_shape(int64_t shape){
    Pulse[Pulse.size() - 1].shape = (double(*)(double *))shape; Pulse[Pulse.size() - 1].checklist[8] = true;
};

void Field_loop_cb(int *ind, double *r, double *E, double *B, double *dataDouble, double *dataInt){
    for(int i = 0; i < int(Pulse.size()); i++) Pulse[i].setField(ind, E, B);
};

// return the callback
int64_t field_loop_cb(){
    return (int64_t)Field_loop_cb;
};

namespace py = pybind11;
PYBIND11_MODULE(_focused_pulse, object) {
    object.attr("name") = name;
    object.def("set_box", &set_box, py::arg("nx"), py::arg("ny"), py::arg("nz"), py::arg("xmin"), py::arg("ymin"), py::arg("zmin"), py::arg("xmax"), py::arg("ymax"), py::arg("zmax"));
    object.def("set_center", &set_center, py::arg("x"), py::arg("y"), py::arg("z"));
    object.def("set_path", &set_path, py::arg("x"), py::arg("y"), py::arg("z"));
    object.def("set_e_axis", &set_e_axis, py::arg("x"), py::arg("y"), py::arg("z"));
    object.def("set_theta_max", &set_theta_max, py::arg("theta_max"));
    object.def("set_l_size", &set_l_size, py::arg("l_size"));
    object.def("set_shape", &set_shape, py::arg("shape"));
    object.def("add_pulse", &add_pulse);
    object.def("get_pulse_count", &get_pulse_count);
    object.def("clear_pulse", &clear_pulse);
    object.def("field_loop_cb", &field_loop_cb);
}
