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
// Description of the file: Implementation of an electrostatic 1D solver.


#include "interfaces.h"

struct ES1DFieldSolver: public field_solver // spectral solver for electromagnetic field evolution with zero current
{
    string type;
    vector<double> Ex; // electric field at the nodes
    vector<double> Jx; // current density at the nodes
    ES1DFieldSolver(simulationBox box): field_solver(box){
        type = "ES1D";
        Ex.resize(box.n.x, 0.0); // initialize electric field
        Jx.resize(box.n.x, 0.0); // initialize current density
    }

    void advance(double timeStep){
        //advance Ex
        for(size_t ix = 0; ix < Ex.size(); ix++){
            Ex[ix] += -timeStep*4*M_PI*Jx[ix];
        }
    }

    // CIC weighting function (not mandatory method for field solver interface)
    double CIC(double x, int i, double dx, double xmin){
        return ((i + 0.5)*dx + xmin - (x - dx/2))/dx;
    }

    void getEB(double3 r, double3 &E, double3 &B)
    {
        // get electric field at the point r
        int ix = int((r.x - box.min.x)/box.step.x);
        double w_left = CIC(r.x, ix, box.step.x, box.min.x);
        double w_right = 1 - w_left;

        E.x = w_left * Ex[ix];
        E.x += w_right * Ex[(ix + 1)%box.n.x];
        E.y = 0; 
        E.z = 0;
        B = {0, 0, 0}; // magnetic field is zero in this solver
    }

    // interface for setting/modifying field state:
    // makes a loop over all nodes and calls a function handler(ind, r, E, B, dataDouble_, dataInt_), which takes:
    // three-index (ind[0], ind[1], ind[2]),
    // field component code: ind[4] = 0 for Ex, 1 for Ey, 2 for Ez, 3 for Bx, 4 for By, 5 for Bz, 6 for all)
    // coordinate of the node in question,
    // field values (whatever is applicable according to ind[4]), and refernces to data of double and int type.

    void fieldLoop(int64_t handler, int64_t dataDouble = 0, int64_t dataInt = 0, bool useOmp = false){
        void(*handler_)(int*, double*, double*, double*, double*, int*) = (void(*)(int*, double*, double*, double*, double*, int*))handler;
        double* dataDouble_ = nullptr; if(dataDouble != 0) dataDouble_ = (double*)dataDouble;
        int* dataInt_ = nullptr; if(dataInt != 0) dataInt_ = (int*)dataInt;
        if(useOmp){
            #pragma omp parallel for collapse(3)
            for(int iz = 0; iz < box.n.z; iz++)
            for(int iy = 0; iy < box.n.y; iy++)
            for(int ix = 0; ix < box.n.x; ix++)
                fieldLoopBody({ix, iy, iz}, handler_, dataDouble_, dataInt_);
        } else {
            for(int iz = 0; iz < box.n.z; iz++)
            for(int iy = 0; iy < box.n.y; iy++)
            for(int ix = 0; ix < box.n.x; ix++)
                fieldLoopBody({ix, iy, iz}, handler_, dataDouble_, dataInt_);
        }
    }

    void fieldLoopBody(int3 i, void(*handler_)(int*, double*, double*, double*, double*, int*), double* dataDouble_, int* dataInt_){
        double3 r3 = box.min;
        r3.x += i.x*box.step.x;
        r3.y += i.y*box.step.y;
        r3.z += i.z*box.step.z;

        double r[3] = {r3.x, r3.y, r3.z};
        int ind[4] = {i.x, i.y, i.z, 6};
        double E[3], B[3];
        E[0] = Ex[i.x]; B[0] = 0;
        E[1] = 0; B[1] = 0;
        E[2] = 0; B[2] = 0;
        handler_(ind, r, E, B, dataDouble_, dataInt_);
        Ex[i.x] = E[0];
    }

    
    //read-only interface for accessing fields in a set of points:
    void customFieldLoop(int numberOfIterations, int64_t it2coord, int64_t field2data, int64_t dataDouble = 0, int64_t dataInt = 0){
        //customFieldLoop_via_getEB<fourierSolver>(this, numberOfIterations, it2coord, field2data, dataDouble, dataInt);
    } 
      
    
    ~ES1DFieldSolver(){}

    
    inline void cellSetField(cellInterface &CI, int3 i){
        setGridType(CI, 0); 

        // F_data is an array containing the field at the nodes surrounding the cell
        // the structure is [Ex1,Ey1,Ez1,Bx1,By1,Bz1, Ex2,Ey2,Ez2,Bx2,By2,Bz2,...,Ex8,Ey8,Ez8,Bx8,By8,Bz8]
        // where 1-8 is the enumeration of the nodes. 
        if(getCI_F_Data(CI) == nullptr) getCI_F_Data(CI) = new double[48];
        double* F_data = getCI_F_Data(CI);

        int cig[2] = {i.x,i.x+1};
        if (cig[1] > box.n.x){cig[1] = 0;}

        for(int j = 0; j < (1 << box.dim); j++){
            F_data[6*j] = Ex[cig[j]];
        }
    }
};