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
// Description: Here main base classes and interfaces are introduced.

#ifndef INTERFACES_H
#define INTERFACES_H

#include "primitives.h"

struct particle 
{
    double3 r, p; // coordinate and momentum
    double w; // weight
    unsigned long long int id; // id 
    // If necessary add variables to allocate memory for additional attributes (must be of type double).
    // The code must be rebuilt for a specific number of attributes.
    // Use of attributes should be coordinated between extensions.
};

struct particleType
{
    string name;
    double charge, mass;
    particleType(string name, double charge, double mass) : name(name), charge(charge), mass(mass) {} 
};

struct ensemble; // forward declaration of data storying structures
struct field_solver;
struct cellHandler;

struct cellInterface // main structure for developing extensions; provides access to particles and/or field values for all cells
{
    //access to basic parameters:
    const int3 &i, &n; // three-index of the cell and three-size of the grid
    const int &dim, &aNum; // dimensionality and the number of additional attributes of particles
    const int &PType, &PSize; // particle type being processed and their number in the cell being processed
    const int &NPcapacity; // capacity of the buffer for new particles
    int &NPSize; // a variable to indicate the number of particles to be added (must stay <= NPcapacity)
    const int &gridType; // indicates the type of the grid and defines the way of interpolation
    const int &threadNum; // indicates the number of active thread, use 1 to avoid clashes (0 sometimes is not called)
    const double3 &globalMin, &globalMax; // the limits of computational region and cell being processed
    const double3 &step, &invStep; // step size and its inverse;
    const double &timeStep; // time step;
    const double &PCharge, &PMass; // charge and mass of particles being processed
    const double3 cellMin(){return {globalMin.x + i.x*step.x, globalMin.y + i.y*step.y, globalMin.z + i.z*step.z};}
    const double3 cellMax(){return {globalMin.x + (i.x+1)*step.x, globalMin.y + (i.y+1)*step.y, globalMin.z + (i.z+1)*step.z};}
   
    //access to particles being processed and new to be added (indices must be within limits):
    particle* P(int i){return (particle*)(P_data + i*(8 + I[7]));} // pointer to i-th particle to be processed (use of attributes requires manual compartibility control)
    //Warning: It is foriden to change particles' coordinates P(i)->r
    //To remove a particle set its weight P(i).w to zero.
    particle* NP(int i){return (particle*)(NP_data + i*(8 + I[7]));} // pointer to i-th particle in the buffer for new particles
    //The type is to be placed into id (by convention).
    double& a(particle *P, int ia){ return *((double*)P + 8 + ia);} // access to ia-th attribute; ia must be < aNum

    //access to field:
    void interpolateField(double3 r, double3 &E, double3 &B); // interpolates E and B field at point r 
    
    // initialization:
    cellInterface(int *I, double *D, double *F, double *P, double *NP):
    I(I), D(D), F_data(F), P_data(P), NP_data(NP),
    i(*((int3*)I)), n(*((int3*)I + 1)), 
    dim(I[6]), aNum(I[7]), PType(I[8]), PSize(I[9]), NPcapacity(I[10]), NPSize(I[11]), 
    gridType(I[12]), threadNum(I[13]),
    globalMin(*(double3*)(D)), globalMax(*(double3*)(D + 3)), 
    step(*(double3*)(D + 6)), invStep(*(double3*)(D + 9)),
    timeStep(D[12]), PCharge(D[13]), PMass(D[14])
    {}
    private:
    // data (to be accessed via interfaces):
    int *I; // integer parameters
    double *D; // double parameters
    double *F_data; // em-field at the coners of the cell
    double *P_data; // data for storing particles at the cell of the type being processed
    double *NP_data; // data for the buffer array of particles to be added    
    
    friend struct ensemble;
    friend struct field_solver;
    friend struct cellHandler;
};

struct simulationBox
{
    double3 min, max, size, invSize, step, invStep; // physical limits of the computational region, auxiliary vectors
    int3 n; // number of cells alongs each dimension; must be powers of two; for 2D set n.z = 1, for 1D set n.z = n.y = 1 
    unsigned int ng; // total number of cells
    int dim; // problem dimensionality
    simulationBox(int3 n, double3 min, double3 max): n(n), min(min), max(max), ng(((unsigned int)n.x)*n.y*n.z), 
    size(max.x - min.x, max.y - min.y, max.z - min.z),
    invStep(n.x/(max.x - min.x), n.y/(max.y - min.y), n.z/(max.z - min.z)),
    invSize(1/(max.x - min.x), 1/(max.y - min.y), 1/(max.z - min.z)),
    step((max.x - min.x)/n.x, (max.y - min.y)/n.y, (max.z - min.z)/n.z)
    {
        dim = 3;
        if(n.z == 1) dim = 2; 
        if((n.y == 1)&&(n.z == 1)) dim = 1; 
    }
    unsigned int ig(int3 i) // conversion from three-index to global index
    {
        return i.x + (i.y + ((unsigned int)i.z)*n.y)*n.x; // must be consistent with fftw plans
    }
};

struct field_solver
{
    string type;
    simulationBox box;
    field_solver(simulationBox box): box(box){}

    virtual void advance(double timeStep) = 0;

    // interface for setting/modifying field state:
    // makes a loop over all nodes and calls a function handler(ind, r, E, B, dataDouble_, dataInt_), which takes:
    // three-index (ind[0], ind[1], ind[2]), 
    // field component code: ind[4] = 0 for Ex, 1 for Ey, 2 for Ez, 3 for Bx, 4 for By, 5 for Bz, 6 for all)
    // coordinate of the node in question,
    // field values (whatever is applicable according to ind[4]), and frefernces to data of double and int type.
    virtual void fieldLoop(int64_t handler, int64_t dataDouble = 0, int64_t dataInt = 0, bool useOmp = false) = 0;
    //read-only interface for accessing fields in a set of points:
    virtual void customFieldLoop(int numberOfIterations, int64_t it2coord, int64_t field2data, int64_t dataDouble = 0, int64_t dataInt = 0) = 0;
    
    //functions neded for cellInterface
    inline void setGridType(cellInterface &CI, int gridType){CI.I[12] = gridType;};
    inline double*& getCI_F_Data(cellInterface &CI){return CI.F_data;};
};

struct pic_solver
{
    field_solver *Field;
    ensemble *Ensemble;
    virtual void advance(double timeStep) = 0;
    //optional functions:
    void preLoop(){};
    void postLoop(){};
    virtual ~pic_solver(){}
};

void cellInterface::interpolateField(double3 r, double3 &E, double3 &B){
    if(gridType == 0){ // collocated grid: E and B are defined at the corners of each cell
        double c[8];
        double3 cmin = cellMin();
        c[1] = (r.x - cmin.x)*invStep.x; c[0] = 1 - c[1];
        if(dim > 1){
            double wy = (r.y - cmin.y)*invStep.y; double wy_ = 1 - wy;
            for(int i = 0; i < 2; i++) c[i+2] = c[i]*wy;
            for(int i = 0; i < 2; i++) c[i] *= wy_;
        }
        if(dim > 2){
            double wz = (r.z - cmin.z)*invStep.z; double wz_ = 1 - wz;
            for(int i = 0; i < 4; i++) c[i+4] = c[i]*wz;
            for(int i = 0; i < 4; i++) c[i] *= wz_;
        }
        E = {0, 0, 0}; B = {0, 0, 0};
        for(int i = 0; i < (1 << dim); i++){
            E.x += c[i]*F_data[6*i];
            E.y += c[i]*F_data[6*i+1];
            E.z += c[i]*F_data[6*i+2];
            B.x += c[i]*F_data[6*i+3];
            B.y += c[i]*F_data[6*i+4];
            B.z += c[i]*F_data[6*i+5];
        }
        return;
    }
    cout << "pi-PIC cellInterface error: unknown gridType(" << gridType << "). Update interface.h" << endl;
};

#endif