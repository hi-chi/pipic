/* Description: This file implements an extension for absorbing boundary condition for a particle-in-cell simulation. 
It allows particles to be added or removed based on their position relative to the boundaries of the simulation box and 
the mask function. The main purpose of this script is to show how to build extensions with piPIC it is there fore a downscaled 
version of the absorbing boundaries that comes with piPIC. */

#include "interfaces.h"
#include "ensemble.h"
#include "services.h"
#include <pybind11/pybind11.h>
#include "pybind11/stl.h"
#include <pybind11/operators.h>

const string name = "extension";
static double boundarySize;
static double density; // density of the particles
static int particleTypeIndex = 0; // type of particles to be removed 
static double fall = 1e-1; // shape parameter for the boundary
static int ax = 2; // axis along which the boundary is applied, 0 - x, 1 - y, 2 - z
static double gmin;
static double gmax;
static double temperature;
// Struct to hold the cell interface data.
static cellContainer ***cell;


// struct for managing the thread-specific data.
struct threadHandler{
    mt19937 rng;
    std::uniform_real_distribution<double> U1;
    std::normal_distribution<double> N1;
    threadHandler(): U1(0, 1.0), N1(0, 1.0) {}
    double random() {return U1(rng);} // returns a random number from [0, 1)
    double nrandom() {return N1(rng);} // returns a normal random number from [0, 1)
};
static vector<threadHandler> Thread;


// Function for adding new particles to the cell interface.
void addParticle(cellInterface &CI, particle &P){
    if(CI.particleBufferSize < CI.particleBufferCapacity){ // checking if the buffer permits adding a particle 
        *CI.newParticle(CI.particleBufferSize) = P; // copy particle to a new particle (buffer)
        CI.particleBufferSize++;
    } else {
        pipic_log.message("pi-PIC error: particle buffer overflow.", true);
    };
};

// Function for removing particles according to some probability rate.
void removeParticles(cellContainer* C, threadHandler &cthread, double rate){ 
    
    int particles = int(C->P.size());
    // iterate over particles in cell. The C->endShift particles from the end of the list 
    // came from other cells during the paricle push and should thus not be processed again.
    for (int ip = 0; ip < (particles - C->endShift); ip++){
        if (cthread.random() > rate) {
            C->P[ip].w = 0.; // set particle weight to zero to mark it for removal
        }
    }
};


// absorbing boundary
double mask(double x, double fall){
    return exp(fall*(cos(x*pi/2) - 1/cos(x*pi/2)));
};


// Field handler which is applied to each field grid node in the simulation box.
void fieldHandler(int* ind, double *r, double *E, double *B, double *dataDouble, int *dataInt){
    if(r[ax] < gmin + boundarySize){
        double rate = mask((-r[ax] + gmin + boundarySize)/boundarySize, fall);
        E[0] = rate*E[0];
        E[1] = rate*E[1];
        E[2] = rate*E[2];
        B[0] = rate*B[0];
        B[1] = rate*B[1];
        B[2] = rate*B[2];
    } else if (r[ax] > gmax - boundarySize){ 
        double rate = mask((r[ax] - gmax + boundarySize)/boundarySize, fall);
        E[0] = rate*E[0];
        E[1] = rate*E[1];
        E[2] = rate*E[2];
        B[0] = rate*B[0];
        B[1] = rate*B[1];
        B[2] = rate*B[2];
    };
};


// Particle handler which is applied to each cell in the simulation box.
void Handler(int *I, double *D, double *F, double *P, double *NP, double *dataDouble, int *dataInt){
    cellInterface CI(I, D, F, P, NP); // interface for manipulating with the content of a cell


    double cellMax[3] = {CI.cellMax().x, CI.cellMax().y, CI.cellMax().z};
    double cellMin[3] = {CI.cellMin().x, CI.cellMin().y, CI.cellMin().z};
    Thread.resize(omp_get_max_threads());

    // check if the cell is in the boundary region
    if ((cellMax[ax] < gmin + boundarySize) || (cellMin[ax] > gmax - boundarySize)) {   
        
        // thread safe handling
        threadHandler &cthread(Thread[CI.threadNum]); 
        cthread.rng.seed(CI.rngSeed); 

        // get rate depending on upper or lower boundary
        double rate;
        if (cellMax[ax] < gmin + boundarySize){
            rate = mask((-cellMax[ax]+ gmin + boundarySize)/boundarySize, fall);
        } else {
            rate = mask((cellMin[ax] - gmax + boundarySize)/boundarySize, fall);
        };

        // The handler can be applied to both cells and particles (depending on subject 
        // in the add_handler arguments). If the current cell interface points only to a cell then
        // CI.ParticleTypeIndex is -1,  otherwise it is the index of the particle type. Here we add particles during 
        // the handler call if CI.particleTypeIndex is -1, and remove particles if it is equal to particleTypeIndex.
        if (CI.particleTypeIndex==-1){ 
            // calculate number of particles to be generated in the cell
            double nb_particles = (1-rate)*density*CI.step.x*CI.step.y*CI.step.z; // number of particles to be generated in the cell
            double weight = nb_particles/10.; // ten particles per cell 

            for(int ip = 0; ip < 10; ip++){ 
                particle P;
                // generate position
                P.r.x = cellMin[0] + (cthread.random())*CI.step.x;
                P.r.y = cellMin[1] + (cthread.random())*CI.step.y;
                P.r.z = cellMin[2] + (cthread.random())*CI.step.z;
                        
                // generate momentum
                double3 p = {cthread.nrandom(), cthread.nrandom(), cthread.nrandom()};
                p = sqrt(electronMass*temperature)*p;
                p = sqrt(1 + 0.25*p.norm2()/sqr(electronMass*lightVelocity))*p;
                P.p = p;

                P.w = weight;
                P.id = particleTypeIndex; // set particle type index
                addParticle(CI, P);
            };

        } else if (CI.particleTypeIndex==particleTypeIndex){ // remove particles 
            // calculate global cell index
            int ig = CI.i.x + (CI.i.y + CI.i.z*CI.n.y)*CI.n.x; 
            // check if the cell is not empty
            if(cell[ig] != nullptr){
                // check if the cell contains particles of the type to be removed
                if(cell[ig][particleTypeIndex] != nullptr){
                    removeParticles(cell[ig][particleTypeIndex], cthread, rate); // remove particles according to the rate
                };
            };
        };
    };  
};


// particle handler initialization
int64_t handler(int64_t ensembleData, // adress of the ensemble data
                int64_t simbox,  // adress of the simulation box
                double _density, 
                double boundary_size,
                double _temperature){ 
    cell = (cellContainer***)ensembleData;
    boundarySize = boundary_size; // size of the boundary region
    density = _density; // density of the particles
    temperature = _temperature; // temperature of the particles
    gmin = ((simulationBox*)simbox)->min.z; // minimum z coordinate of the simulation box
    gmax = ((simulationBox*)simbox)->max.z; // maximum z coordinate of the simulation box
    ax = 2; // axis along which the boundary is applied, 0 - x, 1 - y, 2 - z
    Thread.resize(omp_get_max_threads());
    return (int64_t)Handler;
};

// field handler initialization
int64_t field_handler(int64_t simbox,
                      double boundary_size){ 
    gmin = ((simulationBox*)simbox)->min.z; // minimum z coordinate of the simulation box
    gmax = ((simulationBox*)simbox)->max.z; // maximum z coordinate of the simulation box
    ax = 2; // axis along which the boundary is applied, 0 - x, 1 - y, 2 - z
    boundarySize = boundary_size; // size of the boundary region
    return (int64_t)fieldHandler;
};

namespace py = pybind11;
PYBIND11_MODULE(extension, object) {
    object.attr("name") = name;
    object.def("handler", &handler, 
                py::arg("ensemble"), 
                py::arg("simulation_box"), 
                py::arg("density"),
                py::arg("boundary_size"),
                py::arg("temperature") = 0.0);
    object.def("field_handler", &field_handler,
                py::arg("simulation_box"),
                py::arg("boundary_size"));
}
