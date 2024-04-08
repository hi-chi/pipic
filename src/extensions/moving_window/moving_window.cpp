#include "interfaces.h"
#include "services.h"
#include <pybind11/pybind11.h>
#include "pybind11/stl.h"
#include <pybind11/operators.h>

const string name = "moving_window";
static double _temperature; 
static int _thickness; 
static int ppc;
static int64_t _density_profile;

void addParticle(cellInterface &CI, particle &P){
    if(CI.particleBufferSize < CI.particleBufferCapacity){ // checking if the buffer permits adding a particle 
        *CI.newParticle(CI.particleBufferSize) = P; // copy particle to a new particle (buffer)
        CI.particleBufferSize++;
    } else {
        pipic_log.message("pi-PIC error: particle buffer overflow.", true);
    };
}

struct threadHandler{
    mt19937 rng;
    std::uniform_real_distribution<double> U1;
    std::normal_distribution<double> N1;
    threadHandler(): U1(0, 1.0), N1(0, 1.0) {}
    double random() {return U1(rng);} // returns a random number from [0, 1)
    double nrandom() {return N1(rng);} // returns a normal random number from [0, 1)
};

static vector<threadHandler> Thread;

// function that is called to process particles in a given cell
void Handler(int *I, double *D, double *F, double *P, double *NP, double *dataDouble, int *dataInt){
    cellInterface CI(I, D, F, P, NP); // interface for manipulating with the content of a cell
    
    threadHandler &cthread(Thread[CI.threadNum]); // cthread (current thread) is to run in a thread-safe way
    cthread.rng.seed(CI.rngSeed);

    int rollback = floor(dataInt[0]*CI.timeStep*lightVelocity/CI.step.z);

    if(rollback%(_thickness/2)==0){        

        // initial r_rel is at the end of the cell and marks the end of the cleaning region 
        double r_rel = CI.globalMin.z + CI.step.z*(rollback%CI.n.z); 
        // r_min marks the beggining of the cleaning region
        double r_min = r_rel - _thickness*CI.step.z;
        double3 cell_min = CI.cellMin();
        double3 cell_max = CI.cellMax();
        double eps = CI.step.z/10;

        if ((cell_min.z+eps >= r_min and cell_max.z-eps <= r_rel) || 
            (cell_max.z+eps >= CI.globalMax.z - (CI.globalMin.z - r_min)) || 
            (cell_min.z-eps <= CI.globalMin.z + (r_rel - CI.globalMax.z))){
            
            if (CI.particleTypeIndex==0) {
                for(int ip = 0; ip < CI.particleSubsetSize; ip++){
                    CI.Particle(ip)->w = 0;
                }
            } else if (CI.particleTypeIndex==-1){
                double3 r = (cell_min + 0.5*CI.step);
                double R[3]; 
                R[0] = r.x;
                R[1] = r.y;
                
                // The position of the front of the window in 'real' coordinates
                double z_real = dataInt[0]*CI.timeStep*lightVelocity + CI.globalMax.z;
                if (r.z > r_rel){
                    R[2] = z_real - (r_rel - CI.globalMin.z) - (CI.globalMax.z - r.z);
                } else {
                    R[2] = z_real - (r_rel - r.z);
                };
                               
                double(*density_profile)(double*, double*, int*) = (double(*)(double*, double*, int*))_density_profile;
                double _density = density_profile(R, dataDouble, dataInt);

                double nb_particles = _density*CI.step.x*CI.step.y*CI.step.z;
                double weight = nb_particles/(double)ppc;   
                //double expectedNumber = nb_particles/weight;
                //int numberToGenerate = int(expectedNumber) + (cthread.random() < (expectedNumber - int(expectedNumber)));
                if(ppc > 0){
                    for(int ip = 0; ip < ppc; ip++){
                        particle P;
                        // generate position
                        P.r.x = cell_min.x + (cthread.random())*CI.step.x;
                        P.r.y = cell_min.y + (cthread.random())*CI.step.y;
                        P.r.z = cell_min.z + (cthread.random())*CI.step.z;
                        
                        // generate momentum
                        double3 p = {cthread.nrandom(), cthread.nrandom(), cthread.nrandom()};
                        p = sqrt(electronMass*_temperature)*p;
                        p = sqrt(1 + 0.25*p.norm2()/sqr(electronMass*lightVelocity))*p;
                        P.p = p;
                        
                        P.w = weight;
                        P.id = 0;
                        addParticle(CI, P);
                    };
                };
            };
        };
    };
};


// extension initialization
int64_t handler(int thickness, int particles_per_cell, double temperature, int64_t density){
    _thickness = thickness;
    ppc = particles_per_cell;
    _temperature = temperature;
    _density_profile = density;
    Thread.resize(omp_get_max_threads());
    return (int64_t)Handler;
};


namespace py = pybind11;
PYBIND11_MODULE(_moving_window, object) {
    object.attr("name") = name;
    object.def("handler", &handler, py::arg("thickness"), py::arg("particles_per_cell"), py::arg("temperature"),py::arg("density"));
}

