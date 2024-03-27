#include "interfaces.h"
#include "services.h"
#include <pybind11/pybind11.h>
#include "pybind11/stl.h"
#include <pybind11/operators.h>

const string name = "moving_window";
static double _temperature, _thickness; 
static int nbp;

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
    int rollback_ = floor((dataInt[0]-1)*CI.timeStep*lightVelocity/CI.step.z);

    if(rollback!=rollback_){        

        double r_rel = CI.globalMin.z + CI.step.z*(rollback%CI.n.z); 
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

                double nb_particles = dataDouble[0]*CI.step.x*CI.step.y*CI.step.z;
                double weight = nb_particles/(double)nbp;   
                double expectedNumber = dataDouble[0]*CI.step.x*CI.step.y*CI.step.z/weight;
                int numberToGenerate = int(expectedNumber) + (cthread.random() < (expectedNumber - int(expectedNumber)));
                if(numberToGenerate > 0){
                    for(int ip = 0; ip < numberToGenerate; ip++){
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
int64_t handler(double thickness, int number_of_particles, double temperature){
    //_density = density;
    _thickness = thickness;
    nbp = number_of_particles;
    _temperature = temperature;
    Thread.resize(omp_get_max_threads());
    return (int64_t)Handler;
};


namespace py = pybind11;
PYBIND11_MODULE(moving_window, object) {
    object.attr("name") = name;
    object.def("handler", &handler, py::arg("thickness"), py::arg("number_of_particles"), py::arg("temperature"));
}

