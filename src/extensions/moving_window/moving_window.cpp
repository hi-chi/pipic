#include "interfaces.h"
#include "services.h"
#include <pybind11/pybind11.h>
#include "pybind11/stl.h"
#include <pybind11/operators.h>

const string name = "moving_window";
static double _density, _temperature, _thickness; 
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

    double ts = CI.timeStep;
    double dx = CI.step.x;

    int rollback = floor(dataInt[0]*ts*lightVelocity/dx);
    int rollback_ = floor((dataInt[0]-1)*ts*lightVelocity/dx);

    if(rollback!=rollback_){        
        int nx = CI.n.x;
        double dy = CI.step.y;
        double dz = CI.step.z;
        double xmin = CI.globalMin.x; 
        double xmax = CI.globalMax.x;

        double r_rel = xmin + dx*(rollback%nx); 
        double r_min = r_rel - _thickness*dx;
        double r_max = r_rel;
        double3 cell_min = CI.cellMin();
        double3 cell_max = CI.cellMax();
        double eps = dx/10;

        if ((cell_min.x+eps >= r_min and cell_max.x-eps <= r_max) || 
                (cell_max.x+eps >= xmax - (xmin - r_min)) || 
                (cell_min.x-eps <= xmin + (r_max - xmax))){

            int ti = CI.particleTypeIndex;
            if (ti==0) {
                for(int ip = 0; ip < CI.particleSubsetSize; ip++){
                    CI.Particle(ip)->w = 0;
                }
            } else if (ti==-1){ 
                double nb_particles = _density*dx*dy*dz;
                double weight = nb_particles/(double)nbp;   
                double expectedNumber = _density*dx*dy*dz/weight;
                int numberToGenerate = int(expectedNumber) + (cthread.random() < (expectedNumber - int(expectedNumber)));
                if(numberToGenerate > 0){
                    for(int ip = 0; ip < numberToGenerate; ip++){
                        particle P;
                        // generate position
                        P.r.x = cell_min.x + (cthread.random())*dx;
                        P.r.y = cell_min.y + (cthread.random())*dy;
                        P.r.z = cell_min.z + (cthread.random())*dz;
                        
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
int64_t handler(double density, double thickness, int number_of_particles, double temperature){
    _density = density;
    _thickness = thickness;
    nbp = number_of_particles;
    _temperature = temperature;
    Thread.resize(omp_get_max_threads());
    return (int64_t)Handler;
};


namespace py = pybind11;
PYBIND11_MODULE(moving_window, object) {
    object.attr("name") = name;
    object.def("handler", &handler, py::arg("density"), py::arg("thickness"), py::arg("number_of_particles"), py::arg("temperature"));
}

