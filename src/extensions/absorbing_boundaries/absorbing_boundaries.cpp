#include "interfaces.h"
#include "ensemble.h"
#include "services.h"
#include <pybind11/pybind11.h>
#include "pybind11/stl.h"
#include <pybind11/operators.h>

const string name = "absorbing_boundaries";
static double boundarySize; 
static double _timeStep; // time step of the simulation
static double _temperature; // temperature of the particles
static double velocity; // velocity of the moving window
static int64_t densityProfile; // density of the particles
static double ppc; // particles per cell
static int ax;
static int mdir; // moving direction, 0 - x, 1 - y, 2 - z
static int n[3]; // number of cells in each direction
static double step[3]; // step in each direction
static double gmin[3]; // minimum coordinate in each direction
static double gmax[3]; // maximum coordinate in each direction
static double _fall; // fall rate for the field
static int pmf = 10; // particle removal frequency, every pmf-th time step particles are removed and replenished
bool staticsHasBeenSet = false; // flag to check if static variables have been set
simulationBox* _simbox; // cell interface for manipulating with the content of a cell
static cellContainer ***cell; // cellContainer is defined in ensemble.h, it contains particles of all types in a given cell

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

double mask(double x, double fall, double timeStep){
    return exp(fall*(cos(x*pi/2) - 1/cos(x*pi/2))*timeStep*1e16);
};

// add angle dependence
void fieldHandler(int* ind, double *r, double *E, double *B, double *dataDouble, int *dataInt){
    if(r[ax] < gmin[ax] + boundarySize){
        double rate = mask((-r[ax] + gmin[ax] + boundarySize)/boundarySize, _fall, _timeStep);
        E[0] = rate*E[0];
        E[1] = rate*E[1];
        E[2] = rate*E[2];
        B[0] = rate*B[0];
        B[1] = rate*B[1];
        B[2] = rate*B[2];
    } else if (r[ax] > gmax[ax] - boundarySize){ 
        double rate = mask((r[ax] - gmax[ax] + boundarySize)/boundarySize, _fall, _timeStep);
        E[0] = rate*E[0];
        E[1] = rate*E[1];
        E[2] = rate*E[2];
        B[0] = rate*B[0];
        B[1] = rate*B[1];
        B[2] = rate*B[2];
    };
};

void moving_r(double *r, int *dataInt, double timeStep){

    int rollback = floor(dataInt[0]*timeStep*velocity/step[mdir]);
    double r_rel = gmin[mdir] + step[mdir]*(rollback%n[mdir]); 

    // The position of the front of the window in 'real' coordinates
    double z_real = dataInt[0]*timeStep*velocity + gmax[mdir];// - gmax[ndir]*tan(_angle);
    if (r[mdir] > r_rel){
        r[mdir] = z_real - (r_rel - gmin[mdir]) - (gmax[mdir] - r[mdir]);
    } else {
        r[mdir] = z_real - (r_rel - r[mdir]);
    };
}

// function that is called to process particles in a given cell
void Handler(int *I, double *D, double *F, double *P, double *NP, double *dataDouble, int *dataInt){
    cellInterface CI(I, D, F, P, NP); // interface for manipulating with the content of a cell

    if (dataInt[0]%pmf==0){
        double cell_min[3] = {CI.cellMin().x, CI.cellMin().y, CI.cellMin().z};
        double cell_max[3] = {CI.cellMax().x, CI.cellMax().y, CI.cellMax().z};
        if (CI.particleTypeIndex==0){ // remove particles when electrons are called to count this time separately
            int ig, ix = CI.i.x, iy = CI.i.y, iz = CI.i.z; 
            ig = ix + (iy + iz*CI.n.y)*CI.n.x; 
            if ((cell_max[ax] < gmin[ax] + boundarySize) || (cell_min[ax] > gmax[ax] - boundarySize)) {

                double r[3] = {CI.cellMin().x + CI.step.x/2, 
                               CI.cellMin().y + CI.step.y/2, 
                               CI.cellMin().z + CI.step.z/2}; // position in the cell
                double rate;
                if (r[ax] < gmin[ax] + boundarySize){
                    rate = mask((-r[ax] + gmin[ax] + boundarySize)/boundarySize, _fall, CI.timeStep);
                } else {
                    rate = mask((r[ax] - gmax[ax] + boundarySize)/boundarySize, _fall, CI.timeStep);
                };                

                // removing particles randomly 
                if(cell[ig] != nullptr){
                    int it = 0;
                    if(cell[ig][it] != nullptr){
                        threadHandler &cthread(Thread[CI.threadNum]); // cthread (current thread) is to run in a thread-safe way
                        cthread.rng.seed(CI.rngSeed);

                        cellContainer *C = cell[ig][it];
                        int particlesToRemove = 0;
                        int particles = int(C->P.size());
                        // iterate over particles in cell
                        for (int ip = 0; ip < (particles - cell[ig][it]->endShift - particlesToRemove); ip++){
                            if (cthread.random() > rate) {
                                // move particle in the back of the list (just before processed particles and those that should be removed) to the place of this particle
                                C->P[ip] = C->P[particles - cell[ig][it]->endShift - 1 - particlesToRemove]; 
                                ip--; // adjust index 
                                particlesToRemove += 1; // increment the counter of particles to remove
                            };
                        };
                        // move processed particles to the positions of the particles to be removed
                        for (int ip = 0; ip < particlesToRemove; ip++){
                            if (cell[ig][it]->endShift > ip){
                                C->P[particles - cell[ig][it]->endShift - particlesToRemove + ip] = C->P[C->P.size()- 1 - ip]; // move last particle to current position
                            } else { 
                            break; 
                            };
                        };
                        // remove last particles
                        C->P.resize(C->P.size() - particlesToRemove); 
                    };
                };
            };
        } else if (CI.particleTypeIndex==-1 && densityProfile != -1){ // add particles when cell is called
            double r[3] = {CI.cellMin().x + CI.step.x/2, 
                           CI.cellMin().y + CI.step.y/2, 
                           CI.cellMin().z + CI.step.z/2}; // position in the cell
            if((r[ax] < gmin[ax] + boundarySize) || (r[ax] > gmax[ax] - boundarySize) ){ 
                threadHandler &cthread(Thread[CI.threadNum]); // cthread (current thread) is to run in a thread-safe way
                cthread.rng.seed(CI.rngSeed);
                double rate;


                if (r[ax] < gmin[ax] + boundarySize){
                    rate = mask((-r[ax] + gmin[ax] + boundarySize)/boundarySize, _fall, CI.timeStep);
                } else {
                    rate = mask((r[ax] - gmax[ax] + boundarySize)/boundarySize, _fall, CI.timeStep);
                };

                double(*getDensity)(double*, double*, int*) = (double(*)(double*, double*, int*))densityProfile;
                if (velocity != 0.0){
                    moving_r(r, dataInt, CI.timeStep); // adjust position of the particles according to the moving window
                };
                double density = getDensity(r, dataDouble, dataInt);

                double nb_particles = (1-rate)*density*step[0]*step[1]*step[2]; // number of particles to be generated in the cell
                int _ppc = int(ppc) + (cthread.random() < (ppc - int(ppc))); // particles per cell, it can be fractional
                double weight = nb_particles/double(_ppc); 

                if(_ppc > 0){
                    for(int ip = 0; ip < _ppc; ip++){ 
                        particle P;
                        // generate position
                        P.r.x = cell_min[0] + (cthread.random())*CI.step.x;
                        P.r.y = cell_min[1] + (cthread.random())*CI.step.y;
                        P.r.z = cell_min[2] + (cthread.random())*CI.step.z;

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


void set_statics(int64_t simbox, double boundary_size, char axis, double fall){
    _simbox = (simulationBox*)simbox;
    _fall = fall;

    if (axis == 'x') ax = 0;
    else if (axis == 'y') ax = 1;
    else if (axis == 'z') ax = 2;
    else {
        pipic_log.message("pi-PIC error: absorbing_boundaries set_statics(): axis must be 'x', 'y' or 'z'.", true);
    };

    n[0] = _simbox->n.x;
    n[1] = _simbox->n.y;
    n[2] = _simbox->n.z;
    step[0] = _simbox->step.x;
    step[1] = _simbox->step.y;
    step[2] = _simbox->step.z;
    gmin[0] = _simbox->min.x;
    gmin[1] = _simbox->min.y;
    gmin[2] = _simbox->min.z;
    gmax[0] = _simbox->max.x;
    gmax[1] = _simbox->max.y;
    gmax[2] = _simbox->max.z;

    if (boundary_size == -1.){
        boundarySize = (gmax[ax]-gmin[ax])/8; // default boundary size is equal to the cell size in x direction
        pipic_log.message("pi-PIC warning: absorbing_boundaries set_statics(): boundary_size is not set, using default value of (max-min)/8.", false);
    }else{
        boundarySize = boundary_size;
    }    
}



// extension initialization
int64_t handler(int64_t ensembleData, 
                int64_t simbox,  
                int64_t density_profile, 
                double boundary_size, 
                char axis, 
                double fall, 
                double temperature, 
                double particles_per_cell, 
                int remove_particles_every,
                double moving_window_velocity,
                char moving_window_direction){ 
    cell = (cellContainer***)ensembleData;
    _temperature = temperature;
    densityProfile = density_profile;
    ppc = particles_per_cell;
    pmf = remove_particles_every;
    velocity = moving_window_velocity;

    if (moving_window_direction == 'x') mdir = 0;
    else if (moving_window_direction == 'y') mdir = 1;
    else if (moving_window_direction == 'z') mdir = 2;
    else {
        pipic_log.message("pi-PIC error: absorbing_boundaries handler(): moving_window_direction must be 'x', 'y' or 'z'.", true);
    };

    if(!staticsHasBeenSet){
        set_statics(simbox,boundary_size, axis, fall);
        staticsHasBeenSet = true;
    };
    if (boundarySize <= 0){
        pipic_log.message("pi-PIC error: absorbing_boundaries handler(): boundary_size must larger than 0.", true);
    };
    Thread.resize(omp_get_max_threads());
    return (int64_t)Handler;
};


// extension initialization
int64_t field_handler(int64_t simbox, double timestep, double boundary_size, char axis, double fall){ 
    _timeStep = timestep;
    if(!staticsHasBeenSet){
        set_statics(simbox,boundary_size, axis, fall);
        staticsHasBeenSet = true;
    };
    if (boundarySize <= 0){
        pipic_log.message("pi-PIC error: absorbing_boundaries handler(): boundary_size must larger than 0.", true);
    };
    return (int64_t)fieldHandler;
};

namespace py = pybind11;
PYBIND11_MODULE(_absorbing_boundaries, object) {
    object.attr("name") = name;
    object.def("handler", &handler, 
               py::arg("ensemble"), 
               py::arg("simulation_box"), 
               py::arg("density_profile")=-1,
               py::arg("boundary_size")=-1.0,
               py::arg("axis")='x',
               py::arg("fall")=0.01,
               py::arg("temperature") = 0.0,
               py::arg("particles_per_cell") = 1.0,
               py::arg("remove_particles_every") = 10,
               py::arg("moving_window_velocity") = 0.0,
               py::arg("moving_window_direction") = 'x');

    object.def("field_handler", &field_handler, 
               py::arg("simulation_box"), 
               py::arg("timestep"),
               py::arg("boundary_size")=-1.0,
               py::arg("axis")='x',
               py::arg("fall")=0.01);
}
