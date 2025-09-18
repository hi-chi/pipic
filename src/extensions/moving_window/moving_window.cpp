#include "interfaces.h"
#include "services.h"
#include <pybind11/pybind11.h>
#include "pybind11/stl.h"
#include <pybind11/operators.h>

const string name = "moving_window";
static double _temperature; 
static int _thickness; 
static double ppc;
static int64_t _density_profile;
static double _velocity; // velocity of the moving window
static double _angle; // angle of the moving window
static double _timeStep; // time step of the simulation
static int dir; // principal direction of the moving window
static int ndir;
static int n[3]; // number of cells in each direction
static double step[3]; // step in each direction
static double gmin[3]; // minimum coordinate in each direction
static double gmax[3]; // maximum coordinate in each direction
bool staticsHasBeenSet = false; // flag to check if static variables have been set
simulationBox* _simbox; // cell interface for manipulating with the content of a cell
double *_time; // current simulation time

struct cellContainer
{
    vector<particle> P;
    int endShift; 
    cellContainer(): endShift(0) {}
};
static cellContainer ***cell;

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


// add angle dependence
void fieldHandler(int* ind, double *r, double *E, double *B, double *dataDouble, int *dataInt){
    int rollback = floor(_time[0] *_velocity/step[dir]);
    if(rollback%(_thickness/2)==0 && rollback > 0){ 
        int rollback_prev = floor((_time[0] - _timeStep)*_velocity/step[dir]);
        if(rollback_prev!=rollback){
            // initial r_rel is at the end of the cell and marks the end of the cleaning region 
            double r_rel = gmin[dir] + step[dir]*(rollback%n[dir]); 
            // r_min marks the beggining of the cleaning region
            double r_min = r_rel - _thickness*step[dir];

            r_rel -= (r[ndir])*tan(_angle); 
            r_min -= (r[ndir])*tan(_angle); 

            // if r_min or r_rel is smaller than z_min we move it to the other side of the window
            if (r_min < gmin[dir]) r_min = gmax[dir] - (gmin[dir] - r_min);
            if (r_rel < gmin[dir]) r_rel = gmax[dir] - (gmin[dir] - r_rel);

            double eps = 0;
            if ((r[dir]+eps >= r_min and r[dir]-eps <= r_rel) ||
                ( (r_min > r_rel) and ((r[dir] + eps >= r_min) || (r[dir] - eps <= r_rel)))){ 
                E[0] = 0;
                E[1] = 0;
                E[2] = 0;
                B[0] = 0;
                B[1] = 0;
                B[2] = 0;
            }
        }
    }
}


// function that is called to process particles in a given cell
void Handler(int *I, double *D, double *F, double *P, double *NP, double *dataDouble, int *dataInt){
    cellInterface CI(I, D, F, P, NP); // interface for manipulating with the content of a cell
    int rollback = floor(_time[0] * _velocity/step[dir]);
    if(rollback%(_thickness/2)==0 && rollback > 0){ 
        int rollback_prev = floor((_time[0]-_timeStep)*_velocity/step[dir]);
        if(rollback_prev!=rollback){
            double cell_min[3] = {CI.cellMin().x, CI.cellMin().y, CI.cellMin().z};
            double cell_max[3] = {CI.cellMax().x, CI.cellMax().y, CI.cellMax().z};


            // initial r_rel is at the end of the cell and marks the end of the cleaning region 
            double r_rel = gmin[dir] + step[dir]*(rollback%n[dir]); 
            // r_min marks the beggining of the cleaning region
            double r_min = r_rel - _thickness*step[dir];

            r_rel -= ( cell_min[ndir])*tan(_angle); 
            r_min -= ( cell_min[ndir])*tan(_angle); 

            // if r_min or r_rel is smaller than z_min we move it to the other side of the window
            if (r_min < gmin[dir]) r_min = gmax[dir] - (gmin[dir] - r_min);
            if (r_rel < gmin[dir]) r_rel = gmax[dir] - (gmin[dir] - r_rel);

            // if r_min or r_rel is larger than z_max we move it to the other side of the window
            if (r_min > gmax[dir]) r_min = gmin[dir] + (r_min - gmax[dir]);
            if (r_rel > gmax[dir]) r_rel = gmin[dir] + (r_rel - gmax[dir]);

            double eps = step[dir]/10;
            
            int ig, ix = CI.i.x, iy = CI.i.y, iz = CI.i.z; 
            ig = ix + (iy + iz*CI.n.y)*CI.n.x; 
            if ((cell_min[dir]+eps >= r_min and cell_max[dir]-eps <= r_rel) ||
                ( (r_min > r_rel) and ((cell_min[dir] + eps >= r_min) || (cell_max[dir] - eps <= r_rel)))){ 

                if (CI.particleTypeIndex==0){ // remove particles when electrons are called to count this time separately
                    // removing particles 
                    if(cell[ig] != nullptr){
                        int it = 0;
                        if(cell[ig][it] != nullptr){
                            cellContainer *C = cell[ig][it];
                            if(C->endShift > 0) memcpy(&C->P[0], &C->P[C->P.size() - cell[ig][it]->endShift], sizeof(particle)*C->endShift);
                            C->P.resize(C->endShift);
                        };
                    };
                };
                if (CI.particleTypeIndex==-1){
                    threadHandler &cthread(Thread[CI.threadNum]); // cthread (current thread) is to run in a thread-safe way
                    cthread.rng.seed(CI.rngSeed);
                    
                    // adding particles
                    double3 r = (CI.cellMin() + 0.5*CI.step);
                    double R[3]; 
                    R[0] = r.x;
                    R[1] = r.y;
                    R[2] = r.z;
                    
                    // The position of the front of the window in 'real' coordinates
                    double z_real = _time[0]*_velocity + gmax[dir];// - gmax[ndir]*tan(_angle);
                    if (R[dir] > r_rel){
                        R[dir] = z_real - (r_rel - gmin[dir]) - (gmax[dir] - R[dir]);
                    } else {
                        R[dir] = z_real - (r_rel - R[dir]);
                    };
                                    
                    double(*density_profile)(double*, double*, int*) = (double(*)(double*, double*, int*))_density_profile;
                    double _density = density_profile(R, dataDouble, dataInt);
                    double nb_particles = _density*CI.step.x*CI.step.y*CI.step.z; // number of particles to be generated in the cell
                    
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
};


void set_statics(int64_t simbox, int thickness, double velocity, double angle, char direction){
    _simbox = (simulationBox*)simbox;
    _angle = angle;
    _velocity = velocity/cos(angle); // adjust velocity according to the angle
    
    if (direction == 'x') dir = 0;
    else if (direction == 'y') dir = 1;
    else if (direction == 'z') dir = 2;
    else {
        pipic_log.message("pi-PIC error: moving_window handler(): direction must be 'x', 'y' or 'z'.", true);
    };

    if (_simbox->dim == 1 && angle != 0){
        pipic_log.message("pi-PIC error: moving_window handler(): angle is not supported in 1D simulations.", true);
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
    _time = &_simbox->time;

    if (dir==2){
        ndir=0;
    }else{ ndir = dir+1; };

    if (thickness == -1){
        _thickness = int(n[dir]/8); // default thickness is equal to the cell size in the principal direction
        pipic_log.message("pi-PIC warning: moving_window handler(): thickness is not set, using default value of n/8.", false);
    }else{
        _thickness = thickness;
    };
};



// extension initialization
int64_t handler(int64_t ensembleData, int64_t simulation_box, double particles_per_cell, double temperature, int64_t density_profile, int thickness, double velocity, char axis, double angle){ 
    cell = (cellContainer***)ensembleData;
    ppc = particles_per_cell;
    _temperature = temperature;
    _density_profile = density_profile;
    if(!staticsHasBeenSet){
        set_statics(simulation_box, thickness, velocity, angle, axis);
        staticsHasBeenSet = true;
    }else{
        pipic_log.message("pi-PIC warning: moving_window handler(): "
        "static variables have already been set, using values for simulation_box, thickness, "
        "velocity, angle and direction set in field_handler (or defaults).", false);
    };
    Thread.resize(omp_get_max_threads());
    return (int64_t)Handler;
};


// extension initialization
int64_t field_handler(int64_t simulation_box, double timestep, int thickness, double velocity, char axis, double angle){ 
    _timeStep = timestep;
    if(!staticsHasBeenSet){
        set_statics(simulation_box, thickness, velocity, angle, axis);
        staticsHasBeenSet = true;
    }else{
        pipic_log.message("pi-PIC warning: moving_window field_handler(): "
        "static variables have already been set, using values for simulation_box, thickness," 
        "velocity, angle and direction set in handler (or defaults).", false);
    };
    if (timestep > step[dir]/lightVelocity){
        pipic_log.message("pi-PIC error: moving_window field_handler(): time step must be smaller than the cell size divided by the light velocity.", true);
    };
    return (int64_t)fieldHandler;
};

namespace py = pybind11;
PYBIND11_MODULE(_moving_window, object) {
    object.attr("name") = name;
    object.def("handler", &handler, 
               py::arg("ensemble"), 
               py::arg("simulation_box"), 
               py::arg("particles_per_cell"), 
               py::arg("temperature"),
               py::arg("density_profile"), 
               py::arg("thickness")=-1, 
               py::arg("velocity")=lightVelocity,
               py::arg("axis")='x',
               py::arg("angle")=0);

    object.def("field_handler", &field_handler, 
               py::arg("simulation_box"), 
               py::arg("timestep"),
               py::arg("thickness")=-1, 
               py::arg("velocity")=lightVelocity,
               py::arg("axis")='x',
               py::arg("angle")=0);
};
