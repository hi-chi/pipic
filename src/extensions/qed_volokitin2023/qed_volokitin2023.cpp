/*-------------------------------------------------------------------------------------------------------
This file is an implementaion of an extension 'qed_volokitin2023', which uses Fast_QED module of 
hi-chi/pyHiChi (Dec 2023) and ports it to the interface of pi-PIC.

qed_volokitin2023, Copyright 2024 Joel Magnusson
---------------------------------------------------------------------------------------------------------
qed_volokitin2023 is free software: you can redistribute it and/or modify it under the terms of the GNU 
General Public License as published by the Free Software Foundation, either version 3 of the License, or 
(at your option) any later version.

qed_volokitin2023 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General 
Public License for more details.

You should have received a copy of the GNU General Public License along with pi-PIC. If not, se
<https://www.gnu.org/licenses/>.
---------------------------------------------------------------------------------------------------------
Website: https://github.com/hi-chi/pipic
Contact: arkady.gonoskov@gu.se.
-------------------------------------------------------------------------------------------------------*/
// Description: The extension enables QED event generation using the minimal possible number of rate 
// computations per QED event [V. Volokitin et al. JCS (2023) https://doi.org/10.1016/j.jocs.2023.102170].

#include "interfaces.h"
#include "compton.h"
#include "breit_wheeler.h"

static pfc::Compton compton;
static pfc::Breit_wheeler breit_wheeler;

const string name = "qed_volokitin2023";
static int electronType;
static int positronType;
static int photonType;

// Simple class for gathering particles, particle-time and particle-type into a single container
class AvalancheContainer {
    public:
    vector<particle> particles;
    vector<double> times;
    vector<int> types;

    void clear() {
        particles.clear();
        times.clear();
        types.clear();
    }
};

struct threadHandler{
    mt19937 rng;
    std::uniform_real_distribution<double> U1;
    AvalancheContainer container; // Vectors for storing particles during processing
    threadHandler(): U1(0, 1.0) {}
    double random() {return U1(rng);} // returns a random number from [0, 1)
    double getDt(double rate, double chi, double gamma) {
        double r = -log(random());
        r *= (gamma) / (chi);
        return r / rate;
    }
    double photonGenerator(double chi) {
        double r = random();
        return compton.inv_cdf(r, chi);
    }
    double pairGenerator(double chi) {
        double r = random();
        if (r < 0.5)
            return breit_wheeler.inv_cdf(r, chi);
        else
            return 1.0 - breit_wheeler.inv_cdf(1.0 - r, chi);
    }
};

static vector<threadHandler> Thread;

void handleElectronOrPositron(threadHandler &cthread, int ip, const double3 &E, const double3 &B, double timeStep) {
    double mc = constants::electronMass * constants::lightVelocity;
    double time = cthread.container.times[ip];
    while (time < timeStep){
        double3 p = cthread.container.particles[ip].p;
        double gamma = sqrt(1 + sqr(p)/sqr(mc));
        double3 beta = p / (gamma*mc);
        double H_eff = sqr(E + cross(beta, B)) - sqr(dot(E, beta));
        if (H_eff < 0) H_eff = 0;
        H_eff = sqrt(H_eff);
        double chi = gamma * H_eff / pfc::schwingerField;
        
        // Compute rate and subtimestep (dt)
        double rate = 0.0, dt = timeStep; // temporal setting
        if (chi > 0.0) {
            rate = compton.rate(chi);
            dt = cthread.getDt(rate, chi, gamma);
        }

        if (dt + time > timeStep) {
            // Move particle to end of timeStep. No particle push being carried out in this part of code
            time = timeStep;
        } else {
            // Move particle by a subtimestep (dt). No particle push being carried out in this part of code
            time += dt;

            // determine new particle energy
            double delta = cthread.photonGenerator(chi);

            // Create new particle (photon)
            particle new_particle;
            new_particle.w = cthread.container.particles[ip].w;
            new_particle.r = cthread.container.particles[ip].r;
            new_particle.p = delta * p;

            // Add new particle to container for later processing
            cthread.container.particles.push_back(new_particle);
            cthread.container.times.push_back(time);
            cthread.container.types.push_back(photonType);
            // Change current particle momentum
            cthread.container.particles[ip].p = (1 - delta)*p;
        }
    }
    cthread.container.times[ip] = timeStep; // not needed, but for keeping things clear.
};

void handlePhoton(threadHandler &cthread, int ip, const double3 &E, const double3 &B, double timeStep) {
    double mc = constants::electronMass * constants::lightVelocity;
    double time = cthread.container.times[ip];
    while (time < timeStep) {
        double3 p = cthread.container.particles[ip].p;
        double pNorm = p.norm();
        if(pNorm == 0) return;
        double gamma = pNorm/mc; // gamma = E/mc2 = |p|/mc
        double3 k_ = p/pNorm; // normalized wave vector
        double H_eff = sqrt(sqr(E + cross(k_, B)) - sqr(dot(E, k_)));
        double chi = gamma * H_eff / pfc::schwingerField;
        
        // Compute rate and subtimestep (dt)
        double rate = 0.0, dt = timeStep;
        if (chi > 0.0) {
            rate = breit_wheeler.rate(chi);
            dt = cthread.getDt(rate, chi, gamma);
        }

        if (dt + time > timeStep) {
            // Move particle to end of timeStep. No particle push being carried out in this part of code
            time = timeStep;
        } else {
            time += dt; // Move photon by a subtimestep (dt).
            
            // determine new particle energy
            double delta = cthread.pairGenerator(chi);

            // Create new particle (electron)
            particle new_particle;
            new_particle.w = cthread.container.particles[ip].w;
            new_particle.r = cthread.container.particles[ip].r;
            new_particle.p = delta * p;

            // Add new particle to container for later processing
            cthread.container.particles.push_back(new_particle);
            cthread.container.times.push_back(time);
            cthread.container.types.push_back(electronType);

            // Create new particle (positron)
            new_particle.p = (1-delta) * p;

            // Add new particle to container for later processing
            cthread.container.particles.push_back(new_particle);
            cthread.container.times.push_back(time);
            cthread.container.types.push_back(positronType);

            // Change current particle momentum
            cthread.container.particles[ip].w = 0.0; // Marks particle for deletion
            return;
        }
    }
}

void RunAvalanche(threadHandler &cthread, const double3 E, const double3 B, double timeStep, AvalancheContainer &avalancheContainer) {
    for (int ip = 0; ip < int(avalancheContainer.particles.size()); ip++) {
        int typeId = avalancheContainer.types[ip];
        // Process a single particle a full timeStep
        if (typeId == electronType || typeId == positronType) {
            handleElectronOrPositron(cthread, ip, E, B, timeStep);
        } else if (typeId == photonType) {
            handlePhoton(cthread, ip, E, B, timeStep);
        }
    }
}

// function that is called to process particles in a given cell
void Handler(int *I, double *D, double *F, double *P, double *NP, double *dataDouble, int *dataInt) {
    // interface for manipulating with the content of a cell
    cellInterface CI(I, D, F, P, NP);
    if((CI.particleTypeIndex == electronType)||(CI.particleTypeIndex == positronType)||(CI.particleTypeIndex == photonType))
    for(int ip = 0; ip < CI.particleSubsetSize; ip++) {
        threadHandler &cthread(Thread[CI.threadNum]); // cthread (current thread) is to run in a thread-safe way
        if(ip == 0) cthread.rng.seed(CI.rngSeed);
        particle *P = CI.Particle(ip); // Grab a single particle from CellInterface

        double3 E, B; 
        CI.interpolateField(P->r, E, B); // Compute the field at particle position (E, B are passed by reference)

        cthread.container.clear(); // Clear particle container and place there the particle
        cthread.container.particles.push_back(*P);
        cthread.container.times.push_back(0.0);
        cthread.container.types.push_back(CI.particleTypeIndex);

        RunAvalanche(cthread, E, B, CI.timeStep, cthread.container); // Run a full timeStep for the particle (including for its daughter particles)
        // Update particle momentum and weight (in case of removal)
        P->p = cthread.container.particles[0].p;
        P->w = cthread.container.particles[0].w;

        // Place all created particles in the buffer
        for (int k = 1; k < int(cthread.container.particles.size()); k++) {
            // Make sure that the particle isn't marked for deletion
            if (cthread.container.particles[k].w > 0.0) {
                if(CI.particleBufferSize < CI.particleBufferCapacity){ // Check if the buffer permits adding a particle
                    *CI.newParticle(CI.particleBufferSize) = cthread.container.particles[k]; // Copy particle to the new-particle buffer
                    CI.newParticle(CI.particleBufferSize)->id = cthread.container.types[k]; // Set new-particle type id
                    CI.particleBufferSize++;
                }
            }
        }
    }
};

// extension initialization
int64_t handler(int electronType_, int positronType_, int photonType_){
    electronType = electronType_;
    positronType = positronType_;
    photonType = photonType_;
    Thread.resize(omp_get_max_threads());
    return (int64_t)Handler;
};

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/operators.h"

namespace py = pybind11;
PYBIND11_MODULE(_qed_volokitin2023, object) {
    object.attr("name") = name;
    object.def("handler", &handler, py::arg("electron_type"), py::arg("positron_type"), py::arg("photon_type"));
}
