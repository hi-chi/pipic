/*-------------------------------------------------------------------------------------------------------
This file is an implementation of an extension qed_gonoskov2015, which is compartible with pi-PIC.
qed_gonoskov2015, Copyright 2024 Arkady Gonoskov
---------------------------------------------------------------------------------------------------------
qed_gonoskov2015 is free software: you can redistribute it and/or modify it under the terms of the GNU 
General Public License as published by the Free Software Foundation, either version 3 of the License, or 
(at your option) any later version.

qed_gonoskov2015 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General 
Public License for more details.

You should have received a copy of the GNU General Public License along with qed_gonoskov2015. If not, se
<https://www.gnu.org/licenses/>.
---------------------------------------------------------------------------------------------------------
Website: https://github.com/hi-chi/pipic
Contact: arkady.gonoskov@gu.se.
-------------------------------------------------------------------------------------------------------*/
// Description: a simple implementation of QED event generation based on rejection sampling and 
// subcycling [A. Gonoskov et al. PRE (2015), arXiv:1412.6426]

#include "interfaces.h"

const string name = "qed_gonoskov2015";
static int electronType;
static int positronType;
static int photonType;
static double probabilityThreshold; // threshold below which events are rerified getting correspondingly increased weight factor
static double probabilitySubcycle; // maximal estimated probability for a QED event within a single time substep 

//implementation of approximations for the first and second synchrotron functions following [Fouka and Ouichaoui, arXiv:1301.6908]
const double F1 = 2.149528241534479;
const double F2 = 1.2533141373155001;
const double G1 = 1.0747641207672396;
const double G2 = 1.2533141373155001; 
float K_2_3(float x){ // relative error < 0.0054
    float x13 = pow(x, 1/3.0), sqrtx = sqrt(x);
    float H1 = -1.3746667760953621*x + 0.44040512552162292*sqrtx - 0.15527012012316799*x13;
    float H2 = -0.33550751062084*x;
    return G1*(1/sqr(x13))*exp(H1) + G2*exp(-x)*(1/sqrtx)*(1 - exp(H2));
}
float synchFunc2(float x){
    return x*K_2_3(x);
}
float synchFunc1(float x){ // relative error < 0.0026
    float x13 = pow(x, 1/3.0), sqrtx = sqrt(x);
    float H1 = -0.97947838884478688*x - 0.83333239129525072*sqrtx + 0.15541796026816246*x13;
    float H2 = -0.0469247165562628882*x - 0.70055018056462881*sqrtx + 0.0103876297841949544*x13;  
    return F1*x13*exp(H1) + F2*exp(-x)*sqrtx*(1 - exp(H2));
}

//qed probabilities (LCFA)
const double plancksConstant = 6.62607015e-27;
const double hbar = plancksConstant/(2*pi);
const double schwingerField = sqr(electronMass*lightVelocity)*lightVelocity/(-electronCharge*hbar);
const double const1 = (2*pi/137.035999084)*(sqr(hbar)/(sqr(electronCharge)*electronMass*lightVelocity));

inline double Synchrotron_I(double chi, double gamma, double d)
{
    double H_eff = (chi/gamma)*(sqr(electronMass*lightVelocity)*lightVelocity)/(electronCharge*hbar);
    double z = (2/3.0)*(1/chi)*d/(1 - d);
    if((z < 700)&&(z > 0))
        return (sqrt(3)/(2*pi))*((sqr(electronCharge)*electronCharge*H_eff)/(electronMass*sqr(lightVelocity)))*(1 - d)*(synchFunc1(z) + (3/2.0)*d*chi*z*synchFunc2(z));
    else
        return 0;
}
inline double Photon_probability(double chi, double gamma, double d)
{
    double z = (2/3.0)*(1/chi)*d/(1 - d);
    if((z < 700)&&(z > 0))
        return (sqrt(3)/(2*pi))*(chi/gamma)*((1 - d)/d)*(synchFunc1(z) + (3/2.0)*d*chi*z*synchFunc2(z));
    else
        return 0;
}
inline double Pair_probability(double chi, double gamma, double d)
{
    double z_p = (2/3.0)/(chi*(1 - d)*d);
    if((z_p < 700)&&(z_p > 0))
        return (sqrt(3)/(2*pi))*(chi/gamma)*(d - 1)*d*(synchFunc1(z_p) - (3/2.0)*chi*z_p*synchFunc2(z_p));
    else
        return 0;
}
inline double Photon_Generator(double Factor, double chi, double gamma, double dt, double r1, double r2) //returns photon energy in mc2gamma in case of generation, !doesn't change gamma
{
    double factor = Factor*dt*sqr(electronCharge)*electronMass*lightVelocity/sqr(hbar);
    if(r2 < factor*Photon_probability(chi, gamma, r1))
        return r1;
    else
        return 0;
}
inline double Pair_Generator(double Factor, double chi, double gamma, double dt, double r1, double r2) //returns photon energy in mc2gamma in case of generation.
{
    double factor = Factor*dt*sqr(electronCharge)*electronMass*lightVelocity/sqr(hbar);
    if(r2 < factor*Pair_probability(chi, gamma, r1))
        return r1;
    else
        return 0;
}
inline double Photon_MGenerator(double Factor, double chi, double gamma, double dt, double r0, double r2) //Modified event generator: returns photon energy in mc2gamma in case of generation, !doesn't change gamma
{
    double r1 = r0*r0*r0;
    double factor = Factor*dt*sqr(electronCharge)*electronMass*lightVelocity/sqr(hbar);
    if(r2 < factor*Photon_probability(chi, gamma, r1)*3*r0*r0)
        return r1;
    else
        return 0;
}

void addParticle(cellInterface &CI, particle &P){
    if(CI.particleBufferSize < CI.particleBufferCapacity){ // checking if the buffer permits adding a particle 
        *CI.newParticle(CI.particleBufferSize) = P; // copy particle to a new particle (buffer)
        CI.particleBufferSize++;
    }
}

struct cascadeParticle{double weight, gamma; int type;};

struct threadHandler{
    mt19937 rng;
    std::uniform_real_distribution<double> U1;
    vector<cascadeParticle> cParticle;
    threadHandler(): U1(0, 1.0) {}
    double random() {return U1(rng);} // returns a random number from [0, 1)
    void handleCascade(double E_E_cr, double timeStep, int seed_ip, cellInterface &CI){ // E_E_cr = H_eff/schwingerField
        double gamma0_1 = 1/cParticle[0].gamma;
        double timeSubStep = probabilitySubcycle*14*const1*pow(E_E_cr, -2/3.0);
        int numberOfSubSteps = 1 + int(timeStep/timeSubStep);
        timeSubStep = timeStep/double(numberOfSubSteps);
        for(int is = 0; is < numberOfSubSteps; is++){
            int nParticles = int(cParticle.size());
            for(int ip = 0; ip < nParticles; ip++){
                if(cParticle[ip].type == photonType){
                    // handling photon
                    if(cParticle[ip].weight > 0){
                        double delta = Pair_Generator(1.0, E_E_cr*cParticle[ip].gamma, cParticle[ip].gamma, timeSubStep, random(), random());
                        if(delta > 0){
                            cParticle.push_back({cParticle[ip].weight, cParticle[ip].gamma*delta, electronType});
                            cParticle.push_back({cParticle[ip].weight, cParticle[ip].gamma*(1 - delta), positronType});
                            cParticle[ip].weight = 0; // removing the photon
                        }
                    }
                } else {
                    // handling electron or positron
                    double delta = Photon_MGenerator(1.0, E_E_cr*cParticle[ip].gamma, cParticle[ip].gamma, timeSubStep, random(), random());
                    if(delta > 0){
                        cParticle.push_back({cParticle[ip].weight, cParticle[ip].gamma*delta, photonType});
                        cParticle[ip].gamma *= (1 - delta);
                    }
                }
            }
        }
        for(int ip = 1; ip < int(cParticle.size()); ip++)
        if(cParticle[ip].weight > 0)
        {
            addParticle(CI, *CI.Particle(seed_ip));
            CI.newParticle(CI.particleBufferSize-1)->p = (cParticle[ip].gamma*gamma0_1)*CI.Particle(seed_ip)->p;
            CI.newParticle(CI.particleBufferSize-1)->w = cParticle[ip].weight;
            CI.newParticle(CI.particleBufferSize-1)->id = cParticle[ip].type;
        }
        CI.Particle(seed_ip)->p = (cParticle[0].gamma*gamma0_1)*CI.Particle(seed_ip)->p;
        CI.Particle(seed_ip)->w = cParticle[0].weight;
    }
    void handlePhoton(int ip, cellInterface &CI){
        double3 k = CI.Particle(ip)->p;
        double p_norm = CI.Particle(ip)->p.norm();
        if(p_norm > 0) k = (1/p_norm)*k;
        double3 E, B;
        CI.interpolateField(CI.Particle(ip)->r, E, B);
        double H_eff = sqrt(sqr(E + cross(k, B)) - sqr(dot(E, k)));
        double E_E_cr = H_eff/schwingerField;
        double estimatedProbability = CI.timeStep*E_E_cr/(130*const1);
        if(estimatedProbability < probabilityThreshold){
            //handle single event
            if(random() > estimatedProbability/probabilityThreshold) return; 
            else {
                double factor = probabilityThreshold/estimatedProbability;
                double gamma = p_norm/(electronMass*lightVelocity);
                double chi = gamma*E_E_cr;
                double delta = Pair_Generator(factor, chi, gamma, CI.timeStep, random(), random());
                if(delta > 0){
                    addParticle(CI, *CI.Particle(ip));
                    CI.newParticle(CI.particleBufferSize-1)->p = delta * CI.Particle(ip)->p;
                    CI.newParticle(CI.particleBufferSize-1)->id = electronType;
                    addParticle(CI, *CI.Particle(ip));
                    CI.newParticle(CI.particleBufferSize-1)->p = (1 - delta) * CI.Particle(ip)->p;
                    CI.newParticle(CI.particleBufferSize-1)->id = positronType;
                    CI.Particle(ip)->w = 0;
                }
            }
        } else {
            //handle cascade
            cParticle.clear();
            double gamma = p_norm/(electronMass*lightVelocity);
            cParticle.push_back({CI.Particle(ip)->w, gamma, photonType});
            handleCascade(E_E_cr, CI.timeStep, ip, CI);
        }
    }
    void handleElectronOrPositron(int ip, int type, cellInterface &CI){
        double gamma = sqrt(1 + sqr(CI.Particle(ip)->p)/sqr(electronMass*lightVelocity));
        double3 v = (1/(electronMass*gamma))*CI.Particle(ip)->p;
        double3 E, B;
        CI.interpolateField(CI.Particle(ip)->r, E, B);
        double H_eff = sqr(E + (1/lightVelocity)*cross(v, B)) - (1/sqr(lightVelocity))*sqr(dot(E, v));
        if(H_eff < 0) H_eff = 0;
        H_eff = sqrt(H_eff);
        double E_E_cr = H_eff/schwingerField;
        double estimatedProbability = CI.timeStep/(14*const1*pow(E_E_cr, -2/3.0));
        if(estimatedProbability < probabilityThreshold){
            //handle single event
            if(random() > estimatedProbability/probabilityThreshold) return; 
            else {
                double factor = probabilityThreshold/estimatedProbability;
                double chi = gamma*E_E_cr;
                double delta = Photon_MGenerator(factor, chi, gamma, CI.timeStep, random(), random());
                if(delta > 0){
                    addParticle(CI, *CI.Particle(ip));
                    CI.newParticle(CI.particleBufferSize-1)->p = delta * CI.Particle(ip)->p;
                    CI.newParticle(CI.particleBufferSize-1)->id = photonType;
                    CI.Particle(ip)->p = (1 - delta) * CI.Particle(ip)->p;
                }
            }
        } else {
            //handle cascade
            cParticle.clear();
            cParticle.push_back({CI.Particle(ip)->w, gamma, type});
            handleCascade(E_E_cr, CI.timeStep, ip, CI);
        }
    }
};
static vector<threadHandler> Thread;

void Handler(int *I, double *D, double *F, double *P, double *NP, double *dataDouble, int *dataInt){
    cellInterface CI(I, D, F, P, NP);
    threadHandler &cthread(Thread[CI.threadNum]); // cthread (current thread) is to run in a thread-safe way
    cthread.rng.seed(CI.rngSeed);
    if(CI.particleTypeIndex == photonType) for(int ip = 0; ip < CI.particleSubsetSize; ip++) cthread.handlePhoton(ip, CI);
    if(CI.particleTypeIndex == electronType) for(int ip = 0; ip < CI.particleSubsetSize; ip++) cthread.handleElectronOrPositron(ip, electronType, CI);
    if(CI.particleTypeIndex == positronType) for(int ip = 0; ip < CI.particleSubsetSize; ip++) cthread.handleElectronOrPositron(ip, positronType, CI);
};

// extension initialization
int64_t handler(int electronType_, int positronType_, int photonType_, 
                double probabilityThreshold_, double probabilitySubcycle_){
    electronType = electronType_;
    positronType = positronType_;
    photonType = photonType_;
    probabilityThreshold = probabilityThreshold_;
    probabilitySubcycle = probabilitySubcycle_;
    Thread.resize(omp_get_max_threads());
    return (int64_t)Handler;
};

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/operators.h"

namespace py = pybind11;
PYBIND11_MODULE(_qed_gonoskov2015, object) {
    object.attr("name") = name;
    object.def("handler", &handler, py::arg("electron_type"), py::arg("positron_type"), py::arg("photon_type"),
                py::arg("probability_threshold") = 1e-3, py::arg("probability_subcycle") = 0.1);
    object.def("synch_func_1", &synchFunc1, py::arg("x"));
    object.def("synch_func_2", &synchFunc2, py::arg("x"));
}
