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
Contact: frida.brogren@gu.se.
-------------------------------------------------------------------------------------------------------*/
// Electrostatic 1D implicit PIC solver implementation.

#include "ensemble.h"
#include "ESI1D_field_solver.h"
#include "interfaces.h"
#include <iomanip>

struct ESI1DPicSolver: public pic_solver
{
    double timeStep;
    ESI1DFieldSolver *field;
    ESI1DFieldSolver *field_intermediate; // auxiliary pointer to store the intermediate field for the implicit scheme
    ESI1DFieldSolver *field_last; // auxiliary pointer to store the field at the previous iteration for the implicit scheme
    ensemble *Ensemble_avg_intermediate; // auxiliary ensemble copy for the intermediate implicit state
    ensemble *Ensemble_avg_last; // auxiliary ensemble copy for the previous implicit state
    double3 Ep;
    double3 Ep_intermediate;
    double3 B;

    ESI1DPicSolver(simulationBox box){
        field = new ESI1DFieldSolver(box);
        field_intermediate = new ESI1DFieldSolver(box);
        field_last = new ESI1DFieldSolver(box);
        Field = field;
        Ensemble = new ensemble(Field->box);
        Ensemble_avg_intermediate = new ensemble(Field->box);
        Ensemble_avg_last = new ensemble(Field->box);
        name = "electrostatic_1d_implicit";
    }
    ~ESI1DPicSolver(){
        delete Ensemble;
        delete Ensemble_avg_intermediate;
        delete Ensemble_avg_last;
        delete field;
        delete field_intermediate;
        delete field_last;
    }

    void advance(double _timeStep){
        // Run one PIC step using the ensemble single-loop method.
        timeStep = _timeStep;
        Ensemble->advance_singleLoop<ESI1DPicSolver, ESI1DFieldSolver>(this, timeStep);
    }

    void updateImplicitParticle(particle &P_orig, particle &P_intermediate, particle &P_last, double charge, double mass, double dx, double &epsMomentum){
        P_intermediate.r = P_orig.r + timeStep * P_last.p / mass / 2;
        if(P_intermediate.r.x < field->box.min.x) P_intermediate.r.x += (field->box.max.x - field->box.min.x);
        if(P_intermediate.r.x >= field->box.max.x) P_intermediate.r.x -= (field->box.max.x - field->box.min.x);

        field->getEB(P_intermediate.r, Ep, B);
        field_last->getEB(P_last.r, Ep_intermediate, B);
        P_intermediate.p = P_orig.p + (charge / 4) * (Ep_intermediate + Ep) * timeStep;

        int ix = int((P_intermediate.r.x - field->box.min.x) / field->box.step.x);
        double w_left = field->CIC(P_intermediate.r.x, ix, field->box.step.x, field->box.min.x);
        double w_right = 1 - w_left;
        field_intermediate->Jx[ix] += charge * P_intermediate.w * P_intermediate.p.x / mass * w_left / dx;
        field_intermediate->Jx[(ix + 1) % field_intermediate->box.n.x] += charge * P_intermediate.w * P_intermediate.p.x / mass * w_right / dx;

        epsMomentum += sqrt(sqr(P_intermediate.p.x - P_last.p.x));
    }

    void postLoop()
    {
        double dx = field->box.step.x; // grid spacing

        // copy the current state of the ensemble and field to the auxiliary variables for the implicit iteration
        Ensemble_avg_intermediate->transferParticleState(*Ensemble, true);
        Ensemble_avg_last->transferParticleState(*Ensemble, true);

        for (size_t it = 0; it < field->Ex.size(); it++){
            field_intermediate->Ex[it] = field -> Ex[it]; // initialize the intermediate field with the old field
            field_last -> Ex[it] = field -> Ex[it]; // initialize the last field with the old field
        }

        unsigned long long int particleCount = Ensemble->totalNumberOfParticles;
        double particleCountDenominator = particleCount > 0 ? double(particleCount) : 1.0;
        int nb_types = Ensemble->type.size();
        double eps_field; // convergence criterion for the field
        double eps_momentum; // convergence criterion for the momentum
        
        // implicit iteration
        size_t iterations = 0; // number of iterations for the implicit scheme 
        double eps=1;
        double tolerance = 1e-12; // convergence tolerance for the implicit scheme
        while(eps > tolerance and iterations < 100){ // maximum number of iterations to avoid infinite loop in case of non-convergence
            eps_field = 0;
            eps_momentum = 0;
            for (size_t it=0; it < field->Ex.size(); it++){
                field_intermediate->Jx[it] = 0; // reset current for the next iteration
            }

            // change pointer of this state to last state to store the data of the previous iteration
            cellContainer ***tmp_cell = Ensemble_avg_last->cell;
            Ensemble_avg_last->cell = Ensemble_avg_intermediate->cell;
            Ensemble_avg_intermediate->cell = tmp_cell;
            (field_last -> Ex).swap(field_intermediate->Ex); 

            // iterate through cells and update particle momentum, position and current with the new field and calculate the convergence criterion for the momentum
            for(intg ig = 0; ig < Ensemble->box.ng; ig++){
                for(int typeIndex = 0; typeIndex < nb_types; typeIndex++){
                    for(size_t ip = 0; ip < Ensemble->cell[ig][typeIndex]->P.size(); ip++){
                        updateImplicitParticle(Ensemble->cell[ig][typeIndex]->P[ip], 
                                               Ensemble_avg_intermediate->cell[ig][typeIndex]->P[ip], 
                                               Ensemble_avg_last->cell[ig][typeIndex]->P[ip], 
                                               Ensemble->type[typeIndex].charge, 
                                               Ensemble->type[typeIndex].mass, 
                                               dx, 
                                               eps_momentum);
                    }  
                }
            }

            // update electric field with the new current and calculate the convergence criterion for the field
            for (size_t ix = 0; ix < field->Ex.size(); ix++){
                field_intermediate->Ex[ix] = field -> Ex[ix] - timeStep*4*M_PI*field_intermediate->Jx[ix]; 
                eps_field += sqrt(sqr(field_intermediate->Ex[ix] - field_last->Ex[ix]));
            }
            eps = max(eps_field/field->Ex.size(), eps_momentum/particleCountDenominator); // overall convergence criterion is the maximum of field and momentum convergence criteria
            iterations++;
        }
        if (iterations == 100) cout << "Warning: maximum number of iterations reached without convergence, eps = " << eps << "." << endl;

        // Convert the time-centered implicit state back to the stored particle state.
        for(intg ig = 0; ig < Ensemble->box.ng; ig++){
            for(int typeIndex = 0; typeIndex < nb_types; typeIndex++){
                size_t particleCount = Ensemble->cell[ig][typeIndex]->P.size();
                double mass = Ensemble->type[typeIndex].mass;
                for(size_t ip = 0; ip < particleCount; ip++){
                    particle &targetParticle = Ensemble->cell[ig][typeIndex]->P[ip];
                    particle &sourceParticle = Ensemble_avg_intermediate->cell[ig][typeIndex]->P[ip];
                    targetParticle.r += sourceParticle.p * timeStep / mass;
                    targetParticle.p = 2 * sourceParticle.p - targetParticle.p;
                }
            }
        }

        // After convergence particleLoop can be used to redistribute particles migrated particles into respective cells. The method particleLoop does this in accordance with the particle storage strategy.
        Ensemble->particleLoop(0, "all", 0, 0, Ensemble->advanceWithOmp);

        for (size_t it = 0; it < field->Ex.size(); it++){
            field->Ex[it] = field_intermediate->Ex[it]; // update the field with the converged field from the implicit iteration
            field->Jx[it] = 2*field_intermediate->Jx[it] - field->Jx[it]; // update the current with the time-averaged current from the implicit iteration for the next iteration. Not really needed, but can be useful for diagnostic purposes.
        }
    }



    // Hooks unused by this solver variant.
    void preStep(double timeStep){}
    void preLoop(){}
    void processParticle(particle &P, double charge, double mass, double timeStep){}
    void postStep(double timeStep){}
    void startSubLoop(int3 i3, double charge, double mass, double timeStep){}
    void endSubLoop(){}

};






