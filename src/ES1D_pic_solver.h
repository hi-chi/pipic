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
// Electrostatic 1D PIC solver implementation.

#include "ensemble.h"
#include "ES1D_field_solver.h"

struct ES1DPicSolver: public pic_solver
{
    double timeStep;
    ES1DFieldSolver *field;

    ES1DPicSolver(simulationBox box){
        field = new ES1DFieldSolver(box);
        Field = field;
        Ensemble = new ensemble(Field->box);
        name = "electrostatic_1d";
    }
    ~ES1DPicSolver(){
        delete Ensemble;
        delete field;
    }

    void advance(double _timeStep){
        // Run one PIC step using the ensemble single-loop method.
        timeStep = _timeStep;
        Ensemble->advance_singleLoop<ES1DPicSolver, ES1DFieldSolver>(this, timeStep);
    }

    void halfstep(particle &P, double charge, double mass, double timeStep){
        double3 E;
        double3 B;
        field->getEB(P.r, E, B); // get electric field at the particle position
        P.p.x += timeStep*charge*E.x/2; // pull back/move forward 1/2 timestep to initalize/remove leapfrog scheme
    }

    void preStep(double timeStep){
        // Pull momenta back by dt/2 before the first full push.
        for (int it = 0; it < int(Ensemble->type.size()); it++){
            for(ensemble::nonOmpIterator iP = Ensemble->begin(it); iP < Ensemble->end(); iP++){
                particle *P = &*iP;
                halfstep(*P, Ensemble->type[it].charge, Ensemble->type[it].mass, -timeStep);
            }
        }
    }

    void preLoop()
    {c
        // Advance field with current from previous iteration, then clear current.
        field->advance(timeStep);
        for(size_t ix = 0; ix < field->Jx.size(); ix++){
            field->Jx[ix] = 0;
        }
    }

    void processParticle(particle &P, double charge, double mass, double timeStep){
        simulationBox &box(field->box);
        double3 E;
        double3 B;
        field->getEB(P.r, E, B); // get electric field at the particle position
        P.p.x += timeStep*charge*E.x; // advance particle momentum
        P.r.x += timeStep*P.p.x/(mass); // advance particle position
        if (P.r.x < box.min.x) P.r.x += (box.max.x - box.min.x); // periodic boundary conditions
        else if (P.r.x > box.max.x) P.r.x += (-box.max.x + box.min.x); // periodic boundary conditions
        // Momentum is stored per physical particle, not per macroparticle.

        // Deposit Jx to adjacent nodes using linear (CIC) weights.
        int indx = int((P.r.x - box.min.x)/box.step.x);
        double w_left = field->CIC(P.r.x, indx, box.step.x, box.min.x);
        double w_right = 1 - w_left;
        
        field->Jx[indx] += w_left*P.w*charge*(P.p.x/mass)/box.step.x; 
        field->Jx[(indx + 1) % box.n.x] += w_right*P.w*charge*(P.p.x/mass)/box.step.x;
    }

    void postStep(double timeStep){
        // Push momenta forward by dt/2 to leave leapfrog-centered representation.
        for (int it = 0; it < int(Ensemble->type.size()); it++){
            for(ensemble::nonOmpIterator iP = Ensemble->begin(it); iP < Ensemble->end(); iP++){
                particle *P = &*iP;
                halfstep(*P, Ensemble->type[it].charge, Ensemble->type[it].mass, timeStep);
            }
        }
    }

    // Hooks unused by this solver variant.
    void postLoop(){}
    void startSubLoop(int3 i3, double charge, double mass, double timeStep){}
    void endSubLoop(){}

};






