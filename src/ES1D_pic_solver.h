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
Contact: arkady.gonoskov@gu.se.
-------------------------------------------------------------------------------------------------------*/
// Description: Implementation of Boris pusher with Fourier solver.

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
    
    void preLoop()
    {
        field->advance(timeStep);
        cout << "Jx before reset " << field->Jx[0] << endl;
        for(size_t ix = 0; ix < field->Jx.size(); ix++){
            field->Jx[ix] = 0;
        }
        cout << "Jx after reset " << field->Jx[0] << endl;
    }

    void postLoop(){}
    void startSubLoop(int3 i3, double charge, double mass, double timeStep){}
    void endSubLoop(){}
    
    void processParticle(particle &P, double charge, double mass, double timeStep){
        simulationBox &box(field->box);
        double3 E;
        double3 B;
        field->getEB(P.r, E, B); // get electric field at the particle position
        P.p.x += timeStep*charge*E.x; // advance particle momentum
        P.r.x += timeStep*P.p.x/(mass); // advance particle position

        // deposit current
        int indx = int(P.r.x*box.step.x + 0.5);
        field->Jx[indx] += P.w*electronCharge*(P.p.x/electronMass)/box.step.x; 
    }

    void advance(double _timeStep){
        timeStep = _timeStep;
        Ensemble->advance_singleLoop<ES1DPicSolver, ES1DFieldSolver>(this, timeStep);
    }
};






