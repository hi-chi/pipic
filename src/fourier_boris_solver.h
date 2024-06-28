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

struct fourier_boris_solver: public pic_solver
{
    int3 n;
    fourierSolver *field;
    bool overStepMigration;
    double3 *dtJ; // array for storing the cumulative current multiplied by timeStep
    // auxiliary variables for optimization:
    double invCellVolume;
    vector<fieldSubMap64> subField; // local field state, thread-local
    vector<double> inv_mc_, qdt_2_; // thread-local
    fourier_boris_solver(simulationBox box): n(box.n),
    subField(omp_get_max_threads()), inv_mc_(omp_get_max_threads()), qdt_2_(omp_get_max_threads())
    {
        field = new fourierSolver(box, FFTW_PATIENT);
        Field = field;
        Ensemble = new ensemble(Field->box);
        Ensemble->shuffle = false;
        dtJ = new double3[Field->box.ng];
        memset((void*)&dtJ[0], 0, sizeof(double3)*intg(n.x)*n.y*n.z); // set cuttent to zero everywhere
        invCellVolume = 1/(Field->box.step.x*Field->box.step.y*Field->box.step.z);
        field->enableDivergenceCleaning();
        name = "fourier_boris";
    }
    ~fourier_boris_solver(){
        delete Ensemble;
        delete field;
        delete []dtJ;
    }
    void preLoop()
    {
        #pragma omp parallel for collapse(3)
        for(int iz = 0; iz < n.z; iz++)
        for(int iy = 0; iy < n.y; iy++)
        for(int ix = 0; ix < n.x; ix++)
        {
            field->Ex(field->box.ig({ix, iy, iz})) -= 2*pi*dtJ[field->box.ig({ix, iy, iz})].x;
            field->Ey(field->box.ig({ix, iy, iz})) -= 2*pi*dtJ[field->box.ig({ix, iy, iz})].y;
            field->Ez(field->box.ig({ix, iy, iz})) -= 2*pi*dtJ[field->box.ig({ix, iy, iz})].z;
        }
        memset((void*)&dtJ[0], 0, sizeof(double3)*intg(n.x)*n.y*n.z); // set cuttent to zero everywhere
        if(field->divergenceCleaning) field->setRhoToZero();
        overStepMigration = false;
    }
    void postLoop() // here we deposits the second half of the current, i.e. E -= 2*pi*timeStep*J
    {
        if(overStepMigration){
            pipic_log.message("fourier_boris_solver warning: restricting move by more than one space step. Consider reducing the time step.");
            overStepMigration = false;
        }
        #pragma omp parallel for collapse(3)
        for(int iz = 0; iz < n.z; iz++)
        for(int iy = 0; iy < n.y; iy++)
        for(int ix = 0; ix < n.x; ix++)
        {
            field->Ex(field->box.ig({ix, iy, iz})) -= 2*pi*dtJ[field->box.ig({ix, iy, iz})].x;
            field->Ey(field->box.ig({ix, iy, iz})) -= 2*pi*dtJ[field->box.ig({ix, iy, iz})].y;
            field->Ez(field->box.ig({ix, iy, iz})) -= 2*pi*dtJ[field->box.ig({ix, iy, iz})].z;
        }
    }
    void startSubLoop(int3 i3, double charge, double mass, double timeStep){
        inv_mc_[omp_get_thread_num()] = 1/(mass*lightVelocity);
        qdt_2_[omp_get_thread_num()] = 0.5*charge*timeStep;
        subField[omp_get_thread_num()].downloadField(i3, *field);
    }
    void endSubLoop(){
        subField[omp_get_thread_num()].uploadField(*field);
    }
    void CIC(double3 r, double *c, int *cil, double3 &E, double3 &B){
        subField[omp_get_thread_num()].CIC_c(r, c, cil);
        subField[omp_get_thread_num()].CIC(c, cil, E, B);
    }
    void processParticle(particle &P, double charge, double mass, double timeStep){
        fieldSubMap64 &map(subField[omp_get_thread_num()]);
        double &inv_mc(inv_mc_[omp_get_thread_num()]);
        double &qdt_2(qdt_2_[omp_get_thread_num()]);
        simulationBox &box(field->box);

        // computation of coupling weights
        double c[8]; int cil[8];
        double3 E({0, 0, 0}), B({0, 0, 0});
        CIC(P.r, c, cil, E, B);

        // Boris push algorithm
        double3 u0 = P.p + qdt_2*E;
        double gamma = sqrt(1 + u0.norm2()*sqr(inv_mc));
        double3 t = qdt_2*(inv_mc/gamma)*B;
        double3 u1 = u0 + cross(u0, t);
        double3 s = (2/(1 + t.norm2()))*t;
        P.p = u0 + cross(u1, s) + qdt_2*E;
        gamma = sqrt(1 + P.p.norm2()*sqr(inv_mc));
        double3 dtv = (timeStep/(mass*gamma))*P.p;

        if(unlikely(abs(dtv.x) > 0.99999*box.step.x)){ dtv = 0.99999*(box.step.x/abs(dtv.x))*dtv; overStepMigration = true;}
        if(box.dim > 1)if(unlikely(abs(dtv.y) > 0.99999*box.step.y)){ dtv = 0.99999*(box.step.y/abs(dtv.y))*dtv; overStepMigration = true;}
        if(box.dim > 2)if(unlikely(abs(dtv.z) > 0.99999*box.step.z)){ dtv = 0.99999*(box.step.z/abs(dtv.z))*dtv; overStepMigration = true;}

        P.r += dtv;

        //CIC weighting, deposition of density
        if(field->divergenceCleaning){
            intg cig[8];
            double3 r = P.r;
            r += (-0.5)*dtv; // half-step backward shift
            if(r.x < box.min.x) r.x += box.max.x - box.min.x; else if(r.x > box.max.x) r.x -= box.max.x - box.min.x;
            if(r.y < box.min.y) r.y += box.max.y - box.min.y; else if(r.y > box.max.y) r.y -= box.max.y - box.min.y;
            if(r.z < box.min.z) r.z += box.max.z - box.min.z; else if(r.z > box.max.z) r.z -= box.max.z - box.min.z;

            int cix = floor((r.x - box.min.x)*box.invStep.x); c[1] = (r.x - box.min.x)*box.invStep.x - cix; c[0] = 1 - c[1];
            int cix_ = cix + 1; if(cix_ == box.n.x) cix_ = 0;
            int ciy_ = 0, ciz_ = 0;
            int ciy = map.i3.y, ciz = map.i3.z;
            if(box.dim > 1){
                ciy = floor((r.y - box.min.y)*box.invStep.y); double wy = (r.y - box.min.y)*box.invStep.y - ciy;
                for(int i = 0; i < 2; i++) c[i+2] = c[i]*wy;
                for(int i = 0; i < 2; i++) c[i] *= 1 - wy;
                ciy_ = ciy + 1; if(ciy_ > box.n.y - 1) ciy_ -= box.n.y;
            }
            if(box.dim > 2){
                int ciz = floor((r.z - box.min.z)*box.invStep.z); double wz = (r.z - box.min.z)*box.invStep.z - ciz;
                for(int i = 0; i < 4; i++) c[i+4] = c[i]*wz;
                for(int i = 0; i < 4; i++) c[i] *= 1 - wz;
                ciz_ = ciz + 1; if(ciz_ > box.n.z - 1) ciz_ -= box.n.z;
            }
            cig[0] = field->box.ig({cix , ciy, ciz});
            cig[1] = field->box.ig({cix_, ciy, ciz});
            if(box.dim > 1){
                cig[2] = field->box.ig({cix , ciy_, ciz});
                cig[3] = field->box.ig({cix_, ciy_, ciz});
            }
            if(box.dim > 2){
                cig[4] = field->box.ig({cix , ciy , ciz_});
                cig[5] = field->box.ig({cix_, ciy , ciz_});
                cig[6] = field->box.ig({cix , ciy_, ciz_});
                cig[7] = field->box.ig({cix_, ciy_, ciz_});
            }
            for(int i = 0; i < (1 << box.dim); i++) field->Rho(cig[i]) += c[i]*P.w*charge*invCellVolume;
        }
        { //CIC weighting, deposition of current
            intg cig[8];
            double3 r = P.r;
            if(r.x < box.min.x) r.x += box.max.x - box.min.x; else if(r.x > box.max.x) r.x -= box.max.x - box.min.x;
            if(r.y < box.min.y) r.y += box.max.y - box.min.y; else if(r.y > box.max.y) r.y -= box.max.y - box.min.y;
            if(r.z < box.min.z) r.z += box.max.z - box.min.z; else if(r.z > box.max.z) r.z -= box.max.z - box.min.z;

            int cix = floor((r.x - box.min.x)*box.invStep.x); c[1] = (r.x - box.min.x)*box.invStep.x - cix; c[0] = 1 - c[1];
            int cix_ = cix + 1; if(cix_ == box.n.x) cix_ = 0;
            int ciy_ = 0, ciz_ = 0;
            int ciy = map.i3.y, ciz = map.i3.z;
            if(box.dim > 1){
                ciy = floor((r.y - box.min.y)*box.invStep.y); double wy = (r.y - box.min.y)*box.invStep.y - ciy;
                for(int i = 0; i < 2; i++) c[i+2] = c[i]*wy;
                for(int i = 0; i < 2; i++) c[i] *= 1 - wy;
                ciy_ = ciy + 1; if(ciy_ > box.n.y - 1) ciy_ -= box.n.y;
            }
            if(box.dim > 2){
                int ciz = floor((r.z - box.min.z)*box.invStep.z); double wz = (r.z - box.min.z)*box.invStep.z - ciz;
                for(int i = 0; i < 4; i++) c[i+4] = c[i]*wz;
                for(int i = 0; i < 4; i++) c[i] *= 1 - wz;
                ciz_ = ciz + 1; if(ciz_ > box.n.z - 1) ciz_ -= box.n.z;
            }
            cig[0] = field->box.ig({cix , ciy, ciz});
            cig[1] = field->box.ig({cix_, ciy, ciz});
            if(box.dim > 1){
                cig[2] = field->box.ig({cix , ciy_, ciz});
                cig[3] = field->box.ig({cix_, ciy_, ciz});
            }
            if(box.dim > 2){
                cig[4] = field->box.ig({cix , ciy , ciz_});
                cig[5] = field->box.ig({cix_, ciy , ciz_});
                cig[6] = field->box.ig({cix , ciy_, ciz_});
                cig[7] = field->box.ig({cix_, ciy_, ciz_});
            }
            for(int i = 0; i < (1 << box.dim); i++) dtJ[cig[i]] += c[i]*P.w*charge*invCellVolume*dtv;
        }
    }
    void advance(double timeStep){
        Ensemble->advance_singleLoop<fourier_boris_solver, fourierSolver>(this, timeStep);
    }
};
