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
// Description: Implementation of second-order energy-conserving pic solver.

#include "ensemble.h"
#include "fourier_solver.h"

struct ec2_solver: public pic_solver //energy-conserving solver
{
    fourierSolver *field;
    vector<unsigned long long int> overStepMove; // counter of overStepMigrations
    // auxiliary variables for optimization:
    struct threadData{double val[8];};
    vector<threadData> data; // thread-local
    double invCellVolume;
    vector<fieldSubMap64> subField;
    int en_corr_type;
    
    ec2_solver(simulationBox box): overStepMove(omp_get_max_threads(), 0),
    data(omp_get_max_threads()), subField(omp_get_max_threads())
    {
        field = new fourierSolver(box, FFTW_PATIENT);
        Field = field;
        Ensemble = new ensemble(Field->box);
        Ensemble->shuffle = false;
        invCellVolume = 1/(Field->box.step.x*Field->box.step.y*Field->box.step.z);
        name = "ec2";
        en_corr_type = 2;
    }
    ~ec2_solver(){
        delete field;
        delete Ensemble;
    }
    void preLoop(int loopNumber){
        if(loopNumber == 0){
            if(field->divergenceCleaning) field->setRhoToZero();
            for(int i = 0; i < omp_get_max_threads(); i++) overStepMove[i] = 0;
        }
    }
    void postLoop(int loopNumber){
        if(loopNumber == 1){
            unsigned long long int totalOverStepMoves = 0;
            for(int i = 0 ; i < omp_get_max_threads(); i++) totalOverStepMoves += overStepMove[i];
            if(totalOverStepMoves != 0){
                string warning = "ec2_solver warning: restricting overstep move in " + to_string(totalOverStepMoves);
                warning += " (" + to_string(100*0.5*totalOverStepMoves/double(Ensemble->totalNumberOfParticles)) + "%) occasions (two checks ber update). ";
                warning +=  "Consider reducing time step.";
                pipic_log.message(warning);
            }
        }
    }
    inline void moveCap(double3 &dr, double3 step){
        if(unlikely(abs(dr.x) > step.x)){dr = (step.x/abs(dr.x))*dr; overStepMove[omp_get_thread_num()]++;}
        //if(dim > 1)if(unlikely(abs(dr.y) > step.y)){dr = (step.y/abs(dr.y))*dr; overStepMove[omp_get_thread_num()]++;}
        //if(dim > 2)if(unlikely(abs(dr.z) > step.z)){dr = (step.z/abs(dr.z))*dr; overStepMove[omp_get_thread_num()]++;}
    }
    void startSubLoop(int3 i3, double charge, double mass, double timeStep, int loopNumber){
        subField[omp_get_thread_num()].downloadField(i3, *field);
        double Vg = field->box.step.x*field->box.step.y*field->box.step.z;
        int thread = omp_get_thread_num();
        data[thread].val[0] = 1/(mass*lightVelocity); //inv_mc
        data[thread].val[1] = 0.25*charge*timeStep*data[thread].val[0]; // qdt_4mc
        data[thread].val[2] = 4*pi*sqr(charge)/(mass*Vg); // _4piq2_mVg
        data[thread].val[3] = charge*data[thread].val[0]; // q_mc
        data[thread].val[4] = mass*lightVelocity/charge; // mc_q
        data[thread].val[5] = Vg/(8*pi*mass*sqr(lightVelocity)); //Vg_8pimc2
        data[thread].val[6] = mass*lightVelocity; //mc
        data[thread].val[7] = Vg/(4*pi*charge); //Vg_4piq
    }
    void endSubLoop(int loopNumber){
        subField[omp_get_thread_num()].uploadField(*field);
    }
    void CIC(double3 r, double *c, int *cil, double3 &E, double3 &B)
    {
        subField[omp_get_thread_num()].CIC_c(r, c, cil);
        subField[omp_get_thread_num()].CIC(c, cil, E, B);
    }
    void processParticle(particle &P, double charge, double mass, double timeStep, int loopNumber){
        if(timeStep == 0) return;
        fieldSubMap64 &map(subField[omp_get_thread_num()]);
        simulationBox &box(field->box);
        int thread = omp_get_thread_num();

        double &inv_mc(data[thread].val[0]);
        double &qdt_4mc(data[thread].val[1]);
        double &_4piq2_mVg(data[thread].val[2]);
        double &q_mc(data[thread].val[3]);
        double &mc_q(data[thread].val[4]);
        double &Vg_8pimc2(data[thread].val[5]);
        double &mc(data[thread].val[6]);
        double &Vg_4piq(data[thread].val[7]);

        double p2_ = P.p.norm2()*sqr(inv_mc);
        double gamma = sqrt(1 + p2_);
        double inv_gamma = 1/gamma;

        double3 vdt_4 = 0.25*timeStep*inv_gamma*lightVelocity*inv_mc*P.p;
        moveCap(vdt_4, 0.4999999*box.step);

        double c[8]; int cil[8];
        double3 E({0, 0, 0}), B({0, 0, 0});
        CIC(P.r + vdt_4, c, cil, E, B);

        //Boris rotation:
        if(loopNumber == 0){
            double3 t = (qdt_4mc*inv_gamma)*B;
            double3 u1 = P.p + cross(P.p, t);
            double3 s = (2/(1 + t.norm2()))*t;
            P.p += cross(u1, s);
        }

        double3 p = inv_mc*P.p; // dimensionless momentum, like in the paper
        double3 p_b(p);
        // E-p coupling
        double xi = 0;
        int dim = box.dim;
        for(int i = 0; i < (1 << dim); i++) xi += sqr(c[i]);
        double sqrt_kappa = sqrt(_4piq2_mVg*P.w*inv_gamma*xi);
        double b = 0.5*sqrt_kappa*timeStep;
        double cosb_1, sinb;
        if(abs(b) > 1e-5){
            cosb_1 = cos(b) - 1;
            sinb = sin(b);
        } else {
            cosb_1 = -0.5*sqr(b) + (1/24.0)*sqr(sqr(b));
            sinb = b - (1/6.0)*b*sqr(b) + (1/120.0)*b*sqr(sqr(b));
        }
        double3 pt = q_mc*E;
        double3 pt_diff = cosb_1*pt - sqrt_kappa*sinb*p; // pt_diff = pt_ - pt
        p = (1 + cosb_1)*p + (1/sqrt_kappa)*sinb*pt;
        double3 dE = (1/xi)*(mc_q*pt_diff);
        
        double eta = 0;
        if(en_corr_type == 1){ // achieves machine accuracy for energy conservation by multiplying p by a number close to 1; keeps dE unchanged
            double E2_diff = 0;
            for(int i = 0; i < (1 << dim); i++)
                E2_diff += -c[i]*dE.x*(2*map.E_(cil[i]).x + c[i]*dE.x) + 
                            -c[i]*dE.y*(2*map.E_(cil[i]).y + c[i]*dE.y) + 
                            -c[i]*dE.z*(2*map.E_(cil[i]).z + c[i]*dE.z);

            if(P.w*p.norm2() > 0){
                double alpha = E2_diff*Vg_8pimc2;
                double p__2 = 1/(sqr(P.w)*p.norm2());
                double h = p__2*(alpha*(alpha + 2*P.w*gamma) + sqr(P.w)*((p_b.x - p.x)*(p_b.x + p.x) + (p_b.y - p.y)*(p_b.y + p.y) + (p_b.z - p.z)*(p_b.z + p.z)));
                if(h > 1e-7) eta = sqrt(1 + h) - 1; 
                else eta = (1/2.0)*h - (1/8.0)*h*h + (1/16.0)*h*h*h;
            }
        }
        if(en_corr_type == 2){ // (default) achieves machine accuracy for energy conservation by multiplying dE by a number close to 1; keeps p unchanged
            double C = P.w*dot(p - p_b, p_b + p)/(gamma + sqrt(1 + p.norm2()));
            double alpha = Vg_8pimc2;
            double b = alpha*2*dot(E, dE);
            double a = alpha*xi*dE.norm2();
            double sigma = 1;
            if(a != 0){
                double D = b*b - 4*a*C; if(D > 0) D = sqrt(D); else D = 0;
                double a2_ = 1/(2*a);
                double sigma1 = (-b + D)*a2_;
                double sigma2 = (-b - D)*a2_;
                if(abs(sigma1 - 1) < abs(sigma2 - 1)) sigma = sigma1; else sigma = sigma2;
            }
            dE = sigma*dE;
        }
        P.p = ((1 + eta)*mc)*p;
        for(int i = 0; i < (1 << dim); i++) map.E_(cil[i]) += c[i]*dE;
        
        double3 dr = (-Vg_4piq/P.w)*dE;
        moveCap(dr, 0.4999999*box.step);
        P.r += dr;

        //Boris rotation:
        if(loopNumber == 1){
            //double c[8]; int cil[8];
            //double3 E({0, 0, 0}), B({0, 0, 0});
            //CIC(CI, map, P.r, c, cil, E, B);
            double p2_ = P.p.norm2()*sqr(inv_mc);
            double gamma = sqrt(1 + p2_);
            double inv_gamma = 1/gamma;
            double3 t = (qdt_4mc*inv_gamma)*B;
            double3 u1 = P.p + cross(P.p, t);
            double3 s = (2/(1 + t.norm2()))*t;
            P.p += cross(u1, s);
        }
    }
    void advance(double timeStep){
        Ensemble->advance_doubleLoop<ec2_solver, fourierSolver>(this, timeStep);
    }
};
