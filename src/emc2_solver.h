/*-------------------------------------------------------------------------------------------------------
This file is part of pi-PIC.
pi-PIC, Copyright 2024 Arkady Gonoskov
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

struct mat4{
    complex<double> v[4][4];
    void print(){
        for(int i = 0; i < 4; i++){
            for(int j = 0; j < 4; j++) cout << v[i][j] << char(9);
            cout << endl;
        }
    }
    void set_rnd(){
        for(int i = 0; i < 4; i++)
        for(int j = 0; j < 4; j++)
            v[i][j] = rand_int(-1, 1);
    }
    void setToNull(){
        for(int i = 0; i < 4; i++)
        for(int j = 0; j < 4; j++)
            v[i][j] = 0;
    }
};

void setJSS_(double b, mat4& J, mat4 &S, mat4 &S_){
    double a = sqrt(b*b + 4);
    double p = b*b + 2 + b*a;
    double q = b*b + 2 - b*a;
    double o1 = sqrt(0.5*p), o2 = sqrt(0.5*q);
    J.setToNull(); S.setToNull(); S_.setToNull();
    J.v[0][0].imag(-o1);
    J.v[1][1].imag(o1);
    J.v[2][2].imag(-o2);
    J.v[3][3].imag(o2);
    
    S.v[0][0] = p/(b+a); S.v[0][1] = S.v[0][0]; S.v[0][2] = q/(b-a); S.v[0][3] = S.v[0][2];
    S.v[1][0].imag(-o1); S.v[1][1].imag(o1); S.v[1][2].imag(-o2); S.v[1][3].imag(o2);
    S.v[2][0].imag(2*o1/(b+a)); S.v[2][1].imag(-S.v[2][0].imag()); S.v[2][2].imag(2*o2/(b-a)); S.v[2][3].imag(-S.v[2][2].imag());
    S.v[3][0] = 1; S.v[3][1] = 1; S.v[3][2] = 1; S.v[3][3] = 1;

    S_.v[0][0] = 1/(2*a); S_.v[0][1].imag((a+b)/(4*a*o1)); S_.v[0][2].imag(-1/(2*a*o1)); S_.v[0][3]=0.25 - b/(4*a);
    S_.v[1][0] = S_.v[0][0]; S_.v[1][1].imag(-S_.v[0][1].imag()); S_.v[1][2].imag(-S_.v[0][2].imag()); S_.v[1][3] = S_.v[0][3];
    S_.v[2][0] = -S_.v[0][0]; S_.v[2][1].imag(-(b-a)/(4*a*o2)); S_.v[2][2].imag(1/(2*a*o2)); S_.v[2][3] = 0.25*(b/a + 1);
    S_.v[3][0] = -S_.v[0][0]; S_.v[3][1].imag(-S_.v[2][1].imag()); S_.v[3][2].imag(-S_.v[2][2].imag()); S_.v[3][3] = S_.v[2][3];
}

void multiply(mat4 &AB, mat4 &A, mat4 &B){
    for(int i = 0; i < 4; i++)
    for(int j = 0; j < 4; j++){
        AB.v[i][j] = 0;
        for(int k = 0; k < 4; k++)
            AB.v[i][j] += A.v[i][k]*B.v[k][j];
    }
};

void multiply(vector<complex<double>> &Ax, mat4 &A, vector<complex<double>> &x){
    for(int i = 0; i < 4; i++){
        Ax[i] = 0;
        for(int j = 0; j < 4; j++)
            Ax[i] += A.v[i][j]*x[j];
    }
};

template<typename vec>
void hSolver(vec &x, vec &xt, double dt){ // advances the state of harmonic oscillator xtt + x = 0 by time dt
    double b = dt;
    double cosb_1, sinb;
    if(abs(b) > 1e-5){
        cosb_1 = cos(b) - 1;
        sinb = sin(b);
    } else {
        cosb_1 = -0.5*sqr(b) + (1/24.0)*sqr(sqr(b));
        sinb = b - (1/6.0)*b*sqr(b) + (1/120.0)*b*sqr(sqr(b));
    }
    vec tmp = xt + cosb_1*xt - sinb*x;
    x = (1 + cosb_1)*x + sinb*xt;
    xt = tmp;
};

struct magnetoharmonicSolver{ // solves ptt = -p + cross(pt, B)
    // should be allocated for each thread
    vector<complex<double>> u, v, v_;
    mat4 J, S, S_;
    magnetoharmonicSolver(): u(4), v(4), v_(4)
    {}
    void solve(double3 &p, double3 &pt, double3 B, double dt){
        //introduce an orthonormal basis (b0, b1, b2) with b2 along B 
        double B_norm = B.norm();
        double3 b2(1, 0, 0);
        if(abs(dt*B_norm) > 1e-15) b2 = (1/B_norm)*B;
        double3 b0 = cross(b2, double3(1, 0, 0)); 
        if(b0.norm2() < 1e-20) b0 = cross(b2, double3(0, 1, 0)); // in case b2 is perpendcular to x
        if(b0.norm2() < 1e-20) b0 = cross(b2, double3(0, 0, 1)); //
        b0.normalize();
        double3 b1 = cross(b2, b0);
        
        double p0, p1, p2, pt0, pt1, pt2; // 2 refers to component along B
        p0 = dot(p, b0); pt0 = dot(pt, b0);
        p1 = dot(p, b1); pt1 = dot(pt, b1);
        p2 = dot(p, b2); pt2 = dot(pt, b2);

        hSolver<double>(p2, pt2, dt);    
        setJSS_(B_norm, J, S, S_);

        u[0] = pt0; u[1] = pt1; u[2] = p0; u[3] = p1;
        multiply(v_, S_, u);

        complex<double> exp0(cos(dt*J.v[0][0].imag()), sin(dt*J.v[0][0].imag()));
        complex<double> exp1(exp0.real(), -exp0.imag());
        complex<double> exp2(cos(dt*J.v[2][2].imag()), sin(dt*J.v[2][2].imag()));
        complex<double> exp3(exp2.real(), -exp2.imag());

        v[0] = exp0*v_[0];
        v[1] = exp1*v_[1];
        v[2] = exp2*v_[2];
        v[3] = exp3*v_[3];

        multiply(u, S, v);

        pt0 = u[0].real();
        pt1 = u[1].real();
        p0 = u[2].real();
        p1 = u[3].real();
        
        p = p0*b0 + p1*b1 + p2*b2;
        pt = pt0*b0 + pt1*b1 + pt2*b2;
    }
}; 

struct emc2_solver: public pic_solver //energy-conserving solver
{
    fourierSolver *field;
    vector<unsigned long long int> overStepMove; // counter of overStepMigrations
    // auxiliary variables for optimization:
    struct threadData{
        double val[8]; 
        magnetoharmonicSolver mhSolver;
    };
    vector<threadData> data; // thread-local
    double invCellVolume;
    vector<fieldSubMap64> subField;
    int en_corr_type;

    emc2_solver(simulationBox box): overStepMove(omp_get_max_threads(), 0),
    data(omp_get_max_threads()), subField(omp_get_max_threads())
    {
        field = new fourierSolver(box, FFTW_PATIENT);
        Field = field;
        Ensemble = new ensemble(Field->box);
        Ensemble->shuffle = false;
        invCellVolume = 1/(Field->box.step.x*Field->box.step.y*Field->box.step.z);
        name = "emc2";
        en_corr_type = 2;
    }
    ~emc2_solver(){
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
        //double &qdt_4mc(data[thread].val[1]);
        double &_4piq2_mVg(data[thread].val[2]);
        double &q_mc(data[thread].val[3]);
        double &mc_q(data[thread].val[4]);
        //double &Vg_8pimc2(data[thread].val[5]);
        double &mc(data[thread].val[6]);
        double &Vg_4piq(data[thread].val[7]);
        
        double p2_ = P.p.norm2()*sqr(inv_mc);
        double gamma = sqrt(1 + p2_);
        double inv_gamma = 1/gamma;
        double3 vdt_4 = 0.25*timeStep*inv_gamma*lightVelocity*inv_mc*P.p;
        moveCap(vdt_4, 0.4999999*box.step);

        double c[8]; int cil[8];
        double3 E_({0, 0, 0}), B_({0, 0, 0});
        CIC(P.r + vdt_4, c, cil, E_, B_);
        double xi = 0;
        int dim = box.dim;
        for(int i = 0; i < (1 << dim); i++) xi += sqr(c[i]);
        double sqrt_kappa = sqrt(_4piq2_mVg*P.w*inv_gamma*xi);

        //double3 E = (q_mc/sqrt_kappa)*E_;
        double3 B = (q_mc/(sqrt_kappa*gamma))*B_;

        double3 p = inv_mc*P.p;
        double3 p_b(p); // before

        double3 pt = (1/sqrt_kappa)*q_mc*E_ + cross(p, B);
        double3 pt0 = pt, p0 = p;
        
        data[thread].mhSolver.solve(p, pt, B, 0.5*sqrt_kappa*timeStep);

        double3 dp = (p - p0);
        double3 dE = mc_q*(1/xi)*sqrt_kappa*(pt - pt0 - cross(dp, B));

        P.p = mc*p;
        for(int i = 0; i < (1 << dim); i++) map.E_(cil[i]) += c[i]*dE;

        double3 dr = (-Vg_4piq/P.w)*dE;
        moveCap(dr, 0.4999999*box.step);
        P.r += dr;
    }
    void advance(double timeStep){
        Ensemble->advance_doubleLoop<emc2_solver, fourierSolver>(this, timeStep);
    }
};