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
// Description of the file: File introduces structures for simulation region and spectral field solver.

#ifndef FOURIER_SOLVER_H
#define FOURIER_SOLVER_H

#include "interfaces.h"

struct fourierSolver: public field_solver // spectral solver for electromagnetic field evolution with zero current
{
    fftw_complex *fdata[3];
    fftw_plan fX, bX, fY, bY, fZ, bZ; // plans for FFT, the first letter standf for forward (f) or backward (b), whereas the second letter indicates the dimension.
    double invNg; // inverse total number of cells

    bool sin2Filter;
    bool divergenceCleaning;
    fftw_complex *rho;
    complex<double> *rhoBackground;
    bool initBackground;

    fourierSolver(simulationBox box, int fftw_flags): field_solver(box),
    invNg(1/double(box.ng)), sin2Filter(false), divergenceCleaning(false), rho(nullptr)
    {
        type = "fourier";
        for(int i = 0; i < 3; i++) fdata[i] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*box.ng);
        for(int i = 0; i < 3; i++) memset(fdata[i], 0, sizeof(fftw_complex)*box.ng);

        fX = fftw_plan_many_dft(1, &box.n.x, 1, fdata[0], &box.n.x, 1, 1, fdata[0], &box.n.x, 1, 1, FFTW_FORWARD, fftw_flags);
        bX = fftw_plan_many_dft(1, &box.n.x, 1, fdata[0], &box.n.x, 1, 1, fdata[0], &box.n.x, 1, 1, FFTW_BACKWARD, fftw_flags);
        fY = fftw_plan_many_dft(1, &box.n.y, 1, fdata[0], &box.n.y, box.n.x, 1, fdata[0], &box.n.y, box.n.x, 1, FFTW_FORWARD, fftw_flags);
        bY = fftw_plan_many_dft(1, &box.n.y, 1, fdata[0], &box.n.y, box.n.x, 1, fdata[0], &box.n.y, box.n.x, 1, FFTW_BACKWARD, fftw_flags);
        fZ = fftw_plan_many_dft(1, &box.n.z, 1, fdata[0], &box.n.z, box.n.x*box.n.y, 1, fdata[0], &box.n.z, box.n.x*box.n.y, 1, FFTW_FORWARD, fftw_flags);
        bZ = fftw_plan_many_dft(1, &box.n.z, 1, fdata[0], &box.n.z, box.n.x*box.n.y, 1, fdata[0], &box.n.z, box.n.x*box.n.y, 1, FFTW_BACKWARD, fftw_flags);
    }
    void fieldLoopBody(int3 i, void(*handler_)(int*, double*, double*, double*, double*, int*), double* dataDouble_, int* dataInt_){
        intg ig = box.ig({i.x, i.y, i.z});
        double3 r3 = box.min;
        r3.x += i.x*box.step.x;
        r3.y += i.y*box.step.y;
        r3.z += i.z*box.step.z;

        double r[3] = {r3.x, r3.y, r3.z};
        int ind[4] = {i.x, i.y, i.z, 6};
        double E[3], B[3];
        E[0] = Ex(ig); B[0] = Bx(ig);
        E[1] = Ey(ig); B[1] = By(ig);
        E[2] = Ez(ig); B[2] = Bz(ig);
        handler_(ind, r, E, B, dataDouble_, dataInt_);
        Ez(ig) = E[2]; Bz(ig) = B[2];
        Ey(ig) = E[1]; By(ig) = B[1];
        Ex(ig) = E[0]; Bx(ig) = B[0];
    }
    void fieldLoop(int64_t handler, int64_t dataDouble = 0, int64_t dataInt = 0, bool useOmp = false){
        void(*handler_)(int*, double*, double*, double*, double*, int*) = (void(*)(int*, double*, double*, double*, double*, int*))handler;
        double* dataDouble_ = nullptr; if(dataDouble != 0) dataDouble_ = (double*)dataDouble;
        int* dataInt_ = nullptr; if(dataInt != 0) dataInt_ = (int*)dataInt;
        if(useOmp){
            #pragma omp parallel for collapse(3)
            for(int iz = 0; iz < box.n.z; iz++)
            for(int iy = 0; iy < box.n.y; iy++)
            for(int ix = 0; ix < box.n.x; ix++)
                fieldLoopBody({ix, iy, iz}, handler_, dataDouble_, dataInt_);
        } else {
            for(int iz = 0; iz < box.n.z; iz++)
            for(int iy = 0; iy < box.n.y; iy++)
            for(int ix = 0; ix < box.n.x; ix++)
                fieldLoopBody({ix, iy, iz}, handler_, dataDouble_, dataInt_);
        }
    }
    void getEB(double3 r, double3 &E, double3 &B){
        if(unlikely((r.x < box.min.x)||(r.x >= box.max.x)||(r.y < box.min.y)||(r.y >= box.max.y)||(r.z < box.min.z)||(r.z >= box.max.z))){
            pipic_log.message("pi-PIC::fourierSolver::getEB_CIC(r, E, B): r is outside computational region.");
            return;
        }
        int ix = floor((r.x - box.min.x)*box.invStep.x);
        int iy = 0, iz = 0;
        double c[8];
        c[1] = (r.x - box.min.x)*box.invStep.x - ix; c[0] = 1 - c[1];
        if(box.dim > 1){
            iy = floor((r.y - box.min.y)*box.invStep.y);
            double wy = (r.y - box.min.y)*box.invStep.y - iy; double wy_ = 1 - wy;
            for(int i = 0; i < 2; i++) c[i+2] = c[i]*wy;
            for(int i = 0; i < 2; i++) c[i] *= wy_;
        }
        if(box.dim > 2){
            iz = floor((r.z - box.min.z)*box.invStep.z);
            double wz = (r.z - box.min.z)*box.invStep.z - iz; double wz_ = 1 - wz;
            for(int i = 0; i < 4; i++) c[i+4] = c[i]*wz;
            for(int i = 0; i < 4; i++) c[i] *= wz_;
        }
        intg cig[8];
        cig[0] = box.ig({ix    , iy, iz});
        cig[1] = box.ig({ix + 1, iy, iz});
        if(box.dim > 1){
            cig[2] = box.ig({ix    , iy + 1, iz});
            cig[3] = box.ig({ix + 1, iy + 1, iz});
            if(box.dim > 2){
                cig[4] = box.ig({ix    , iy    , iz + 1});
                cig[5] = box.ig({ix + 1, iy    , iz + 1});
                cig[6] = box.ig({ix    , iy + 1, iz + 1});
                cig[7] = box.ig({ix + 1, iy + 1, iz + 1});
            }
        }
        E = {0, 0, 0}; B = {0, 0, 0};
        for(int i = 0; i < (1 << box.dim); i++){
            E.x += c[i]*Ex(cig[i]);
            B.x += c[i]*Bx(cig[i]);
            E.y += c[i]*Ey(cig[i]);
            B.y += c[i]*By(cig[i]);
            E.z += c[i]*Ez(cig[i]);
            B.z += c[i]*Bz(cig[i]);
        }
    }
    void customFieldLoop(int numberOfIterations, int64_t it2coord, int64_t field2data, int64_t dataDouble = 0, int64_t dataInt = 0){
        customFieldLoop_via_getEB<fourierSolver>(this, numberOfIterations, it2coord, field2data, dataDouble, dataInt);
    }
    void enableDivergenceCleaning(){
        pipic_log.message("pi-PIC::fourierSolver: enable divergence cleaning.");
        if(!divergenceCleaning){
            divergenceCleaning = true;
            initBackground = false; // indicates that the backGround has to be set; this happens at the first iteration
        }
        if(rho == nullptr) {
            rho = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*box.ng);
            memset(rho, 0, sizeof(fftw_complex)*box.ng);
            rhoBackground = new complex<double>[box.ng];
        }
    }
    void setRhoToZero(){ // sets rho to zero, applicable only in case of divergence cleaning enabled
        if(rho == nullptr) {pipic_log.message("pipic::fourierSolver error: calling nullRho() without enabling divergence cleaning.", true); exit(0);}
        else memset(&rho[0], 0, sizeof(complex<double>)*intg(box.n.x)*box.n.y*box.n.z);
    }
    ~fourierSolver()
    {
        fftw_destroy_plan(fX);
        fftw_destroy_plan(bX);
        fftw_destroy_plan(fY);
        fftw_destroy_plan(bY);
        fftw_destroy_plan(fZ);
        fftw_destroy_plan(bZ);
        for(int i = 0; i < 3; i++) fftw_free(fdata[i]);
        if(rho != nullptr){
             fftw_free(rho);
             delete []rhoBackground;
        }
    }
    double& Ex(intg ig) {return *((double*)fdata[0] + 0 + 2*ig);}
    double& Bx(intg ig) {return *((double*)fdata[0] + 1 + 2*ig);}
    double& Ey(intg ig) {return *((double*)fdata[1] + 0 + 2*ig);}
    double& By(intg ig) {return *((double*)fdata[1] + 1 + 2*ig);}
    double& Ez(intg ig) {return *((double*)fdata[2] + 0 + 2*ig);}
    double& Bz(intg ig) {return *((double*)fdata[2] + 1 + 2*ig);}
    double& Rho(intg ig) {return *((double*)rho + 0 + 2*ig);} // for divergence cleaning only
    void advance(double timeStep)
    {
        #pragma omp parallel for collapse(3)
        for(int id = 0; id < 3; id++)
        for(int iy = 0; iy < box.n.y; iy++)
        for(int ix = 0; ix < box.n.x; ix++)
            fftw_execute_dft(fZ, fdata[id] + box.ig({ix, iy, 0}), fdata[id] +  box.ig({ix, iy, 0}));

        #pragma omp parallel for collapse(3)
        for(int id = 0; id < 3; id++)
        for(int iz = 0; iz < box.n.z; iz++)
        for(int ix = 0; ix < box.n.x; ix++)
            fftw_execute_dft(fY, fdata[id] + box.ig({ix, 0, iz}), fdata[id] +  box.ig({ix, 0, iz}));

        #pragma omp parallel for collapse(3)
        for(int id = 0; id < 3; id++)
        for(int iz = 0; iz < box.n.z; iz++)
        for(int iy = 0; iy < box.n.y; iy++)
            fftw_execute_dft(fX, fdata[id] + box.ig({0, iy, iz}), fdata[id] +  box.ig({0, iy, iz}));

        if(divergenceCleaning){
            #pragma omp parallel for collapse(2)
            for(int iy = 0; iy < box.n.y; iy++)
            for(int ix = 0; ix < box.n.x; ix++)
                fftw_execute_dft(fZ, rho + box.ig({ix, iy, 0}), rho +  box.ig({ix, iy, 0}));

            #pragma omp parallel for collapse(2)
            for(int iz = 0; iz < box.n.z; iz++)
            for(int ix = 0; ix < box.n.x; ix++)
                fftw_execute_dft(fY, rho + box.ig({ix, 0, iz}), rho +  box.ig({ix, 0, iz}));

            #pragma omp parallel for collapse(2)
            for(int iz = 0; iz < box.n.z; iz++)
            for(int iy = 0; iy < box.n.y; iy++)
                fftw_execute_dft(fX, rho + box.ig({0, iy, iz}), rho +  box.ig({0, iy, iz}));
        }

        #pragma omp parallel for collapse(3)
        for(int iz = 0; iz < box.n.z; iz++)
        for(int iy = 0; iy < box.n.y; iy++)
        for(int ix = 0; ix < box.n.x; ix++)
        if((ix != 0)||(iy != 0)||(iz != 0))
        {
            double kx = (ix <= box.n.x >> 1) ? 2*pi*ix*box.invSize.x : -2*pi*(box.n.x - ix)*box.invSize.x;
            double ky = (iy <= box.n.y >> 1) ? 2*pi*iy*box.invSize.y : -2*pi*(box.n.y - iy)*box.invSize.y;
            double kz = (iz <= box.n.z >> 1) ? 2*pi*iz*box.invSize.z : -2*pi*(box.n.z - iz)*box.invSize.z;
            double k = sqrt(kx*kx + ky*ky + kz*kz);

            double filterFactor = 1;
            if(sin2Filter){
                if(abs(kx*box.step.x) > 0.5*pi) filterFactor *= sqr(sin(pi - abs(kx*box.step.x)));
                if(abs(ky*box.step.y) > 0.5*pi) filterFactor *= sqr(sin(pi - abs(ky*box.step.y)));
                if(abs(kz*box.step.z) > 0.5*pi) filterFactor *= sqr(sin(pi - abs(kz*box.step.z)));
            }

            double invk = 1/k; kx *= invk; ky *= invk; kz *= invk; // normalisation of vector k
			double sina = sin(timeStep*k*lightVelocity);
            double cosa = cos(timeStep*k*lightVelocity);
            double R[3][3] = { // rotation matrix
                {cosa + (1 - cosa)*kx*kx,     (1 - cosa)*kx*ky - sina*kz,  (1 - cosa)*kx*kz + sina*ky},
                {(1 - cosa)*ky*kx + sina*kz,  cosa + (1 - cosa)*ky*ky,     (1 - cosa)*ky*kz - sina*kx},
                {(1 - cosa)*kz*kx - sina*ky,  (1 - cosa)*kz*ky + sina*kx,  cosa + (1 - cosa)*kz*kz   }
            };
            complex<double> *F[3]; // pointer to F
            for(int id = 0; id < 3; id++) F[id] = (complex<double>*)fdata[id] + box.ig({ix, iy, iz});
            complex<double> F_plus[3]; // resulting vector F^+
            for(int id = 0; id < 3; id++) F_plus[id] = R[id][0]**F[0] + R[id][1]**F[1] + R[id][2]**F[2];
            if(divergenceCleaning){
                complex<double> Fdotk = F_plus[0]*kx + F_plus[1]*ky + F_plus[2]*kz;
                complex<double> rho_k = *((complex<double>*)rho + box.ig({ix, iy, iz}))*(4*pi*invk);
                if(!initBackground){
                    rhoBackground[box.ig({ix, iy, iz})] = complex<double>(0, 1)*Fdotk - rho_k;
                } else {
                    complex<double> corr = -Fdotk - complex<double>(0, 1)*(rho_k + rhoBackground[box.ig({ix, iy, iz})]);
                    F_plus[0] += corr*kx;
                    F_plus[1] += corr*ky;
                    F_plus[2] += corr*kz;
                }
            }
            for(int id = 0; id < 3; id++) *F[id] = filterFactor*F_plus[id]*invNg; // setting F+ and accounting for the backward FFT factor
        } else
            for(int id = 0; id < 3; id++) *((complex<double>*)fdata[id] + box.ig({0, 0, 0})) *= invNg; // FFT factor for (kx = 0, ky = 0, kz = 0) element

        #pragma omp parallel for collapse(3)
        for(int id = 0; id < 3; id++)
        for(int iz = 0; iz < box.n.z; iz++)
        for(int iy = 0; iy < box.n.y; iy++)
            fftw_execute_dft(bX, fdata[id] + box.ig({0, iy, iz}), fdata[id] + box.ig({0, iy, iz}));

        #pragma omp parallel for collapse(3)
        for(int id = 0; id < 3; id++)
        for(int iz = 0; iz < box.n.z; iz++)
        for(int ix = 0; ix < box.n.x; ix++)
            fftw_execute_dft(bY, fdata[id] + box.ig({ix, 0, iz}), fdata[id] + box.ig({ix, 0, iz}));

        #pragma omp parallel for collapse(3)
        for(int id = 0; id < 3; id++)
        for(int iy = 0; iy < box.n.y; iy++)
        for(int ix = 0; ix < box.n.x; ix++)
            fftw_execute_dft(bZ, fdata[id] + box.ig({ix, iy, 0}), fdata[id] + box.ig({ix, iy, 0}));
        initBackground = true;
    }
    inline void cellSetField(cellInterface &CI, int3 i){
        setGridType(CI, 0); //
        if(getCI_F_Data(CI) == nullptr) getCI_F_Data(CI) = new double[48];
        double* F_data = getCI_F_Data(CI);
        intg cig[8];
        cig[0] = box.ig(i);
        cig[1] = cig[0] + 1; if(unlikely(i.x == box.n.x - 1)) cig[1] -= box.n.x;
        if(box.dim > 1){
            cig[2] = cig[0] + box.n.x; if(unlikely(i.y == box.n.y - 1)) cig[2] -= box.n.y*box.n.x;
            cig[3] = cig[2] + 1; if(unlikely(i.x == box.n.x - 1)) cig[3] -= box.n.x;
            if(box.dim > 2){
                cig[4] = cig[0] + box.n.x*box.n.y; if(unlikely(i.z == box.n.z - 1)) cig[4] -= box.n.z*box.n.y*intg(box.n.x);
                cig[5] = cig[4] + 1; if(unlikely(i.x == box.n.x - 1)) cig[5] -= box.n.x;
                cig[6] = cig[4] + box.n.x; if(unlikely(i.y == box.n.y - 1)) cig[6] -= box.n.y*box.n.x;
                cig[7] = cig[6] + 1; if(unlikely(i.x == box.n.x - 1)) cig[7] -= box.n.x;
            }
        }

        for(int j = 0; j < (1 << box.dim); j++){
            F_data[6*j  ] = Ex(cig[j]);
            F_data[6*j+3] = Bx(cig[j]);
        }
        for(int j = 0; j < (1 << box.dim); j++){
            F_data[6*j  +1] = Ey(cig[j]);
            F_data[6*j+3+1] = By(cig[j]);
        }
        for(int j = 0; j < (1 << box.dim); j++){
            F_data[6*j  +2] = Ez(cig[j]);
            F_data[6*j+3+2] = Bz(cig[j]);
        }
    }
};

struct fieldSubMap64 // structure for a local field buffer (for optimization purposes)
{
    double3 cell0, step, invStep;
    int3 i3;
    int dim;
    double F_data[384]; // em-field at the corners of the cell and other real-value parameters
    intg cig[64];
    fieldSubMap64() {}
    double3& E_(int i){
        return *((double3*)(&F_data[6*i]));
    }
    double3& B_(int i){
        return *((double3*)(&F_data[6*i + 3]));
    }
    //access to em-field: cx, cy, cz are local indices: 0 correspond to cell0, 1 to cell0 + step
    double3& E(int cx = 0, int cy = 0, int cz = 0){return E_((cx+1) + 4*((cy+1) + 4*(cz+1)));}
    double3& B(int cx = 0, int cy = 0, int cz = 0){return B_((cx+1) + 4*((cy+1) + 4*(cz+1)));}

    void downloadField(int3 i, fourierSolver &field){
        i3 = i;
        simulationBox &box(field.box);
        dim = box.dim;
        step = box.step;
        invStep = box.invStep;
        cell0 = {box.min.x + i.x*step.x, box.min.y + i.y*step.y, box.min.z + i.z*step.z};
         if(box.dim == 1){
            for(int sx = -1; sx < 3; sx++){
                int ix = i.x + sx; if(unlikely(ix < 0)) ix += box.n.x; else if(unlikely(ix >= box.n.x)) ix -= box.n.x;
                int ind = sx+1;
                cig[ind] = box.ig({ix, 0, 0});
                E(sx, 0, 0).x = field.Ex(cig[ind]);
                B(sx, 0, 0).x = field.Bx(cig[ind]);
                E(sx, 0, 0).y = field.Ey(cig[ind]);
                B(sx, 0, 0).y = field.By(cig[ind]);
                E(sx, 0, 0).z = field.Ez(cig[ind]);
                B(sx, 0, 0).z = field.Bz(cig[ind]);
            }
        }
        if(field.box.dim == 2){
            for(int sy = -1; sy < 3; sy++)
            {
                int iy = i.y + sy; if(unlikely(iy < 0)) iy += box.n.y; else if(unlikely(iy >= box.n.y)) iy -= box.n.y;
                for(int sx = -1; sx < 3; sx++)
                {
                    int ix = i.x + sx; if(unlikely(ix < 0)) ix += box.n.x; else if(unlikely(ix >= box.n.x)) ix -= box.n.x;
                    int ind = sx+1 + 4*(sy+1);
                    cig[ind] = box.ig({ix, iy, 0});
                    E(sx, sy, 0).x = field.Ex(cig[ind]);
                    B(sx, sy, 0).x = field.Bx(cig[ind]);
                    E(sx, sy, 0).y = field.Ey(cig[ind]);
                    B(sx, sy, 0).y = field.By(cig[ind]);
                    E(sx, sy, 0).z = field.Ez(cig[ind]);
                    B(sx, sy, 0).z = field.Bz(cig[ind]);
                }
            }
        }
        if(field.box.dim == 3){
            for(int sz = -1; sz < 3; sz++)
            {
                int iz = i.z + sz; if(unlikely(iz < 0)) iz += box.n.z; else if(unlikely(iz >= box.n.z)) iz -= box.n.z;
                for(int sy = -1; sy < 3; sy++)
                {
                    int iy = i.y + sy; if(unlikely(iy < 0)) iy += box.n.y; else if(unlikely(iy >= box.n.y)) iy -= box.n.y;
                    for(int sx = -1; sx < 3; sx++)
                    {
                        int ix = i.x + sx; if(unlikely(ix < 0)) ix += box.n.x; else if(unlikely(ix >= box.n.x)) ix -= box.n.x;
                        int ind = sx+1 + 4*((sy+1) + 4*(sz+1));
                        cig[ind] = box.ig({ix, iy, iz});
                        E(sx, sy, sz).x = field.Ex(cig[ind]);
                        B(sx, sy, sz).x = field.Bx(cig[ind]);
                        E(sx, sy, sz).y = field.Ey(cig[ind]);
                        B(sx, sy, sz).y = field.By(cig[ind]);
                        E(sx, sy, sz).z = field.Ez(cig[ind]);
                        B(sx, sy, sz).z = field.Bz(cig[ind]);
                    }
                }
            }
        }
    }
    void uploadField(fourierSolver &field){
        simulationBox &box(field.box);
        if(box.dim == 1)
        for(int sx = -1; sx < 3; sx++){
            int j = sx + 1;
            field.Ex(cig[j]) = E(sx, 0, 0).x;
            field.Bx(cig[j]) = B(sx, 0, 0).x;
            field.Ey(cig[j]) = E(sx, 0, 0).y;
            field.By(cig[j]) = B(sx, 0, 0).y;
            field.Ez(cig[j]) = E(sx, 0, 0).z;
            field.Bz(cig[j]) = B(sx, 0, 0).z;
        }
        if(box.dim == 2)
        for(int sy = -1; sy < 3; sy++)
        for(int sx = -1; sx < 3; sx++){
            int j = sx + 1 + 4*(sy+1);
            field.Ex(cig[j]) = E(sx, sy, 0).x;
            field.Bx(cig[j]) = B(sx, sy, 0).x;
            field.Ey(cig[j]) = E(sx, sy, 0).y;
            field.By(cig[j]) = B(sx, sy, 0).y;
            field.Ez(cig[j]) = E(sx, sy, 0).z;
            field.Bz(cig[j]) = B(sx, sy, 0).z;
        }
        if(box.dim == 3)
        for(int sz = -1; sz < 3; sz++)
        for(int sy = -1; sy < 3; sy++)
        for(int sx = -1; sx < 3; sx++){
            int j = sx + 1 + 4*((sy+1) + 4*(sz+1));
            field.Ex(cig[j]) = E(sx, sy, sz).x;
            field.Bx(cig[j]) = B(sx, sy, sz).x;
            field.Ey(cig[j]) = E(sx, sy, sz).y;
            field.By(cig[j]) = B(sx, sy, sz).y;
            field.Ez(cig[j]) = E(sx, sy, sz).z;
            field.Bz(cig[j]) = B(sx, sy, sz).z;
        }
    }
    void CIC(double3 r, double3 &E, double3 &B); // CIC weighting of E and B field to point r
    void CIC_c(double3 r, double *c, int *cil); // computation of coupling weights for grid nodes of the cell
    void CIC(double *c, int *cil, double3 &E, double3 &B); // CIC weighting of E and B field for a given set of coupling weights
};

void fieldSubMap64::CIC_c(double3 r, double *c, int *cil){
    double3 l_min = cell0;
    int3 l_ind(0, 0, 0);
    if(unlikely(r.x < l_min.x)) {l_min.x -= step.x; l_ind.x = -1;} else if(unlikely(r.x > l_min.x + step.x)) {l_min.x += step.x; l_ind.x = 1;}
    c[1] = (r.x - l_min.x)*invStep.x; c[0] = 1 - c[1];
    if(dim > 1){
        if(unlikely(r.y < l_min.y)) {l_min.y -= step.y; l_ind.y = -1;} else if(unlikely(r.y > l_min.y + step.y)) {l_min.y += step.y; l_ind.y = 1;}
        double wy = (r.y - l_min.y)*invStep.y; double wy_ = 1 - wy;
        for(int i = 0; i < 2; i++) c[i+2] = c[i]*wy;
        for(int i = 0; i < 2; i++) c[i] *= wy_;
    }
    if(dim > 2){
        if(unlikely(r.z < l_min.z)) {l_min.z -= step.z; l_ind.z = -1;} else if(unlikely(r.z > l_min.z + step.z)) {l_min.z += step.z; l_ind.z = 1;}
        double wz = (r.z - l_min.z)*invStep.z; double wz_ = 1 - wz;
        for(int i = 0; i < 4; i++) c[i+4] = c[i]*wz;
        for(int i = 0; i < 4; i++) c[i] *= wz_;
    }
    if(dim == 1){
        cil[0] = (l_ind.x+1) + 4*(1 + 4*1);
        cil[1] = (l_ind.x+2) + 4*(1 + 4*1);
    }
    if(dim == 2){
        cil[0] = (l_ind.x+1) + 4*((l_ind.y+1) + 4*1);
        cil[1] = (l_ind.x+2) + 4*((l_ind.y+1) + 4*1);
        cil[2] = (l_ind.x+1) + 4*((l_ind.y+2) + 4*1);
        cil[3] = (l_ind.x+2) + 4*((l_ind.y+2) + 4*1);
    }
    if(dim == 3){
        cil[0] = (l_ind.x+1) + 4*((l_ind.y+1) + 4*(l_ind.z+1));
        cil[1] = (l_ind.x+2) + 4*((l_ind.y+1) + 4*(l_ind.z+1));
        cil[2] = (l_ind.x+1) + 4*((l_ind.y+2) + 4*(l_ind.z+1));
        cil[3] = (l_ind.x+2) + 4*((l_ind.y+2) + 4*(l_ind.z+1));
        cil[4] = (l_ind.x+1) + 4*((l_ind.y+1) + 4*(l_ind.z+2));
        cil[5] = (l_ind.x+2) + 4*((l_ind.y+1) + 4*(l_ind.z+2));
        cil[6] = (l_ind.x+1) + 4*((l_ind.y+2) + 4*(l_ind.z+2));
        cil[7] = (l_ind.x+2) + 4*((l_ind.y+2) + 4*(l_ind.z+2));
    }
};

void fieldSubMap64::CIC(double *c, int *cil, double3 &E, double3 &B){
    E = {0, 0, 0}; B = {0, 0, 0};
    for(int i = 0; i < (1 << dim); i++){
        E += c[i]*E_(cil[i]);
        B += c[i]*B_(cil[i]);
    }
};

void fieldSubMap64::CIC(double3 r, double3 &E, double3 &B){
    double c[8];
    int cil[8];
    CIC_c(r, c, cil);
    CIC(c, cil, E, B);
};

#endif
