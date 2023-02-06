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

#ifndef EMFIELD_H
#define EMFIELD_H

#include "primitives.h"

struct simulationBox
{
    double3 min, max, invSize, step; // physical limits of the computational region, the vectors of inverse sizes and of space steps
    int3 n; // number of cells alongs each dimension; must be powers of two; for 2D set n.z = 1, for 1D set n.z = n.y = 1 
    size_t ng; // total number of cells
    int dim; // problem dimensionality
    simulationBox(int3 n, double3 min, double3 max): n(n), min(min), max(max), ng(size_t(n.x)*n.y*n.z), 
    invSize(1/(max.x - min.x), 1/(max.y - min.y), 1/(max.z - min.z)),
    step((max.x - min.x)/n.x, (max.y - min.y)/n.y, (max.z - min.z)/n.z)
    {
        dim = 3;
        if(n.z == 1) dim = 2; 
        if((n.y == 1)&&(n.z == 1)) dim = 1; 
    }
    size_t ig(int3 i) // conversion to global index
    {
        return i.x + (i.y + size_t(i.z)*n.y)*n.x; // must be consistent with fftw plans
    }
    int3 ind(size_t ig) // conversion from global index to three indeces
    {
        return {ig%n.x, ig/n.x%n.y, ig/n.x/n.y%n.z};
    } 
    double3 nodeLocation(int3 i) // returns node location
    {
        return {min.x + step.x*(i.x + 0.5), min.y + step.y*(i.y + 0.5), min.z + step.z*(i.z + 0.5)};
    }
};

struct emField // spectral solver for electromagnetic field evolution with zero current
{
    simulationBox box;
    fftw_complex *fdata[3];
    fftw_plan fX, bX, fY, bY, fZ, bZ; // plans for FFT, the first letter standf for forward (f) or backward (b), whereas the second letter indicates the dimension. 
    double invNg; // inverse total number of cells
    
    emField(simulationBox box, int fftw_flags = FFTW_PATIENT): 
    box(box), invNg(1/double(box.ng))
    {        
        for(int i = 0; i < 3; i++) fdata[i] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*box.ng);
        for(int i = 0; i < 3; i++) memset(fdata[i], 0, sizeof(fftw_complex)*box.ng);
        
        fX = fftw_plan_many_dft(1, &box.n.x, 1, fdata[0], &box.n.x, 1, 1, fdata[0], &box.n.x, 1, 1, FFTW_FORWARD, fftw_flags);
        bX = fftw_plan_many_dft(1, &box.n.x, 1, fdata[0], &box.n.x, 1, 1, fdata[0], &box.n.x, 1, 1, FFTW_BACKWARD, fftw_flags);
        fY = fftw_plan_many_dft(1, &box.n.y, 1, fdata[0], &box.n.y, box.n.x, 1, fdata[0], &box.n.y, box.n.x, 1, FFTW_FORWARD, fftw_flags);
        bY = fftw_plan_many_dft(1, &box.n.y, 1, fdata[0], &box.n.y, box.n.x, 1, fdata[0], &box.n.y, box.n.x, 1, FFTW_BACKWARD, fftw_flags);
        fZ = fftw_plan_many_dft(1, &box.n.z, 1, fdata[0], &box.n.z, box.n.x*box.n.y, 1, fdata[0], &box.n.z, box.n.x*box.n.y, 1, FFTW_FORWARD, fftw_flags);
        bZ = fftw_plan_many_dft(1, &box.n.z, 1, fdata[0], &box.n.z, box.n.x*box.n.y, 1, fdata[0], &box.n.z, box.n.x*box.n.y, 1, FFTW_BACKWARD, fftw_flags);
    }
    ~emField()
    {
        fftw_destroy_plan(fX);
        fftw_destroy_plan(bX);
        fftw_destroy_plan(fY);
        fftw_destroy_plan(bY);
        fftw_destroy_plan(fZ);
        fftw_destroy_plan(bZ);        
        for(int i = 0; i < 3; i++) fftw_free(fdata[i]);
    }
    double& Ex(size_t ig) {return *((double*)fdata[0] + 0 + 2*ig);}
    double& Bx(size_t ig) {return *((double*)fdata[0] + 1 + 2*ig);}
    double& Ey(size_t ig) {return *((double*)fdata[1] + 0 + 2*ig);}
    double& By(size_t ig) {return *((double*)fdata[1] + 1 + 2*ig);}
    double& Ez(size_t ig) {return *((double*)fdata[2] + 0 + 2*ig);}
    double& Bz(size_t ig) {return *((double*)fdata[2] + 1 + 2*ig);}
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
			for(int id = 0; id < 3; id++) *F[id] = F_plus[id]*invNg; // setting F+ and accounting for the backward FFT factor
        } else 
            for(int id = 0; id < 3; id++) *((complex<double>*)fdata[id] + box.ig({0, 0, 0})) *= invNg; // backeard FFT factor for (kx = 0, ky = 0, kz = 0) element

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
    }
};
#endif