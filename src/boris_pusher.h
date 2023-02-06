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
// Description of the file: File introduces a basic implementation of handler for Boris pusher.

#ifndef BORIS_PUSHER_H
#define BORIS_PUSHER_H

#include "cic_weighting.h"
#include "particle.h"

struct borisPusher
{
    emField *field; // pointer to the field
    double3 *J; // array for the current
    simulationBox box;
    double timeStep;
    borisPusher(emField *field): field(field), box(field->box)
    {
        J = new double3[size_t(box.n.x)*box.n.y*box.n.z];
    }
    ~borisPusher()
    {
        delete []J;
    }
    void setTimeStep(double timeStep_)
    {
        timeStep = timeStep_;
    }
    void beforeParticleLoop() // sets cuttent to zero everywhere
    {
        memset(J, 0, sizeof(double3)*size_t(box.n.x)*box.n.y*box.n.z);
    }
    void afterParticleLoop() // deposits current E -= 4*pi*timeStep*J
    {
        #pragma omp parallel for collapse(3)
        for(int iz = 0; iz < box.n.z; iz++)
        for(int iy = 0; iy < box.n.y; iy++)
        for(int ix = 0; ix < box.n.x; ix++)
        {
            field->Ex(box.ig({ix, iy, iz})) -= 4*pi*timeStep*J[box.ig({ix, iy, iz})].x;
            field->Ey(box.ig({ix, iy, iz})) -= 4*pi*timeStep*J[box.ig({ix, iy, iz})].y;
            field->Ez(box.ig({ix, iy, iz})) -= 4*pi*timeStep*J[box.ig({ix, iy, iz})].z;
        }
    }
    void processParticle(particle &Particle, particleType &type) // standard Boris pusher
    {

        // computation of weights
        double c[8];
        size_t cig[8];
        CIC(Particle.r, box, c, cig);
        double3 E = {0, 0, 0}, B = {0, 0, 0};
        for(int i = 0; i < (1 << box.dim); i++) // for CIC the size of coupling vector is 2^dim, which is (1 << bos.dim)
        {
            E.x += c[i]*field->Ex(cig[i]);
            B.x += c[i]*field->Bx(cig[i]);
            E.y += c[i]*field->Ey(cig[i]);
            B.y += c[i]*field->By(cig[i]);
            E.z += c[i]*field->Ez(cig[i]);
            B.z += c[i]*field->Bz(cig[i]);
        }

        // Boris push algorithm
        double3 u0 = Particle.p + 0.5*type.charge*timeStep*E;
        double gamma = sqrt(1 + u0.norm2()/sqr(type.mass*lightVelocity));
        double3 t = (type.charge*timeStep/(2*gamma*type.mass*lightVelocity))*B;
        double3 u1 = u0 + cross(u0, t);
        double3 s = (2/(1 + t.norm2()))*t;
        Particle.p = u0 + cross(u1, s) + 0.5*type.charge*timeStep*E;
        gamma = sqrt(1 + Particle.p.norm2()/sqr(type.mass*lightVelocity));
        Particle.r = Particle.r + (timeStep/(type.mass*gamma))*Particle.p;
        
        //add contribution to current
        double Vg = box.step.x*box.step.y*box.step.z;
        for(int i = 0; i < (1 << box.dim); i++) J[cig[i]] = J[cig[i]] + c[i]*Particle.w*(type.charge/(Vg*type.mass*gamma))*Particle.p;
    }
};
#endif