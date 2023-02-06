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
// Description of the file: File introduces an implementation of energy-conserving pusher with coupling
// based on CIC weighting for mid point determined according to CIC weighting.

#ifndef EC_PUSHER_H
#define EC_PUSHER_H

#include "cic_weighting.h"
#include "particle.h"

struct ecPusher
{
    emField *field; // pointer to the field
    simulationBox box;
    double timeStep;
    ecPusher(emField *field): field(field), box(field->box) {}
    void setTimeStep(double timeStep_)
    {
        timeStep = timeStep_;
    }
    void beforeParticleLoop() {}
    void afterParticleLoop() {}
    void processParticle(particle &Particle, particleType &type)
    {
        double m = type.mass*Particle.w;
        double q = type.charge*Particle.w;
        double3 p = (1/(type.mass*lightVelocity))*Particle.p; // dimensionless momentum, like in the paper

        double inv_mc = 1/(m*lightVelocity);
        double gamma = sqrt(1 + p.norm2());
        double inv_gamma = 1/gamma;
        double Vg = box.step.x*box.step.y*box.step.z;
        double inv_Vg = 1/Vg;
    
        // computation of weights and initial field vectors
        double c[8]; size_t cig[8];
        double3 dr = (0.5*timeStep*inv_gamma*lightVelocity)*p;
        double c_dr = 1; // coefficient to set a cap on the location where the field is considered (should be no more than a cell away to run omp)
        if(abs(dr.x) > box.step.x) c_dr = box.step.x/abs(dr.x);
        if(box.n.y > 1)if(abs(dr.y) > box.step.y)if(c_dr > box.step.y/abs(dr.y)) c_dr = box.step.y/abs(dr.y);
        if(box.n.z > 1)if(abs(dr.z) > box.step.z)if(c_dr > box.step.z/abs(dr.z)) c_dr = box.step.z/abs(dr.z);
        double3 r_mid = Particle.r + c_dr*dr;
        CIC(r_mid, box, c, cig);
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

        // rotation in magnetic field:
        double B_norm = B.norm();
        if(B_norm > 0) B = (1/B_norm)*B;
        double a = (timeStep*inv_gamma*q*inv_mc)*B_norm;
        double cosa = cos(a);
        double sina = sin(a);
        p = cosa*p + (1 - cosa)*dot(p, B)*B + sina*cross(p, B);

        // E-p coupling
        double xi = 0;
        for(int i = 0; i < (1 << box.dim); i++) xi += sqr(c[i]);
        double sqrt_kappa = sqrt(4*pi*sqr(q)*lightVelocity*inv_mc*inv_Vg*inv_gamma*xi);
        double b = sqrt_kappa*timeStep;
        double cosb = cos(b);
        double sinb = sin(b);
        double3 pt = q*inv_mc*E;
        double3 pt_ = cosb*pt - sqrt_kappa*sinb*p;
        p = cosb*p + (1/sqrt_kappa)*sinb*pt;
        double3 dE = (1/xi)*((m*lightVelocity/q)*pt_ - E);

        double E2_b = 0; // \sum E^2 before update
        for(int i = 0; i < (1 << box.dim); i++) E2_b += (sqr(field->Ex(cig[i])) + sqr(field->Ey(cig[i])) + sqr(field->Ez(cig[i])));
        //update field
        for(int i = 0; i < (1 << box.dim); i++)
        {
            field->Ex(cig[i]) += c[i]*dE.x;
            field->Ey(cig[i]) += c[i]*dE.y;
            field->Ez(cig[i]) += c[i]*dE.z;
        }
        double E2_a = 0; // \sum E^2 after update
        for(int i = 0; i < (1 << box.dim); i++) E2_a += (sqr(field->Ex(cig[i])) + sqr(field->Ey(cig[i])) + sqr(field->Ez(cig[i])));
        
        double sigma = sqrt((sqr((E2_b - E2_a)*Vg/(8*pi*m*sqr(lightVelocity)) + gamma) - 1)/p.norm2());
        Particle.p = (sigma*type.mass*lightVelocity)*p;

        Particle.r = Particle.r - (Vg/(4*pi*q))*dE;
    }
};

#endif