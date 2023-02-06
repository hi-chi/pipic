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
// Description of the file: File introduces structures for particles, their types and random generation.

#ifndef PARTICLE_H
#define PARTICLE_H

#include "emfield.h"

struct particle 
{
    double3 r, p; // coordinate and momentum
    double w; // weight
    unsigned long long int id; // id
};

struct particleIdGenerator // structure for issuing ids
{
    unsigned long long int threadRange; // range of values for a single thread
    vector<unsigned long long int> counter; // counter for each thread
    particleIdGenerator():
    threadRange(numeric_limits<unsigned long long int>::max()/omp_get_max_threads()),
    counter(omp_get_max_threads(), 0)
    {}
    unsigned long long int issue()
    {
        counter[omp_get_thread_num()]++;
        return omp_get_thread_num()*threadRange + counter[omp_get_thread_num()] - 1;
    }
};

double3 generateMomentum(double mass, double temperature) // thread-safe generation of random momentum for given mass and temperature
{
    static thread_local mt19937* randGen = nullptr;
    if (!randGen) randGen = new mt19937(clock() + omp_get_thread_num());
    std::normal_distribution<double> distribution(0.0, 1.0);
    double3 p = {distribution(*randGen), distribution(*randGen), distribution(*randGen)};
    p = sqrt(2*mass*temperature/3)*p;
    p = sqrt(1 + 0.25*p.norm2()/sqr(mass*lightVelocity))*p; // relativistic correction (\gamma(p_corrected) mc^2 = p^2/2m)
    return p;
}

struct particleType
{
    string name;
    double charge, mass;
    particleType(string name, double charge, double mass) : name(name), charge(charge), mass(mass) {} 
};

#endif