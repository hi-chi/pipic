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
// Description: Here generic privimitive objects and functions are introduced.

#ifndef PRIMITIVES_H
#define PRIMITIVES_H

#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <random>
#include <algorithm>
#include <string>
#include <cstring>
#include <fstream>
#include <stdio.h>
#include <cstdint>
#include <pthread.h>
#include <omp.h>
#include <fftw3.h>
#include <complex>
#include <algorithm>
#include <random>
#include <thread>

// All the physical constants and variables are in CGS units.
namespace constants {
    const double lightVelocity = 29979245800.0;
    const double electronCharge = -4.80320427e-10;
    const double electronMass = 9.10938215e-28;
    const double protonMass = 1.672622964e-24;
    const double pi = 3.14159265358979323846;
    const double hbar = 1.0545716818e-27;
}

using namespace constants;
using namespace std;

#define likely(x)       __builtin_expect((x), 1)
#define unlikely(x)     __builtin_expect((x), 0)

typedef int intg; // type used to enumerate grid cells
//typedef unsigned long long int intg; // use if the grid is too large for int

inline double sgn(double x){
    return (x > 0) - (x < 0);
};

inline double sqr(double x){
	return x*x;
};

struct int3 
{ 
    int x, y, z; 
    int3(int x=0, int y=0, int z=0) : x(x), y(y), z(z) {}; 
};

struct double3
{ 
    double x, y, z; 
    double3(double x=0, double y=0, double z=0) : x(x), y(y), z(z) {};
    void normalize(){
        double r = sqrt(sqr(x) + sqr(y) + sqr(z));
        if(r > 0) {double inv_r = 1/r; x *= inv_r; y *= inv_r; z *= inv_r;}
    };
	double norm(){
		return sqrt(sqr(x) + sqr(y) + sqr(z));
	}
	double norm2(){
		return sqr(x) + sqr(y) + sqr(z);
	}
    double3& operator+=(const double3 a){
        this->x += a.x;
        this->y += a.y;
        this->z += a.z;
        return *this;
    }
	string str(){
        std::stringstream text;
        text << "<" << x << ", " << y << ", " << z << ">";
        return text.str();
	}
};

double3 operator * (double a, const double3& v)
{
    return double3(v.x * a, v.y * a, v.z * a);
};

double3 operator * (const double3& v, double a)
{
    return double3(v.x * a, v.y * a, v.z * a);
};

double3 operator / (const double3& v, double a)
{
    return double3(v.x, v.y, v.z) * (1/a);
};

double3 operator + (const double3& a, const double3& b)
{
    return double3(a.x + b.x, a.y + b.y, a.z + b.z);
};

double3 operator - (const double3& a, const double3& b)
{
    return double3(a.x - b.x, a.y - b.y, a.z - b.z);
};

double3 cross(const double3& a, const double3& b)
{
    return {a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x};
};

double dot(const double3& a, const double3& b)
{
    return a.x*b.x + a.y*b.y + a.z*b.z;
};

double sqr(const double3& a)
{
    return a.x*a.x + a.y*a.y + a.z*a.z;
};

template <class Type>
void print(vector<Type> x, string name = "")
{
    if(name != "") cout << name << " = ";
    cout << "{";
    if(x.size() > 0) cout << x[0];
    if(x.size() > 1) for(int i = 1; i < x.size(); i++) cout << ", " << x[i];
    cout << "}";  
};

string to_string(double x, int precision){
    std::stringstream text;
    text.precision(precision); text << x;
    return text.str();
}

string to_string(complex<double> a)
{
    return "(" + to_string(a.real()) + " + i*(" + to_string(a.imag()) + "))";
};

#endif