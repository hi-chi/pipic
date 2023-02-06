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
// Description of the file: File introduces basic elements used in pi-PIC.

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

// All the physical constants and variables are in CGS units.
const double lightVelocity = 29979245800.0;
const double electronCharge = -4.80320427e-10;
const double electronMass = 9.10938215e-28;
const double protonMass = 1.672622964e-24;
#define pi M_PI
using namespace std;

double sgn(double x) 
{
    return (x > 0) - (x < 0);
};

double sqr(double x)
{
	return x*x;
};

struct int3
{ 
    int x, y, z; 
    int3(int x=0, int y=0, int z=0) : x(x), y(y), z(z) {}; 
};

struct double3 // 3D vector with components being double numbers 
{ 
    double x, y, z; 
    double3(double x=0, double y=0, double z=0) : x(x), y(y), z(z) {};
    void normalize()
    {
        double r = sqrt(sqr(x) + sqr(y) + sqr(z));
        if(r > 0) {x /= r; y /= r; z /= r;}
    };
	double norm()
	{
		return sqrt(sqr(x) + sqr(y) + sqr(z));
	}
	double norm2()
	{
		return sqr(x) + sqr(y) + sqr(z);
	}
	string str()
	{
        std::stringstream text;
        text << "<" << x << ", " << y << ", " << z << ">";
        return text.str();
	}
};

double3 operator * (double a, const double3& v)
{
    return double3(v.x * a, v.y * a, v.z * a);
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

template <class myType>
void print(vector<myType> x, string name = "")
{
    if(name != "") cout << name << " = ";
    cout << "{";
    if(x.size() > 0) cout << x[0];
    if(x.size() > 1) for(int i = 1; i < x.size(); i++) cout << ", " << x[i];
    cout << "}";  
};

string to_string(complex<double> a)
{
    return "(" + to_string(a.real()) + " + i*(" + to_string(a.imag()) + "))";
};

double rand_double() // thread-safe generator of a double number from [0, 1] with uniform pdf 
{
    static thread_local mt19937* randGen = nullptr;
    if (!randGen) randGen = new mt19937(clock() + omp_get_thread_num());
    std::uniform_real_distribution<double> distribution(0, 1);
    return distribution(*randGen);
};

int rand_int(int max) // thread-safe generator of an integer number from 0 to 'max' (inclusive) with equal probability for each 
{
    static thread_local mt19937* randGen = nullptr;
    if (!randGen) randGen = new mt19937(clock() + omp_get_thread_num());
    std::uniform_int_distribution<> distribution(0, max);
    return distribution(*randGen);
};

struct logHandler{ // basic structure to handle log messages
    bool logToScreen, logToFile;
    string fileName;
    ofstream file;
    logHandler(string fileName): fileName(fileName), logToFile(true), logToScreen(false){
        file.open (fileName, ios::trunc);
        file << "pi-PIC log file" << endl;
        file.close();
    }
    void message(string text, bool error = false){
        if(logToFile){
            file.open (fileName, ios::app);
            file << text << endl;
            file.close();
        }
        if((logToScreen)||(error)) cout << text << endl;
    }
};

static logHandler pipic_log("pipic_log.txt");

#endif