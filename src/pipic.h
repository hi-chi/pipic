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
// Description of the file: In this file all the introduced structures are combined to enable the main 
// PIC loop, as well as auxilary possibilities for running loops from Python.

#ifndef PIPIC_H
#define PIPIC_H

#include "chronometer.h"
#include "ensemble.h"
#include "boris_pusher.h"
#include "ec_pusher.h"

struct piPIC // base class for various implementations
{
    emField *Field;
    ensemble *Ensemble;
    bool reportPerformance;
    int nx, ny, nz;
    double XMin, XMax, YMin, YMax, ZMin, ZMax; 
    piPIC(int3 n, double3 min, double3 max){
        nx = n.x; ny = n.y; nz = n.z;
        XMin = min.x; XMax = max.x; YMin = min.y; YMax = max.y; ZMin = min.z; ZMax = max.z; 
        Field = new emField(simulationBox(n, min, max));
        Ensemble = new ensemble(Field->box);
        reportPerformance = true;
    }
    ~piPIC(){
        delete Field, Ensemble;
    }
    void pyAddParticles(string typeName, int totalNumber, double typeCharge, double typeMass, double temperature, int64_t density, int64_t dataDouble = 0, int64_t dataInt = 0){
        Ensemble->placeParticles(typeName, totalNumber, typeCharge, typeMass, temperature, density, dataDouble, dataInt);
    }
    int getNumberOfParticles()
    {
        return Ensemble->totalNumberOfParticles;
    }
    void pyParticleLoop(string typeName, int64_t handler, int64_t dataDouble, int64_t dataInt)
    {
        void(*handler_)(double*, double*, double*, unsigned long long int*, double*, int*) = (void(*)(double*, double*, double*, unsigned long long int*, double*, int*))handler;
        double* dataDouble_ = nullptr; if(dataDouble != 0) dataDouble_ = (double*)dataDouble;
        int* dataInt_ = nullptr; if(dataInt != 0) dataInt_ = (int*)dataInt;

        int typeIndex = Ensemble->getTypeIndex(typeName);
        double *r; double *p; double *w; unsigned long long int *id;
        for(ensemble::nonOmpIterator iP = Ensemble->begin(typeIndex); iP < Ensemble->end(); iP++)
        {
            particle *P = &*iP;
            handler_(&(P->r.x), &(P->p.x), &(P->w), &(P->id), dataDouble_, dataInt_);
            if(P->w == 0) iP.removeCurrentParticle();
        }
    }
    void pyFieldLoop(int64_t handler, int64_t dataDouble = 0, int64_t dataInt = 0, bool useOmp = false)
    {
        void(*handler_)(int*, double*, double*, double*, double*, int*) = (void(*)(int*, double*, double*, double*, double*, int*))handler;
        double* dataDouble_ = nullptr; if(dataDouble != 0) dataDouble_ = (double*)dataDouble;
        int* dataInt_ = nullptr; if(dataInt != 0) dataInt_ = (int*)dataInt;
        if(useOmp){
            #pragma omp parallel for collapse(3)
            for(int iz = 0; iz < nz; iz++)
            for(int iy = 0; iy < ny; iy++)
            for(int ix = 0; ix < nx; ix++)
            {
                size_t ig = Field->box.ig({ix, iy, iz});
                double3 r3 = Field->box.nodeLocation({ix, iy, iz});
                double r[3] = {r3.x, r3.y, r3.z};
                int ind[3] = {ix, iy, iz};
                double E[3], B[3];
                E[0] = Field->Ex(ig); B[0] = Field->Bx(ig);
                E[1] = Field->Ey(ig); B[1] = Field->By(ig);
                E[2] = Field->Ez(ig); B[2] = Field->Bz(ig);
                handler_(ind, r, E, B, dataDouble_, dataInt_);
                Field->Ez(ig) = E[2]; Field->Bz(ig) = B[2];
                Field->Ey(ig) = E[1]; Field->By(ig) = B[1];
                Field->Ex(ig) = E[0]; Field->Bx(ig) = B[0];
            }
        } else {
            for(int iz = 0; iz < nz; iz++)
            for(int iy = 0; iy < ny; iy++)
            for(int ix = 0; ix < nx; ix++)
            {
                size_t ig = Field->box.ig({ix, iy, iz});
                double3 r3 = Field->box.nodeLocation({ix, iy, iz});
                double r[3] = {r3.x, r3.y, r3.z};
                int ind[3] = {ix, iy, iz};
                double E[3], B[3];
                E[0] = Field->Ex(ig); B[0] = Field->Bx(ig);
                E[1] = Field->Ey(ig); B[1] = Field->By(ig);
                E[2] = Field->Ez(ig); B[2] = Field->Bz(ig);
                handler_(ind, r, E, B, dataDouble_, dataInt_);
                Field->Ez(ig) = E[2]; Field->Bz(ig) = B[2];
                Field->Ey(ig) = E[1]; Field->By(ig) = B[1];
                Field->Ex(ig) = E[0]; Field->Bx(ig) = B[0];
            }
        }
    }
    void pyCustomFieldLoop(int numberOfIterations, int64_t it2coord, int64_t field2data, int64_t dataDouble = 0, int64_t dataInt = 0)
    {
        void(*it2coord_)(int*, double*, double*, int*) = (void(*)(int*, double*, double*, int*))it2coord;
        void(*field2data_)(int*, double*, double*, double*, double*, int*) = (void(*)(int*, double*, double*, double*, double*, int*))field2data;
        double* dataDouble_ = nullptr; if(dataDouble != 0) dataDouble_ = (double*)dataDouble;
        int* dataInt_ = nullptr; if(dataInt != 0) dataInt_ = (int*)dataInt;
        for(int it_ = 0; it_ < numberOfIterations; it_++)
        {
            int it[1]; it[0] = it_;
            double3 r;
            it2coord_(it, (double*)&r, dataDouble_, dataInt_);
            double3 E({0, 0, 0}), B({0, 0, 0});
            double c[8];
            size_t cig[8];
            CIC(r, Field->box, c, cig);
            for(int i = 0; i < (1 << Field->box.dim); i++) // for CIC the size of coupling vector is 2^dim
            {
                E.x += c[i]*Field->Ex(cig[i]);
                B.x += c[i]*Field->Bx(cig[i]);
                E.y += c[i]*Field->Ey(cig[i]);
                B.y += c[i]*Field->By(cig[i]);
                E.z += c[i]*Field->Ez(cig[i]);
                B.z += c[i]*Field->Bz(cig[i]);
            }
            field2data_(it, (double*)&r, (double*)&E, (double*)&B, dataDouble_, dataInt_);
        }
    }
    void pyAdvance(double timeStep_, int numberOfIterations_ = 1)
    {
        advance(timeStep_, numberOfIterations_);
    }
    void pySetShufflingPolicy(double cellSuffleProbability, bool doubleLoop)
    {
        Ensemble->shuffle_cellShuffleProbability = cellSuffleProbability;
        Ensemble->shuffle_doubleLoop = doubleLoop;
    }
    virtual void beforeParticleLoop(double timeStep) = 0;
    virtual void particleLoop() = 0;
    virtual void afterParticleLoop() {};
    void advance(double timeStep, int numberOfIterations = 1)
    {
        for(int iIt = 0; iIt < numberOfIterations; iIt++)
        {
            beforeParticleLoop(timeStep);
            double timeParticles_s = 0, timeField_s = 0; chronometer Chronometer;
            Chronometer.start();
            if(Ensemble->totalNumberOfParticles > 0) particleLoop();
            Chronometer.stop(); timeParticles_s += Chronometer.getTime_s(); Chronometer.reset();
            Chronometer.start();
            afterParticleLoop();
            Field->advance(timeStep);
            Chronometer.stop(); timeField_s += Chronometer.getTime_s(); Chronometer.reset();
            if(reportPerformance) {
                string statement = "iteration completed in " + to_string(timeParticles_s + timeField_s) + " s";
                if(Ensemble->totalNumberOfParticles > 0) statement += + " (particles: " + to_string(int(100*timeParticles_s/(timeParticles_s + timeField_s))) + "\%); " 
                                                                      + to_string((1e+9)*timeParticles_s/Ensemble->totalNumberOfParticles) + " ns per particle";
                statement += "; " + to_string((1e+9)*timeField_s/(Field->box.n.x*Field->box.n.y*Field->box.n.z)) + " ns per cell.";
                pipic_log.message(statement);
            }
        }
    }
    void pyLogPolicy(bool logToFile = true, bool logToScreen = true)
    {
        pipic_log.logToFile = logToFile;
        pipic_log.logToScreen = logToScreen;
    }
};

struct pipic_boris : public piPIC
{
    borisPusher *Pusher;
    pipic_boris(int Nx, double Xmin, double Xmax, int Ny = 0.5, double Ymin = -0.5, double Ymax = 0.5, int Nz = 0.5, double Zmin = -0.5, double Zmax = 0.5) : 
    piPIC({Nx, Ny, Nz}, {Xmin, Ymin, Zmin},{Xmax, Ymax, Zmax})
    {
        Pusher = new borisPusher(Field);
        Ensemble->shuffle_cellShuffleProbability = 0;
        Ensemble->shuffle_doubleLoop = false;
    }
    ~pipic_boris()
    {
        delete Pusher;
    }
    void beforeParticleLoop(double timeStep)
    {
        Pusher->setTimeStep(timeStep); 
        Pusher->beforeParticleLoop();
    }
    void particleLoop()
    {
        Ensemble->ompLoop<borisPusher>(Pusher);
    }
    void afterParticleLoop()
    { 
        Pusher->afterParticleLoop();
    }
};

struct pipic_ec: public piPIC
{
    ecPusher *Pusher;
    pipic_ec(int Nx, double Xmin, double Xmax, int Ny = 0.5, double Ymin = -0.5, double Ymax = 0.5, int Nz = 0.5, double Zmin = -0.5, double Zmax = 0.5) : 
    piPIC({Nx, Ny, Nz}, {Xmin, Ymin, Zmin},{Xmax, Ymax, Zmax})
    {
        Pusher = new ecPusher(Field);
        Ensemble->shuffle_cellShuffleProbability = 1;
        Ensemble->shuffle_doubleLoop = true;
    }
    ~pipic_ec(){
        delete Pusher;
    }
    void beforeParticleLoop(double timeStep)
    {
        Pusher->setTimeStep(timeStep); 
        Pusher->beforeParticleLoop();
    }
    void particleLoop()
    {
        Ensemble->ompLoop<ecPusher>(Pusher);
    }
    void afterParticleLoop()
    { 
        Pusher->afterParticleLoop();
    }
};

#endif