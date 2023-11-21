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
// Description of the file: In this file all the introduced structures are combined.

#ifndef PIPIC_H
#define PIPIC_H

#include "fourier_boris_solver.h"
#include "ec_solver.h"
#include "ec2_solver.h"

struct pipic
{
    int nx, ny, nz;
    double XMin, XMax, YMin, YMax, ZMin, ZMax; 
    simulationBox box;
    field_solver *Field;
    ensemble *Ensemble;
    pic_solver *Solver;
    bool reportPerformance;
    pipic(string solverName, int nx, double XMin, double XMax, int ny = 1, double YMin = -0.5, double YMax = 0.5, int nz = 1, double ZMin = -0.5, double ZMax = 0.5): 
    nx(nx), ny(ny), nz(nz), XMin(XMin), XMax(XMax), YMin(YMin), YMax(YMax), ZMin(ZMin), ZMax(ZMax), 
    box(int3(nx, ny, nz), double3(XMin, YMin, ZMin), double3(XMax, YMax, ZMax)) 
    {
        Solver = nullptr;
        if(solverName == "fourier_boris") Solver = new fourier_boris_solver(box);
        if(solverName == "ec") Solver = new ec_solver(box);
        if(solverName == "ec2") Solver = new ec2_solver(box);
        if(Solver == nullptr){ pipic_log.message("pi-PIC init() error: unknown solver '" + solverName + "'.", true); exit(0);}
        Field = Solver->Field;
        Ensemble = Solver->Ensemble;
        reportPerformance = true;
    }
    ~pipic(){
        delete Solver;
    }
    void pyAddParticles(string typeName, int totalNumber, double typeCharge, double typeMass, double temperature, int64_t density, int64_t dataDouble = 0, int64_t dataInt = 0){
        Ensemble->placeParticles(typeName, totalNumber, typeCharge, typeMass, temperature, density, dataDouble, dataInt);
    }
    int getNumberOfParticles(){
        return Ensemble->totalNumberOfParticles;
    }
    void pyParticleLoop(string typeName, int64_t handler, int64_t dataDouble, int64_t dataInt)
    {
        void(*handler_)(double*, double*, double*, unsigned long long int*, double*, int*) = (void(*)(double*, double*, double*, unsigned long long int*, double*, int*))handler;
        double* dataDouble_ = nullptr; if(dataDouble != 0) dataDouble_ = (double*)dataDouble;
        int* dataInt_ = nullptr; if(dataInt != 0) dataInt_ = (int*)dataInt;

        int typeIndex = Ensemble->getTypeIndex(typeName);
        for(ensemble::nonOmpIterator iP = Ensemble->begin(typeIndex); iP < Ensemble->end(); iP++){
            particle *P = &*iP;
            handler_(&(P->r.x), &(P->p.x), &(P->w), &(P->id), dataDouble_, dataInt_);
        }
    }
    void pyFieldLoop(int64_t handler, int64_t dataDouble = 0, int64_t dataInt = 0, bool useOmp = false){
        Field->fieldLoop(handler, dataDouble, dataInt, useOmp);
    }
    void pyCustomFieldLoop(int numberOfIterations, int64_t it2coord, int64_t field2data, int64_t dataDouble = 0, int64_t dataInt = 0){
        Field->customFieldLoop(numberOfIterations, it2coord, field2data, dataDouble, dataInt);
    }
    void pyAdvance(double timeStep, int numberOfIterations = 1){
        for(int iIt = 0; iIt < numberOfIterations; iIt++) {
            Solver->advance(timeStep);
            Ensemble->Manager.iterationEnd(box.ng);
            pipic_log.saveBuffer();
        }
    }
    void pyLogPolicy(bool logToFile = true, bool logToScreen = true){
        pipic_log.logToFile = logToFile;
        pipic_log.logToScreen = logToScreen;
    }
    void pyFourierSolverSettings(int divergenceCleaning_ = -1, int sin2_kFilter = -1){
        if(Field->type != "fourier") {cout << "pi-PIC error: fourierSolverSettings() is only available when Fourier solver is in use."; return;}
        fourierSolver *field = dynamic_cast<fourierSolver*>(Field);
        if(divergenceCleaning_ == 1) field->enableDivergenceCleaning();
        if(divergenceCleaning_ == 0) field->divergenceCleaning = false;
        if(sin2_kFilter == 1) field->sin2Filter = true;
        if(sin2_kFilter == 0) field->sin2Filter = false;
    }
    void setRngGenSeed(int seed){
        Ensemble->RndGen.setRngSeed(seed);
    }
    int getTypeIndex(string typeName){
        return Ensemble->getTypeIndex(typeName);
    }
    void addHandler(string name, string subject, int64_t handler, int64_t dataDouble = 0, int64_t dataInt = 0){
        double* dataDouble_ = nullptr; if(dataDouble != 0) dataDouble_ = (double*)dataDouble;
        int* dataInt_ = nullptr; if(dataInt != 0) dataInt_ = (int*)dataInt;
        Ensemble->Manager.addCellHandler(name, subject, handler, dataDouble_, dataInt_);
    }
};

#endif