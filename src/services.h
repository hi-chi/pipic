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
// Description: Here essential structures to serve the implementation are defined.

#ifndef SERVICES_H
#define SERVICES_H

#include <chrono>
#include <iomanip>
#include "interfaces.h"

struct chronometer
{
    std::chrono::time_point<std::chrono::high_resolution_clock> t1;
    std::chrono::time_point<std::chrono::high_resolution_clock> t2;
    double total_ms;
    chronometer(): total_ms(0){}
    void start(){
        t1 = std::chrono::high_resolution_clock::now();
    }
    void stop(){
        t2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> ms_double = t2 - t1;
        total_ms += ms_double.count();
    }
    double getTime_s(){
        return 0.001*total_ms;
    }
};

struct pipicLog // basic structure to handle log messages
{
    string fileName;
    bool logToFile, logToScreen;
    ofstream file;
    vector<string> messages;
    pipicLog(string fileName): fileName(fileName), logToFile(true), logToScreen(false){
        file.open (fileName, ios::trunc);
        file << "pi-PIC log file" << endl;
        file.close();
    }
    void message(string text, bool error = false, bool immediately = false){
        if(logToFile) messages.push_back(text);
        if(immediately) saveBuffer();
        if((logToScreen)||(error)) cout << text << endl;
    }
    void saveBuffer(){
        if(logToFile){
            file.open (fileName, ios::app);
            for(int i = 0; i < int(messages.size()); i++) file << messages[i] << endl;
            file.close();
            messages.clear();
        }
    }
};

static pipicLog pipic_log("pipic_log.txt");

static thread_local mt19937* randGen = nullptr;

struct rndGen
{
    mt19937** localGen;
    int rngSeed;
    int nx;
    void setSeeds(){
        mt19937 rng(rngSeed);
        std::uniform_int_distribution<> uniform(std::numeric_limits<int>::min(), std::numeric_limits<int>::max());
        for(int i = 0; i < nx; i++) localGen[i]->seed(uniform(rng));
    }
    rndGen(int nx): nx(nx){
        rngSeed = 0; // can be set to clock() if reproducibility is not needed.
        pipic_log.message("rng_seed = " + to_string(rngSeed));
        localGen = new mt19937*[nx];
        for(int i = 0; i < nx; i++) localGen[i] = new mt19937;
        setSeeds();
    }
    void setRngSeed(int newSeed){
        rngSeed = newSeed;
        pipic_log.message("Setting rng_seed = " + to_string(rngSeed));
        setSeeds();
    }
    ~rndGen(){
        for(int i = 0; i < nx; i++) delete localGen[i];
        delete []localGen;
    }
    void setTreadLocalRng(int ix){
        randGen = localGen[ix];
    }
};

double rand_double(double min, double max){ // generates random double with uniform distribution in [min, max)
    if(unlikely(!randGen)) {cout << "pi-PIC error: rndGen call before assignment" << endl; exit(0);}
    std::uniform_real_distribution<double> distribution(min, max);
    return distribution(*randGen);
};

int rand_int(int min, int max){ // generates random integer with uniform distribution from min to max (including both min and max)
    if(unlikely(!randGen)) {cout << "pi-PIC error: rndGen call before assignment" << endl; exit(0);}
    std::uniform_int_distribution<> distribution(min, max);
    return distribution(*randGen);
};

double3 generateMomentum(double mass, double temperature){ // temperature is defined as two-thirds of the average energy
    if(unlikely(!randGen)) {cout << "pi-PIC error: rndGen call before assignment" << endl; exit(0);}
    std::normal_distribution<double> distribution(0.0, 1.0);
    double3 p = {distribution(*randGen), distribution(*randGen), distribution(*randGen)};
    p = sqrt(mass*temperature)*p;
    p = sqrt(1 + 0.25*p.norm2()/sqr(mass*lightVelocity))*p; // relativistic correction, (\gamma(p_corrected) - 1) mc^2 = p^2/2m
    return p;
};

unsigned long long int generateID(){ // thread-safe generation of a particle ID
    static thread_local unsigned long long int IDcounter = omp_get_thread_num()*(numeric_limits<unsigned long long int>::max()/omp_get_max_threads());
    IDcounter++;
    return IDcounter - 1;
};

inline void placePeriodic(double &v, double vMin, double vMax){
    if(unlikely(v < vMin)){
        v += vMax - vMin;
        if(unlikely(v < vMin)) v -= (vMax - vMin)*floor(v/(vMax - vMin));
    }
    if(unlikely(v > vMax)){
        v -= vMax - vMin;
        if(unlikely(v > vMax)) v -= (vMax - vMin)*floor(v/(vMax - vMin));
    }
};

struct cellHandler
{
    string name;
    string subject;
    bool actOnCell;
    vector<bool> actOn; // vector indicating which particle the handler acts on
    double *dataDouble;
    int *dataInt;
    bool firstRun;
    void(*handler_)(int*, double*, double*, double*, double*, double*, int*);

    vector<chronometer> Chronometer_particles, Chronometer_actOnCell;
    vector<unsigned long long int> latest_numberProcessed;
    double best_partilesTime_ns, best_cellTime_ns, total_partilesTime_s, total_cellTime_s, total_numberProcessed;

    cellHandler(string name, string subject, int64_t handler, double* dataDouble, int* dataInt):
    name(name), subject(subject), dataDouble(dataDouble), dataInt(dataInt), firstRun(true),
    Chronometer_particles(omp_get_max_threads()), Chronometer_actOnCell(omp_get_max_threads()),
    latest_numberProcessed(omp_get_max_threads(), 0),
    best_partilesTime_ns(numeric_limits<double>::max()), best_cellTime_ns(numeric_limits<double>::max()),
    total_partilesTime_s(0), total_cellTime_s(0), total_numberProcessed(0)
    {
        handler_ = (void(*)(int*, double*, double*, double*, double*, double*, int*))handler;
    }
    void handle(cellInterface *CI){// must be thread-safe
        if(CI->particleTypeIndex == -1) Chronometer_actOnCell[omp_get_thread_num()].start();
        else Chronometer_particles[omp_get_thread_num()].start();
        handler_(CI->I, CI->D, CI->F_data, CI->P_data, CI->NP_data, dataDouble, dataInt);
        if(CI->particleTypeIndex == -1) Chronometer_actOnCell[omp_get_thread_num()].stop();
        else {
            Chronometer_particles[omp_get_thread_num()].stop();
            latest_numberProcessed[omp_get_thread_num()] += CI->particleSubsetSize;
        }
    }
    void preLoop(vector<string> typeName){
        actOn.resize(typeName.size());
        for(int i = 0; i < int(actOn.size()); i++) actOn[i] = false;
        actOnCell = false;
        char* pch;
        char* source = new char[subject.size() + 1];
        for(int i = 0; i < int(subject.size()); i++) source[i] = subject[i];
        source[subject.size()] = '\0';
        pch = strtok(source, " ,.-");
        while (pch != NULL){
            bool detected = false;
            if(string(pch) == "all_types") {
                for(int i = 0; i < int(actOn.size()); i++) actOn[i] = true;
                detected = true;
            }
            if(string(pch) == "cells") {
                actOnCell = true;
                detected = true;
            }
            for(int i = 0; i < int(actOn.size()); i++) if(string(pch) == typeName[i].c_str()) {
                actOn[i] = true;
                detected = true;
            }
            if((!detected)&&(firstRun)){
                string message = "Warning: The subject of cellHandler <" + name + "> includes unknown entry <" + string(pch) + ">.";
                pipic_log.message(message);
                cout << "pi-PIC " << message << endl;
            }
            pch = strtok(NULL, " ,.-");
        }
        firstRun = false;
    }
    void postLoop(intg ng){
        double time_actOncell_ms = 0, time_particles_ms = 0;
        unsigned long long int total_latestProcessed = 0;
        for(int iTh = 0; iTh < omp_get_max_threads(); iTh++){
            time_actOncell_ms += Chronometer_actOnCell[iTh].total_ms;
            Chronometer_actOnCell[iTh].total_ms = 0;
            time_particles_ms += Chronometer_particles[iTh].total_ms;
            Chronometer_particles[iTh].total_ms = 0;
            total_latestProcessed += latest_numberProcessed[iTh];
            latest_numberProcessed[iTh] = 0;
        }
        if(actOnCell) {
            best_cellTime_ns = min(best_cellTime_ns, (1e+6)*time_actOncell_ms/double(ng));
            total_cellTime_s += 0.001*time_actOncell_ms;
        }
        if(total_latestProcessed > 0) {
            best_partilesTime_ns = min(best_partilesTime_ns, (1e+6)*time_particles_ms/double(total_latestProcessed));
            total_partilesTime_s += 0.001*time_particles_ms;
            total_numberProcessed += total_latestProcessed;
        }
    }
};


struct handlerManager
{
    vector<cellHandler*> Handler;
    bool reportPerformance;
    // all times are in seconds, unless suffix "_ns" is explicitly added to the name of variable
    double latestEnsembleTime, totalEnsembleTime, latestParticleUpdates, totalParticleUpdates;
    double latestFieldTime, totalFieldTime, totalCellUpdates, numberOfCells;
    double latest_av_ppc, latest_av_cmr; // average number of particles per cell and average inter-cell migration rate
    double best_cellUpdateTime_ns, best_particleUpdateTime_ns;
    handlerManager(double numberOfCells): reportPerformance(true),
    totalEnsembleTime(0), totalParticleUpdates(0),
    totalFieldTime(0), totalCellUpdates(0), numberOfCells(numberOfCells),
    best_cellUpdateTime_ns(numeric_limits<double>::max()),
    best_particleUpdateTime_ns(numeric_limits<double>::max())
    {}
    void iterationEnd(intg ng){
        for(int ih = 0; ih < int(Handler.size()); ih++) Handler[ih]->postLoop(ng);

        best_cellUpdateTime_ns = std::min(best_cellUpdateTime_ns, (1e+9)*latestFieldTime/numberOfCells);
        if(latestParticleUpdates > 0) best_particleUpdateTime_ns = std::min(best_particleUpdateTime_ns, (1e+9)*latestEnsembleTime/latestParticleUpdates);
        totalEnsembleTime += latestEnsembleTime;
        totalFieldTime += latestFieldTime;
        totalParticleUpdates += latestParticleUpdates;
        totalCellUpdates += numberOfCells;
        if(reportPerformance){
            string statement = "iteration completed in " + to_string(latestEnsembleTime + latestFieldTime) + " s";
            if(latestParticleUpdates > 0) statement += + " (particles: " + to_string(int(100*latestEnsembleTime/(latestEnsembleTime + latestFieldTime))) + "\%); "
                                                                    + to_string((1e+9)*latestEnsembleTime/latestParticleUpdates) + " ns per particle "
                                                                    + "(av_ppc=" + to_string(latest_av_ppc, 1) + ", av_cmr=" + to_string(latest_av_cmr, 2) + ")";
            statement += "; " + to_string((1e+9)*latestFieldTime/numberOfCells) + " ns per cell.";
            pipic_log.message(statement);
        }
    }
    void performanceSummary(){
        ofstream file;
        file.open ("pipic_performance.txt", ios::trunc);
        file << "Time for ensemble/field update per particle/cell and the total time taken." << endl;
        file << endl;
        file << left << setw(32) << "Action" << setw(16) << "Average, ns" << setw(16) <<  "Best, ns" << setw(16) << "Total, s" << endl;
        file << string(16*5, '-') << endl;
        if(totalParticleUpdates > 0)
            file << setw(32) << "Ensemble" << setw(16) << (1e+9)*totalEnsembleTime/totalParticleUpdates << setw(16) << best_particleUpdateTime_ns << setw(16) << totalEnsembleTime << endl;
        file << setw(32) << "Field" << setw(16) << (1e+9)*totalFieldTime/totalCellUpdates << setw(16) << best_cellUpdateTime_ns << setw(16) << totalFieldTime << endl;
        file << endl;
        for(int ih = 0; ih < int(Handler.size()); ih++){
            if(Handler[ih]->total_numberProcessed > 0) {
                file << setw(32) << Handler[ih]->name;
                file << setw(16) << (1e+9)*Handler[ih]->total_partilesTime_s/Handler[ih]->total_numberProcessed;
                file << setw(16) << Handler[ih]->best_partilesTime_ns << setw(16) << Handler[ih]->total_partilesTime_s << endl;
            }
            if(Handler[ih]->actOnCell){
                file << setw(32) << Handler[ih]->name + "(cell)";
                file << setw(16) << (1e+9)*Handler[ih]->total_cellTime_s/totalCellUpdates;
                file << setw(16) << Handler[ih]->best_cellTime_ns << setw(16) << Handler[ih]->total_cellTime_s << endl;
            }
        }
        file.close();
    }
    ~handlerManager(){
        performanceSummary();
        for(int i = 0; i < int(Handler.size()); i++) delete Handler[i];
    }
    void addCellHandler(string name, string subject, int64_t handler, double* dataDouble, int* dataInt){
        Handler.push_back(new cellHandler(name, subject, handler, dataDouble, dataInt));
    }
    void preLoop(vector<particleType> types, bool &fieldHandlerExists){
        vector<string> particleNames;
        for(int it = 0; it < int(types.size()); it++) particleNames.push_back(types[it].name);
        for(int ih = 0; ih < int(Handler.size()); ih++) Handler[ih]->preLoop(particleNames);
        for(int ih = 0; ih < int(Handler.size()); ih++) if(Handler[ih]->actOnCell) fieldHandlerExists = true;
    }
};

template <typename fieldSolver>
void customFieldLoop_via_getEB(fieldSolver *field, int numberOfIterations, int64_t it2coord, int64_t field2data, int64_t dataDouble = 0, int64_t dataInt = 0){
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
        field->getEB(r, E, B);
        field2data_(it, (double*)&r, (double*)&E, (double*)&B, dataDouble_, dataInt_);
    }
};

#endif
