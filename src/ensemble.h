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
// Description: Structure for handling particle ensemble.

#ifndef ENSEMBLE_H
#define ENSEMBLE_H

#include "primitives.h"
#include "services.h"
#include "fourier_solver.h"

struct cellContainer
{
    vector<particle> P;
    int endShift; // the number of particles from the end of the vector that have been already processed (came from other cells) during the ongoing OMP loop
    cellContainer(): endShift(0) {}
    void shuffle(){
        if(unlikely(!randGen)) {cout << "pi-PIC error: rndGen call before assignment" << endl; exit(0);}
        if(P.size() - 1 - endShift >= 1) std::shuffle(begin(P), next(begin(P), P.size() - endShift), *randGen);
    }
};

struct loopLayoutStride8x // structure for making a loop with stride 8 along x
{
    vector<int> offset; // i-th stage is to process cells with offset offset[i] and stride 8
    vector<int> index; // one-cell extended lookup table to determine whether a cell has been processed or not: cells with k-th offset are processed at stage index[k+1];
    int stage; // to be varied from 0 to 7;
    loopLayoutStride8x():
    offset(8), index(10)
    {}
    void makePlan(){
        if(unlikely(!randGen)) {cout << "pi-PIC error: rndGen call before assignment" << endl; exit(0);}
        for(int i = 0; i < int(offset.size()); i++) offset[i] = i;
        shuffle(begin(offset), end(offset), *randGen);
        for(int i = 0; i < int(offset.size()); i++) index[offset[i] + 1] = i;
        index[9] = index[1]; index[0] = index[8];
    }
    inline bool unprocessed(int sx){ // returns true if the cell shifted by 'sx' has not been processed
        return index[offset[stage] + sx + 1] > stage;
    }
};

struct P_pointer{intg ig; int it, ip;};

struct threadData
{
    unsigned long long int numMigrated, numDeleted, numCreated; // counters
    vector<int> toRemove; // list of indices of particles to be removed from the list being processed;
    vector<bool> toRemoveLocal; // indicates whether the particle can be removed immediately or should be relocated after the OMP loop
    vector<P_pointer> postOmpMigrationList; // (ig, it, ip), list of particles that were called to be relocated to a cell that can be potentially operated by another thread at the time
    cellInterface* CI; // pointer to cell interface
    vector<particle> NP; // buffer for new particles
    threadData(): NP(128) {
        toRemove.reserve(32);
        toRemoveLocal.reserve(32);
        postOmpMigrationList.reserve(32);
    }
    void reset() {
        numMigrated = 0; numDeleted = 0; numCreated = 0;
        postOmpMigrationList.clear();
    }
};

struct ensemble
{
    simulationBox box;
    vector<particleType> type; // vector of types
    cellContainer ***cell; // cell[ig][it] is a pointer to cellContainer for type $it$ for $ig$-th cell (nullptr if empty)
    loopLayoutStride8x layout;
    vector<threadData> thread;
    unsigned long long int totalNumberOfParticles;
    bool fieldHandlerExists; // true if there is at least one field handler
    bool shuffle;
    handlerManager Manager;
    rndGen RndGen;
    bool advanceWithOmp;

    ensemble(simulationBox box, int stride = 4): box(box), thread(omp_get_max_threads()),
    fieldHandlerExists(false), shuffle(true), Manager(box.ng), RndGen(box.n.x)
    {
        cell = new cellContainer**[box.ng];
        for(intg ig = 0; ig < box.ng; ig++) cell[ig] = nullptr;
        totalNumberOfParticles = 0;

        for(int iTh = 0; iTh < int(thread.size()); iTh++){
            int *I = new int[15];
            I[3] = box.n.x; I[4] = box.n.y; I[5] = box.n.z;
            I[6] = box.dim; I[7] = sizeof(particle)/8 - 8;
            I[10] = thread[iTh].NP.size(); // capacity of the buffer (is extended automatically when particleSubsetSize exceeds its half on the previous call)
            double *D = new double[16];
            D[0] = box.min.x; D[1] = box.min.y; D[2] = box.min.z;
            D[3] = box.max.x; D[4] = box.max.y; D[5] = box.max.z;
            D[6] = box.step.x; D[7] = box.step.y; D[8] = box.step.z;
            D[9] = box.invStep.x; D[10] = box.invStep.y; D[11] = box.invStep.z;
            double* NP = (double*)(&(thread[iTh].NP[0])); // pointer to the buffer
            thread[iTh].CI = new cellInterface(I, D, nullptr, nullptr, NP);
        }
        RndGen.setTreadLocalRng(0);
    }
    ~ensemble()
    {
        for(intg ig = 0; ig < box.ng; ig++) if(cell[ig] != nullptr)
        {
            for(int it = 0; it < int(type.size()); it++) if(cell[ig][it] != nullptr) delete cell[ig][it];
            delete []cell[ig];
        }
        delete []cell;
        for(int iTh = 0; iTh < int(thread.size()); iTh++){
            delete [](thread[iTh].CI->I);
            delete [](thread[iTh].CI->D);
            if(thread[iTh].CI->F_data != nullptr) delete [](thread[iTh].CI->F_data);
            delete thread[iTh].CI;
        }
    }
    void checkInitCell(intg ig, int typeIndex) // checks if the data of a ig-th cell had been allocated and allocates it if not
    {
        if(cell[ig] == nullptr){
            cell[ig] = new cellContainer*[type.size()];
            for(int it = 0; it < int(type.size()); it++) cell[ig][it] = nullptr;
        }
        if(cell[ig][typeIndex] == nullptr){
            cell[ig][typeIndex] = new cellContainer;
        }
    }
    void checkPushBack(intg ig, int typeIndex, particle &P){
        if(unlikely(cell[ig] == nullptr)){
            cell[ig] = new cellContainer*[type.size()];
            for(int it = 0; it < int(type.size()); it++) cell[ig][it] = nullptr;
        }
        if(unlikely(cell[ig][typeIndex] == nullptr)){
            cell[ig][typeIndex] = new cellContainer;
        }
        cell[ig][typeIndex]->P.push_back(P);
    }
    int placeParticles(string typeName, int totalNumber, double typeCharge, double typeMass, double temperature, int64_t density, int64_t dataDouble = 0, int64_t dataInt = 0) // returns typeIndex
    {
        int typeIndex = -1; // index of particle type; initial code "-1" is to be changed if typeName is found in typeName, otherwise new type is added below
        for(int it = 0; it < int(type.size()); it++)
        if(typeName == type[it].name){
            if((typeCharge == type[it].charge)&&(typeMass == type[it].mass))
                typeIndex = it; // found among previously defined types
            else {
                pipic_log.message("ensemble.addParticles() error: same type name has been declared with different charge and/or mass.", true);
                return -1;
            }
        }
        bool newType = (typeIndex == -1); // newType = true if typeName is not found among previously defined types
        if(newType){
            type.push_back(particleType(typeName, typeCharge, typeMass));
            typeIndex = type.size() - 1;
            //extend the vector of types for all non-empty cells
            for(int iz = 0; iz < box.n.z; iz++)
            for(int iy = 0; iy < box.n.y; iy++)
            for(int ix = 0; ix < box.n.x; ix++)
            if(cell[box.ig({ix, iy, iz})] != nullptr)
            {
                cellContainer **newPointers = new cellContainer*[type.size()];
                for(int it = 0; it < int(type.size()) - 1; it++) newPointers[it] = cell[box.ig({ix, iy, iz})][it];
                newPointers[type.size() - 1] = nullptr;
                cellContainer **tmp = cell[box.ig({ix, iy, iz})];
                cell[box.ig({ix, iy, iz})] = newPointers;
                delete []tmp;
            }
        }

        double(*func_density)(double*, double*, int*) = (double(*)(double*, double*, int*))density;
        double* dataDouble_ = nullptr; if(dataDouble != 0) dataDouble_ = (double*)dataDouble;
        int* dataInt_ = nullptr; if(dataInt != 0) dataInt_ = (int*)dataInt;

        double *Density = new double[box.ng];
        vector<double> totalRealParticles_(omp_get_max_threads(), 0); // counter for each thread (estimate)

        #pragma omp parallel for collapse(3)
        for(int iz = 0; iz < box.n.z; iz++)
        for(int iy = 0; iy < box.n.y; iy++)
        for(int ix = 0; ix < box.n.x; ix++)
        {
            double R[3];
            R[0] = box.min.x + (ix + 0.5)*box.step.x;
            R[1] = box.min.y + (iy + 0.5)*box.step.y;
            R[2] = box.min.z + (iz + 0.5)*box.step.z;
            Density[box.ig({ix, iy, iz})] = func_density(R, dataDouble_, dataInt_);
            totalRealParticles_[omp_get_thread_num()] += Density[box.ig({ix, iy, iz})]*box.step.x*box.step.y*box.step.z;
        }

        double totalRealParticles = 0; // total number of real particles
        for(int i = 0; i < int(totalRealParticles_.size()); i++) totalRealParticles += totalRealParticles_[i];
        double weight = totalRealParticles/totalNumber; // estimated weight to be used for all particles

        // estimating and reserving memory for particle allocation
        for(int iz = 0; iz < box.n.z; iz++)
        for(int iy = 0; iy < box.n.y; iy++)
        for(int ix = 0; ix < box.n.x; ix++)
        {
            intg ig = box.ig({ix, iy, iz});
            int numberToReserve = int(1.5*Density[ig]*box.step.x*box.step.y*box.step.z/weight); // the center of cells are shifted by half step relative to node location, but here we use it as an estimate
            if(numberToReserve > 0){
                checkInitCell(ig, typeIndex);
                cell[ig][typeIndex]->P.reserve(cell[ig][typeIndex]->P.size() + numberToReserve);
            }
        }

        //adding particles
        for(int iz = 0; iz < box.n.z; iz++)
        for(int iy = 0; iy < box.n.y; iy++)
        for(int ix = 0; ix < box.n.x; ix++)
        {
            intg ig = box.ig({ix, iy, iz});
            double expectedNumber = Density[ig]*box.step.x*box.step.y*box.step.z/weight;
            int numberToGenerate = int(expectedNumber) + (rand_double(0, 1) < (expectedNumber - int(expectedNumber)));
            if(numberToGenerate > 0)
            for(int ip = 0; ip < numberToGenerate; ip++){
                particle P;
                P.r.x = box.min.x + (ix + rand_double(0, 1))*box.step.x;
                P.r.y = box.min.y + (iy + rand_double(0, 1))*box.step.y;
                P.r.z = box.min.z + (iz + rand_double(0, 1))*box.step.z;
                P.p = generateMomentum(typeMass, temperature);
                P.w = weight;
                P.id = generateID();
                checkPushBack(ig, typeIndex, P);
                totalNumberOfParticles++;
            }
        }
        delete []Density;
        return typeIndex;
    }
    bool checkLocation(double3 r, double3 min, double3 max){
        if(r.x < min.x) return false;
        if(r.y < min.y) return false;
        if(r.z < min.z) return false;
        if(r.x >= max.x) return false;
        if(r.y >= max.y) return false;
        if(r.z >= max.z) return false;
        return true;
    }
    void accommodateNPfromCI(unsigned int ig, string moduleName, bool directOrder){
        threadData &activeThread(thread[omp_get_thread_num()]);
        for(int j = 0; j < activeThread.CI->particleBufferSize; j++){
            int it = activeThread.NP[j].id; // by convention the type of new particles is placed to id;
            activeThread.NP[j].id = generateID();
            if(checkLocation(activeThread.NP[j].r, activeThread.CI->cellMin(), activeThread.CI->cellMax())){
                checkPushBack(ig, it, activeThread.NP[j]);
                activeThread.numCreated++;
                if(directOrder) {
                    if(it >= activeThread.CI->particleTypeIndex) cell[ig][it]->endShift++;
                } else {
                    if(it <= activeThread.CI->particleTypeIndex) cell[ig][it]->endShift++;
                }
            } else
                pipic_log.message("pi-PIC error: ignoring an attempt of module<" + moduleName + "> to add a particle outside current cell.", true);
        }
        if(unlikely(2*activeThread.CI->particleBufferSize > activeThread.CI->particleBufferCapacity)){
            activeThread.NP.resize(2*activeThread.NP.size());
            activeThread.CI->I[10] = activeThread.NP.size();
            activeThread.CI->NP_data = (double*)(&(activeThread.NP[0]));
        }
        activeThread.CI->particleBufferSize = 0;
    }
    template<typename pic_solver, typename field_solver>
    void apply_actOnCellHandlers(pic_solver *Solver, threadData &activeThread, bool &fieldBeenSet, intg ig, bool directOrder){
        fieldBeenSet = false;
        for(int ih = 0; ih < int(Manager.Handler.size()); ih++)
        if(Manager.Handler[ih]->actOnCell){
            if(!fieldBeenSet) ((field_solver*)(Solver->Field))->cellSetField(*activeThread.CI, activeThread.CI->i);
            fieldBeenSet = true;
            Manager.Handler[ih]->handle(activeThread.CI);
            accommodateNPfromCI(ig, Manager.Handler[ih]->name, directOrder);
        }
    }
    template<typename pic_solver, typename field_solver>
    void apply_particleHandlers(int it, pic_solver *Solver, threadData &activeThread, bool &fieldBeenSet, intg ig, bool directOrder){
        for(int ih = 0; ih < int(Manager.Handler.size()); ih++)
        if(Manager.Handler[ih]->actOn[it]){
            if(!fieldBeenSet) ((field_solver*)(Solver->Field))->cellSetField(*activeThread.CI, activeThread.CI->i);
            fieldBeenSet = true;
            Manager.Handler[ih]->handle(activeThread.CI);
            accommodateNPfromCI(ig, Manager.Handler[ih]->name, directOrder);
            activeThread.CI->P_data = (double*)(&(cell[ig][it]->P[0]));
            activeThread.CI->I[9] = cell[ig][it]->P.size() - cell[ig][it]->endShift; // set particleSubsetSize;
            if(activeThread.CI->particleSubsetSize == 0) break;
        }
    }
    void processCharglessParticle(double mass, double timeStep, particle &P){
        if(mass == 0){ // photon
            double pNorm = P.p.norm();
            if(likely(pNorm > 0)) P.r += (timeStep*lightVelocity/pNorm)*P.p;
        } else
            P.r += (timeStep*lightVelocity/sqrt(sqr(mass*lightVelocity) + P.p.norm()))*P.p;
    }
    template<typename pic_solver, typename field_solver>
    void advance_singleLoop(pic_solver *Solver, double timeStep){
        chronometer chronometerCells;
        chronometerCells.start();
        Solver->preLoop();
        chronometerCells.stop();

        layout.makePlan();
        for(int iTh = 0; iTh < int(thread.size()); iTh++){
            thread[iTh].reset();
            thread[iTh].CI->D[12] = timeStep;
        }
        Manager.latestParticleUpdates = totalNumberOfParticles;
        Manager.preLoop(type, fieldHandlerExists);
        chronometer chronometerEnsemble;
        chronometerEnsemble.start();
        for(layout.stage = 0; layout.stage < 8; layout.stage++)
        {
            #pragma omp parallel for collapse(1) if(advanceWithOmp)
            for(int ix = layout.offset[layout.stage]; ix < box.n.x; ix += 8)
            for(int iz = 0; iz < box.n.z; iz += 1) // could be good to shuffle here, but unprocessed() then needs a modification
            for(int iy = 0; iy < box.n.y; iy += 1) // could be good to shuffle here
            {
                RndGen.setTreadLocalRng(ix);
                intg ig = box.ig({ix, iy, iz});
                if((cell[ig] != nullptr)||(fieldHandlerExists)){
                    threadData &activeThread(thread[omp_get_thread_num()]);
                    //download type-independent data to cell interface:
                    activeThread.CI->I[0] = ix; activeThread.CI->I[1] = iy; activeThread.CI->I[2] = iz;
                    activeThread.CI->I[11] = 0;
                    activeThread.CI->I[8] = -1; // indicates "act on cell"
                    activeThread.CI->I[9] = 0;
                    activeThread.CI->I[10] = activeThread.NP.size();
                    activeThread.CI->I[13] = omp_get_thread_num();
                    activeThread.CI->I[14] = rand_int(std::numeric_limits<int>::min(), std::numeric_limits<int>::max());
                    bool fieldBeenSet = false;
                    apply_actOnCellHandlers<pic_solver, field_solver>(Solver, activeThread, fieldBeenSet, ig, true);
                    if(cell[ig] != nullptr)
                    for(int it = 0; it < int(type.size()); it++)
                    if(cell[ig][it] != nullptr){
                        if(cell[ig][it]->P.size() - cell[ig][it]->endShift > 0){
                            activeThread.CI->I[8] = it;
                            activeThread.CI->I[9] = cell[ig][it]->P.size() - cell[ig][it]->endShift;
                            activeThread.CI->D[13] = type[it].charge;
                            activeThread.CI->D[14] = type[it].mass;
                            activeThread.CI->P_data = (double*)(&(cell[ig][it]->P[0]));
                            Solver->startSubLoop(activeThread.CI->i, type[it].charge, type[it].mass, timeStep);
                            if(shuffle) cell[ig][it]->shuffle();
                            apply_particleHandlers<pic_solver, field_solver>(it, Solver, activeThread, fieldBeenSet, ig, true);
                            activeThread.toRemove.clear();
                            activeThread.toRemoveLocal.clear();
                            for(int ip = 0; ip < int(cell[ig][it]->P.size()) - cell[ig][it]->endShift; ip++){
                                if(likely(cell[ig][it]->P[ip].w != 0)){
                                    if(likely(type[it].charge != 0)) Solver->processParticle(cell[ig][it]->P[ip], type[it].charge, type[it].mass, timeStep);
                                    else processCharglessParticle(type[it].mass, timeStep, cell[ig][it]->P[ip]);
                                    bool move, postOmpMove;
                                    checkMove(&(cell[ig][it]->P[ip]), move, postOmpMove);
                                    if(unlikely(move)){
                                        if(likely(!postOmpMove)){
                                            activeThread.toRemove.push_back(ip);
                                            activeThread.toRemoveLocal.push_back(true);
                                            activeThread.numMigrated++; // neighbor cell migration
                                        } else {
                                            activeThread.toRemove.push_back(ip);
                                            activeThread.toRemoveLocal.push_back(false);
                                        }
                                    }
                                } else {
                                    activeThread.toRemove.push_back(ip);
                                    activeThread.toRemoveLocal.push_back(true);
                                    activeThread.numDeleted++;
                                }
                            }
                            compresList(cell[ig][it]->P, activeThread, ig, it);
                            Solver->endSubLoop();
                        }
                        cell[ig][it]->endShift = 0;
                    }
                }
            }
            RndGen.setTreadLocalRng(0);
        }
        int overCellRelocated = 0;
        for(int iTh = 0; iTh < int(thread.size()); iTh++)
        {
            for(int il = thread[iTh].postOmpMigrationList.size() - 1; il >= 0; il--)
            {
                int ig =  thread[iTh].postOmpMigrationList[il].ig;
                int it =  thread[iTh].postOmpMigrationList[il].it;
                int ip =  thread[iTh].postOmpMigrationList[il].ip;
                addParticle_general(cell[ig][it]->P[ip], it, true);
                if(ip < int(cell[ig][it]->P.size()) - 1) memcpy(&(cell[ig][it]->P[ip]), &(cell[ig][it]->P.back()), sizeof(particle));
                cell[ig][it]->P.pop_back();
            }
            overCellRelocated += thread[iTh].postOmpMigrationList.size();
            thread[iTh].postOmpMigrationList.clear();
            totalNumberOfParticles += thread[iTh].numCreated - thread[iTh].numDeleted;
        }
        chronometerEnsemble.stop();
        Manager.latestEnsembleTime = chronometerEnsemble.getTime_s();
        if(overCellRelocated > 0) pipic_log.message("pi-PIC warning: " + to_string(overCellRelocated) + " overcell migrations; consider reducing time step");

        Manager.latest_av_ppc = double(totalNumberOfParticles)/double(box.ng);
        unsigned long long int migrationCounter = 0; for(int iTh = 0; iTh < int(thread.size()); iTh++)migrationCounter += thread[iTh].numMigrated;
        Manager.latest_av_cmr = migrationCounter/double(totalNumberOfParticles);
        chronometerCells.start();
        Solver->postLoop();
        Solver->Field->advance(timeStep);
        chronometerCells.stop();
        Manager.latestFieldTime = chronometerCells.getTime_s();
    }
    template<typename pic_solver, typename field_solver>
    void advance_doubleLoop(pic_solver *Solver, double timeStep)
    {
        chronometer chronometerCells;
        chronometerCells.start();
        Solver->preLoop(0);
        chronometerCells.stop();
        layout.makePlan();
        for(int iTh = 0; iTh < int(thread.size()); iTh++){
            thread[iTh].reset();
            thread[iTh].CI->D[12] = timeStep;
        }
        Manager.latestParticleUpdates = totalNumberOfParticles;
        Manager.preLoop(type, fieldHandlerExists);
        chronometer chronometerEnsemble;
        chronometerEnsemble.start();
        //first loop:
        for(layout.stage = 7; layout.stage >= 0; layout.stage--)
        {
            #pragma omp parallel for collapse(1) if(advanceWithOmp)
            for(int ix = box.n.x - 8 + layout.offset[layout.stage]; ix >= 0; ix -= 8)
            for(int iz = box.n.z - 1; iz >= 0; iz -= 1)
            for(int iy = box.n.y - 1; iy >= 0; iy -= 1)
            {
                RndGen.setTreadLocalRng(ix);
                intg ig = box.ig({ix, iy, iz});
                if((cell[ig] != nullptr)||(fieldHandlerExists)){
                    threadData &activeThread(thread[omp_get_thread_num()]);
                    //download type-independent data to cell interface:
                    activeThread.CI->I[0] = ix; activeThread.CI->I[1] = iy; activeThread.CI->I[2] = iz;
                    activeThread.CI->I[11] = 0;
                    activeThread.CI->I[8] = -1; // indicates "act on cell"
                    activeThread.CI->I[9] = 0;
                    activeThread.CI->I[10] = activeThread.NP.size();
                    activeThread.CI->I[13] = omp_get_thread_num();
                    activeThread.CI->I[14] = rand_int(std::numeric_limits<int>::min(), std::numeric_limits<int>::max());
                    bool fieldBeenSet = false;
                    apply_actOnCellHandlers<pic_solver, field_solver>(Solver, activeThread, fieldBeenSet, ig, true);
                    if(cell[ig] != nullptr)
                    for(int it = type.size() - 1; it >= 0; it--)
                    if(cell[ig][it] != nullptr)
                    if(cell[ig][it]->P.size() - cell[ig][it]->endShift > 0){
                        activeThread.CI->I[8] = it;
                        activeThread.CI->I[9] = cell[ig][it]->P.size() - cell[ig][it]->endShift;
                        activeThread.CI->D[13] = type[it].charge;
                        activeThread.CI->D[14] = type[it].mass;
                        activeThread.CI->P_data = (double*)(&(cell[ig][it]->P[0]));
                        Solver->startSubLoop(activeThread.CI->i, type[it].charge, type[it].mass, timeStep, 0);
                        if(shuffle) cell[ig][it]->shuffle();
                        apply_particleHandlers<pic_solver, field_solver>(it, Solver, activeThread, fieldBeenSet, ig, true);
                        for(int ip = cell[ig][it]->P.size() - cell[ig][it]->endShift - 1; ip >= 0; ip--)
                            if(likely(cell[ig][it]->P[ip].w != 0)){
                                if(likely(type[it].charge != 0)) Solver->processParticle(cell[ig][it]->P[ip], type[it].charge, type[it].mass, timeStep, 0);
                                else processCharglessParticle(type[it].mass, timeStep, cell[ig][it]->P[ip]); // advance for the whole step here, do nothing on the second loop
                            }
                        Solver->endSubLoop(0);
                    }
                }
            }
            RndGen.setTreadLocalRng(0);
        }
        chronometerEnsemble.stop();
        chronometerCells.start();
        Solver->postLoop(0);
        Solver->preLoop(1);
        chronometerCells.stop();
        chronometerEnsemble.start();
        //second loop:
        for(layout.stage = 0; layout.stage < 8; layout.stage++)
        {
            #pragma omp parallel for collapse(1) if(advanceWithOmp)
            for(int ix = layout.offset[layout.stage]; ix < box.n.x; ix += 8)
            for(int iz = 0; iz < box.n.z; iz += 1)
            for(int iy = 0; iy < box.n.y; iy += 1)
            {
                RndGen.setTreadLocalRng(ix);
                intg ig = box.ig({ix, iy, iz});
                if(cell[ig] != nullptr){
                    threadData &activeThread(thread[omp_get_thread_num()]);
                    //download type-independent data to cell interface:
                    activeThread.CI->I[0] = ix; activeThread.CI->I[1] = iy; activeThread.CI->I[2] = iz;
                    activeThread.CI->I[11] = 0;

                    if(cell[ig] != nullptr)
                    for(int it = 0; it < int(type.size()); it++)
                    if(cell[ig][it] != nullptr){
                        if(cell[ig][it]->P.size() - cell[ig][it]->endShift > 0){
                            activeThread.CI->I[8] = it;
                            activeThread.CI->I[9] = cell[ig][it]->P.size() - cell[ig][it]->endShift;
                            activeThread.CI->D[13] = type[it].charge;
                            activeThread.CI->D[14] = type[it].mass;
                            activeThread.CI->P_data = (double*)(&(cell[ig][it]->P[0]));

                            activeThread.toRemove.clear();
                            activeThread.toRemoveLocal.clear();
                            Solver->startSubLoop(activeThread.CI->i, type[it].charge, type[it].mass, timeStep, 1);
                            for(int ip = 0; ip < int(cell[ig][it]->P.size()) - cell[ig][it]->endShift; ip++){
                                if(likely(cell[ig][it]->P[ip].w != 0))
                                {
                                    if(likely(type[it].charge != 0)) Solver->processParticle(cell[ig][it]->P[ip], type[it].charge, type[it].mass, timeStep, 1); // chargeless particles are advanced for the whole step during the first loop
                                    bool move, postOmpMove;
                                    checkMove(&(cell[ig][it]->P[ip]), move, postOmpMove);
                                    if(unlikely(move)){
                                        if(likely(!postOmpMove)){
                                            activeThread.toRemove.push_back(ip);
                                            activeThread.toRemoveLocal.push_back(true);
                                            activeThread.numMigrated++; // neighbor cell migration
                                        } else {
                                            activeThread.toRemove.push_back(ip);
                                            activeThread.toRemoveLocal.push_back(false);
                                        }
                                    }
                                } else {
                                    activeThread.toRemove.push_back(ip);
                                    activeThread.toRemoveLocal.push_back(true);
                                    activeThread.numDeleted++;
                                }
                            }
                            compresList(cell[ig][it]->P, activeThread, ig, it);
                            Solver->endSubLoop(1);
                        }
                        cell[ig][it]->endShift = 0;
                    }
                }
            }
            RndGen.setTreadLocalRng(0);
        }
        int overCellRelocated = 0;
        for(int iTh = 0; iTh < int(thread.size()); iTh++)
        {
            for(int il = thread[iTh].postOmpMigrationList.size() - 1; il >= 0; il--)
            {
                int ig =  thread[iTh].postOmpMigrationList[il].ig;
                int it =  thread[iTh].postOmpMigrationList[il].it;
                int ip =  thread[iTh].postOmpMigrationList[il].ip;
                addParticle_general(cell[ig][it]->P[ip], it, true);
                if(ip < int(cell[ig][it]->P.size()) - 1) memcpy(&(cell[ig][it]->P[ip]), &(cell[ig][it]->P.back()), sizeof(particle));
                cell[ig][it]->P.pop_back();
            }
            overCellRelocated += thread[iTh].postOmpMigrationList.size();
            thread[iTh].postOmpMigrationList.clear();
            totalNumberOfParticles += thread[iTh].numCreated - thread[iTh].numDeleted;
        }
        chronometerEnsemble.stop();
        Manager.latestEnsembleTime = chronometerEnsemble.getTime_s();

        if(overCellRelocated > 0) pipic_log.message("pi-PIC warning: " + to_string(overCellRelocated) + " overcell migrations; consider reducing time step.");

        Manager.latest_av_ppc = double(totalNumberOfParticles)/double(box.ng);
        unsigned long long int migrationCounter = 0; for(int iTh = 0; iTh < int(thread.size()); iTh++)migrationCounter += thread[iTh].numMigrated;
        Manager.latest_av_cmr = migrationCounter/double(totalNumberOfParticles);

        chronometerCells.start();
        Solver->postLoop(1);
        Solver->Field->advance(timeStep);
        chronometerCells.stop();
        Manager.latestFieldTime = chronometerCells.getTime_s();
    }
    inline void compresList(vector<particle> &P, threadData &thread, intg ig, int it)
    {
        int shift = 0;
        for(int k = thread.toRemove.size() - 1; k >= 0; k--){
            if(!thread.toRemoveLocal[k]){
                if(thread.toRemove[k] < int(P.size()) - 1 - shift) swap(P[thread.toRemove[k]], P[P.size() - 1 - shift]);
                shift++;
            } else {
                if(thread.toRemove[k] < int(P.size()) - 1 - shift) memcpy(&P[thread.toRemove[k]], &P[P.size() - 1 - shift], sizeof(particle));
                if(shift > 0) memcpy(&P[P.size() - 1 - shift], &P[P.size() - 1], sizeof(particle));
                P.pop_back();
            }
        }
        for(int ip = int(P.size()) - shift; ip < int(P.size()); ip++) thread.postOmpMigrationList.push_back({ig, it, ip});
    }
    inline void cyclicShift(int &i, int s, int n, bool &limCross){ // optimized for likely cases of abs(s) < n
        i += s;
        if(unlikely(i >= n)) {i -= n; limCross = true;}
            else if(unlikely(i < 0)) {i += n; limCross = true;}
        if(unlikely(i >= n)) {i = i%n; limCross = true;}
            else if(unlikely(i < 0)) {i = n - 1 - (-i - 1)%n; limCross = true;}
    }
    inline void checkMove(particle *P, bool &move, bool &postOmpMove){
        threadData &activeThread(thread[omp_get_thread_num()]);
        int ix = activeThread.CI->i.x, iy = activeThread.CI->i.y, iz = activeThread.CI->i.z, it = activeThread.CI->particleTypeIndex;

        int sx = floor((P->r.x - box.min.x)*box.invStep.x) - ix;
        int sy = 0; if(box.n.y > 1) sy = floor((P->r.y - box.min.y)*box.invStep.y) - iy;
        int sz = 0; if(box.n.z > 1) sz = floor((P->r.z - box.min.z)*box.invStep.z) - iz;

        if(unlikely(P->r.x > box.max.x)) P->r.x -= box.size.x; else if(unlikely(P->r.x < box.min.x)) P->r.x += box.size.x;
        if(unlikely(P->r.y > box.max.y)) P->r.y -= box.size.y; else if(unlikely(P->r.y < box.min.y)) P->r.y += box.size.y;
        if(unlikely(P->r.z > box.max.z)) P->r.z -= box.size.z; else if(unlikely(P->r.z < box.min.z)) P->r.z += box.size.z;

        if(unlikely((sx != 0)||(sy != 0)||(sz != 0))){
            move = true;
            if(likely((abs(sx) <= 1)&&(abs(sy) <= 1)&&(abs(sz) <= 1))){
                int3 newI = {ix, iy, iz};
                bool limCorssed = false;
                cyclicShift(newI.x, sx, box.n.x, limCorssed);
                cyclicShift(newI.y, sy, box.n.y, limCorssed);
                cyclicShift(newI.z, sz, box.n.z, limCorssed);
                {
                    intg newIg = box.ig(newI);
                    checkPushBack(newIg, it, *P);
                    if((layout.unprocessed(sx))||((sx == 0)&&(newI.y + box.n.y*newI.z > iy + box.n.y*iz))) cell[newIg][it]->endShift++;
                    postOmpMove = false;
                }
            } else postOmpMove = true;
        } else move = false;
    }
    void addParticle_general(particle &P, int it, bool migrated = false) //nonOMP
    {
        if(unlikely((P.r.x < box.min.x)||(P.r.x >= box.max.x)||(P.r.y < box.min.y)||(P.r.y >= box.max.y)||(P.r.z < box.min.z)||(P.r.z >= box.max.z))){
            pipic_log.message("pi-PIC error: addParticle_general() ingnors an attempt to add a particle that is outside computation region.", true);
            return;
        }
        int ix = floor((P.r.x - box.min.x)*box.invStep.x);
        int iy = floor((P.r.y - box.min.y)*box.invStep.y);
        int iz = floor((P.r.z - box.min.z)*box.invStep.z);
        intg ig = box.ig({ix, iy, iz});
        checkPushBack(ig, it, P);
        if(!migrated) totalNumberOfParticles++;
    }
    struct nonOmpIterator //forward iterator for use without OMP
    {
        nonOmpIterator(int it, ensemble *Ensemble) : ig(0), ip(-1), it(it), Ensemble(Ensemble) { (*this)++; }
        particle& operator*() const { return *Particle; }
        void removeZeroWeightParticles(vector<particle> &P){ // removes particles in P assigned to be removed (w = 0)
            for(int ip = P.size() - 1; ip >= 0; ip--)
            if(P[ip].w == 0) {
                if(ip < int(P.size()) - 1) memcpy(&(P[ip]), &(P[P.size() - 1]), sizeof(particle));
                P.pop_back();
                Ensemble->totalNumberOfParticles--;
            }
        }
        friend bool operator < (nonOmpIterator& lhs, int const& rhs){
            if(lhs.Particle != nullptr)
                if(lhs.ip == int(lhs.Ensemble->cell[lhs.ig][lhs.it]->P.size()) - 1)
                    lhs.removeZeroWeightParticles(lhs.Ensemble->cell[lhs.ig][lhs.it]->P);
            return (lhs.Particle != nullptr);
        }
        bool operator++(int)
        {
            Particle = nullptr;
            while(Particle == nullptr)
            {
                if(Ensemble->cell[ig] == nullptr) {
                    if(ig < Ensemble->box.ng - 1) ig++;
                        else return false;
                } else {
                    if(Ensemble->cell[ig][it] == nullptr){
                        if(ig < Ensemble->box.ng - 1) ig++;
                        else return false;
                    } else {
                        if(ip + 1 < int(Ensemble->cell[ig][it]->P.size())){
                            ip++;
                            Particle = &(Ensemble->cell[ig][it]->P[ip]);
                            return true;
                        } else {
                            removeZeroWeightParticles(Ensemble->cell[ig][it]->P);
                            ip = -1;
                            if(ig < Ensemble->box.ng - 1) ig++;
                            else return false;
                        }
                    }
                }
            }
            return false;
        }
        intg ig;
        int ip;
        int it;
        particle *Particle;
        ensemble *Ensemble;
    };
    nonOmpIterator begin(int typeIndex) {return nonOmpIterator(typeIndex, this);}
    int end() {return 0;} // this is just to match the standard interface

    int getTypeIndex(string typeName)
    {
        int typeIndex = -1; // error code value
        for(int it = 0; it < int(type.size()); it++) if(typeName == type[it].name) typeIndex = it;
        if(typeIndex == -1){
            pipic_log.message("pipic error: unknown particle type name '" + typeName + "'.", true);
            exit(-1);
        }
        return typeIndex;
    }
};

#endif
