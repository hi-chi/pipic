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
// Description of the file: File introduces structures for storing and making loops over particles.

#ifndef ENSEMBLE_H
#define ENSEMBLE_H

#include "particle.h"
#include "cic_weighting.h"

struct cellContainer // container for particles of the same type
{
    vector<particle> P;
    unsigned int endShift; // parameter indicating how many particles from the end to exclude from the loop (newcoming particles that have already been processed in the current loop) 
    unsigned int doubleLoopMidPoint; // auxilary variable for shuffling: we first process particles after MidPoint and then before in a separate loop 
    cellContainer(): endShift(0), doubleLoopMidPoint(0) {}
    void shuffle()
    {
        static thread_local mt19937* randGen = nullptr;
        if (!randGen) randGen = new mt19937(clock() + omp_get_thread_num());
        if(P.size() - 1 - endShift >= 1) std::shuffle(std::begin(P), std::next(std::begin(P), P.size() - endShift), *randGen);
    }
    int begin(double &cellShuffleProbability, bool &doubleLoop, int &doubleLoopStage) // starts a loop and returns the index of the particle
    {
        if((!doubleLoop)||(doubleLoopStage == 0)) 
            if(cellShuffleProbability > 0) 
                if(rand_double() < cellShuffleProbability) shuffle();
        
        if(!doubleLoop) return 0;
        else {
            if(doubleLoopStage == 0){ // at stage 0 we start from the middle point 
                doubleLoopMidPoint = (P.size() - endShift)/2 + ((P.size() - endShift)%2)*(rand_double() < 0.5);
                return doubleLoopMidPoint;
            } else{ // doubleLoopStage == 1
                endShift = P.size() - doubleLoopMidPoint;
                return 0;
            }
        } 
    }
    bool end(int ip) // return 'true' when ip reaches the size of container
    {
        if(ip < P.size() - endShift) {
            return false;
        } else {
            endShift = 0;
            return true;
        }
    }
    particle& operator[](int ip)
    {
        return P[ip];
    }
    void addParticle(particle &Particle, bool unprocessedCell = false) // uprocessedCell = true neans that the particle should be excluded when considered in the ongoing macroloop
    {
        P.push_back(Particle);
        if(unprocessedCell) endShift++;
    }
    void removeParticle(int &ip)//right after calling this function (in a loop) ip should be icreased and assessed by end(ip)
    {
        if(ip < P.size() - 1) memcpy(&P[ip], &P[P.size() - 1], sizeof(particle));
        P.pop_back();
        if(endShift == 0)    
            ip--;
        else 
            endShift--;
    }
    void sendParticle(int &ip, cellContainer &toCell, bool destinationIsUnprocessedCell = false) // unprocessedToCell = true neans that the particle should be excluded when considered in the ongoing macroloop
    //right after calling this function (in a loop) ip should be icreased and assessed by end(ip)
    {
        toCell.addParticle(P[ip], destinationIsUnprocessedCell);
        removeParticle(ip);
    }
};

struct ensemble
{
    simulationBox box;
    vector<particleType> type; // vector of types 
    cellContainer ***data; // data[ig][it] is a pointer to cellContainer for type $it$ for $ig$-th cell (nullptr if empty)
    // cell[box.ig({ix, iy, iz})] contains particles with r.x \in [box.min.x + (ix-0.5)*box.step.x; box.min.x + (ix+0.5)*box.step.x) and similar for r.y and r.z
    // Note that CIC for particles in cell[box.ig({ix, iy, iz})] concerns nodes with ix, ix-1 and similar for y and z;
    particleIdGenerator idGen;
    vector<int> ompLoopIt;
    vector<bool> ompLoop_flagDeleteCurrent;
    int ompLoopOg;
    vector<int3> ompLoopI; 
    int3 ompLoopOffset;
    bool ompLoopIsRunning;
    vector<int> particleCount;
    int totalNumberOfParticles; // warning: this value is not being updated during ompLoop
    vector<vector<int3>> violationList; // (ig, particleType, ip), list of particles that are relocated for more than a cell during ompLoop 
    bool periodicCoordinates; // true means keeping coordinates within the bounds of box
    
    double shuffle_cellShuffleProbability;
    bool shuffle_doubleLoop;
    int shuffle_doubleLoopStage;

    ensemble(simulationBox box): box(box), 
    ompLoopIsRunning(false), ompLoopIt(omp_get_max_threads(), 0), ompLoop_flagDeleteCurrent(omp_get_max_threads(), false), 
    particleCount(omp_get_max_threads(), 0), totalNumberOfParticles(0), violationList(omp_get_max_threads()), periodicCoordinates(true),
    ompLoopI(omp_get_max_threads(), {0, 0, 0})
    {
        data = new cellContainer**[box.ng];
        for(size_t ig = 0; ig < box.ng; ig++) data[ig] = nullptr;
        shuffle_cellShuffleProbability = 1;
        shuffle_doubleLoop = true;
        shuffle_doubleLoopStage = 0;
    }
    ~ensemble()
    {
        for(size_t ig = 0; ig < box.ng; ig++) if(data[ig] != nullptr)
        {
            for(int it = 0; it < type.size(); it++) if(data[ig][it] != nullptr) delete data[ig][it];
            delete []data[ig];
        } 
        delete []data;
    }
    void checkInitCell(size_t ig, int typeIndex) // checks if the cell had been created and creates it if not
    {
        if(data[ig] == nullptr) 
        {
            data[ig] = new cellContainer*[type.size()];
            for(int it = 0; it < type.size(); it++) data[ig][it] = nullptr;
        }
        if(data[ig][typeIndex] == nullptr)
        {
            data[ig][typeIndex] = new cellContainer;
        }
    }
    int placeParticles(string typeName, int totalNumber, double typeCharge, double typeMass, double temperature, int64_t density, int64_t dataDouble = 0, int64_t dataInt = 0) // returns typeIndex
    {
        int typeIndex = -1; // index of particle type; initial code "-1" is to be changed if typeName is found in typeName, otherwise new type is added below
        for(int it = 0; it < type.size(); it++) 
        if(typeName == type[it].name)
        if((typeCharge == type[it].charge)&&(typeMass == type[it].mass))
            typeIndex = it; // found among previously defined types
        else {
            pipic_log.message("ensemble.addParticles() error: same type name must have same charge and mass.", true); 
            return -1;
        }
        bool newType = (typeIndex == -1); // newType = true if typeName is not found among previously defined types
        if(newType) 
        {
            type.push_back(particleType(typeName, typeCharge, typeMass));
            typeIndex = type.size() - 1;
            //extend the vector of types for all non-empty cells
            for(int iz = 0; iz < box.n.z; iz++)
            for(int iy = 0; iy < box.n.y; iy++)
            for(int ix = 0; ix < box.n.x; ix++)
            if(data[box.ig({ix, iy, iz})] != nullptr)
            {
                cellContainer **newPointers = new cellContainer*[type.size()];
                for(int it = 0; it < type.size() - 1; it++) newPointers[it] = data[box.ig({ix, iy, iz})][it];
                newPointers[type.size() - 1] = nullptr;
                cellContainer **tmp = data[box.ig({ix, iy, iz})];
                data[box.ig({ix, iy, iz})] = newPointers;
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
            double3 R3 = box.nodeLocation({ix, iy, iz});
            double R[3]; R[0] = R3.x; R[1] = R3.y; R[2] = R3.z;
            Density[box.ig({ix, iy, iz})] = func_density(R, dataDouble_, dataInt_);
            totalRealParticles_[omp_get_thread_num()] += Density[box.ig({ix, iy, iz})]*box.step.x*box.step.y*box.step.z;
        }
        
        double totalRealParticles = 0; // total nuber of real particles
        for(int i = 0; i < totalRealParticles_.size(); i++) totalRealParticles += totalRealParticles_[i];
        double weight = totalRealParticles/totalNumber; // estimated weight to be used for all particles
        
        // estimating and reserving memory for particle allocation
        for(int iz = 0; iz < box.n.z; iz++)
        for(int iy = 0; iy < box.n.y; iy++)
        for(int ix = 0; ix < box.n.x; ix++)
        {
            size_t ig = box.ig({ix, iy, iz});
            int numberToReserve = int(1.5*Density[ig]*box.step.x*box.step.y*box.step.z/weight); // the center of cells are shifted by half step relative to node location, but here we use it as an estimate  
            if(numberToReserve > 0)
            {
                checkInitCell(ig, typeIndex);
                data[ig][typeIndex]->P.reserve(data[ig][typeIndex]->P.size() + numberToReserve);
            }
        }
        
        //adding particles
        for(int iz = 0; iz < box.n.z; iz++)
        for(int iy = 0; iy < box.n.y; iy++)
        for(int ix = 0; ix < box.n.x; ix++)
        {
            size_t ig = box.ig({ix, iy, iz});
            double expectedNumber = Density[ig]*box.step.x*box.step.y*box.step.z/weight;
            int numberToGenerate = int(expectedNumber) + (rand_double() < (expectedNumber - int(expectedNumber)));
            if(numberToGenerate > 0)
            for(int ip = 0; ip < numberToGenerate; ip++)
            {
                particle P;
                P.r.x = box.min.x + (ix + rand_double())*box.step.x;
                P.r.y = box.min.y + (iy + rand_double())*box.step.y;
                P.r.z = box.min.z + (iz + rand_double())*box.step.z;
                P.p = generateMomentum(typeMass, temperature);
                P.w = weight;
                P.id = idGen.issue();
                addParticle(P, typeIndex);
            }
        }
        delete []Density;
        return typeIndex;
    }
    int upperCell(double coord, double min, double invSize, int n)
    {
        double u = n*(coord - min)*invSize + 0.5;
        return ((int)floor(u) % n + n) % n;
    }
    int3 int3_coord(double3 &r)
    {
        int ix = 0, iy = 0, iz = 0;
        ix = upperCell(r.x, box.min.x, box.invSize.x, box.n.x);
        if(box.n.y != 1) iy = upperCell(r.y, box.min.y, box.invSize.y, box.n.y);
        if(box.n.z != 1) iz = upperCell(r.z, box.min.z, box.invSize.z, box.n.z);
        return {ix, iy, iz};
    }
    size_t ig_coord(double3 &r) // returns global index
    {
        return box.ig(int3_coord(r));
    }
    int cyclicDist(int i1, int i2, int n) // returns distance in cyclic space
    {
        int d = abs(i1 - i2);
        if(n - d < d) return n - d; else return d;
    }
    void addParticle(particle &P, int typeIndex)
    {
        if(!ompLoopIsRunning)
        {
            size_t ig = ig_coord(P.r);
            checkInitCell(ig, typeIndex);
            data[ig][typeIndex]->P.push_back(P);
            totalNumberOfParticles++;
        } else {
            int3 new_i = int3_coord(P.r);
            int ix = ompLoopI[omp_get_thread_num()].x, iy = ompLoopI[omp_get_thread_num()].x, iz = ompLoopI[omp_get_thread_num()].z;
            size_t ig = box.ig({ix, iy, iz});
            int og = (new_i.x + 4 - ompLoopOffset.x)%box.n.x%4 + ((new_i.y + 4 - ompLoopOffset.y)%box.n.y%4 + ((new_i.z + 4 - ompLoopOffset.z)%box.n.z%4)*4)*4;
            int it = ompLoopIt[omp_get_thread_num()];
            if((ix != new_i.x)||(iy != new_i.y)||(iz != new_i.z))
            {
                if((cyclicDist(ix, new_i.x, box.n.x) > 1)||(cyclicDist(iy, new_i.y, box.n.y) > 1)||(cyclicDist(iz, new_i.z, box.n.z) > 1)){
                    data[ig][it]->addParticle(P, true);
                    violationList[omp_get_thread_num()].push_back({ig, it, data[ig][it]->P.size() - 1});
                } else {
                    bool destinationIsUnprocessedCell = (new_i.x + 4 - ompLoopOffset.x)%box.n.x%4 + ((new_i.y + 4 - ompLoopOffset.y)%box.n.y%4 + ((new_i.z + 4 - ompLoopOffset.z)%box.n.z%4)*4)*4 > ompLoopOg;
                    size_t new_ig = box.ig(new_i);
                    checkInitCell(new_ig, it);
                    data[new_ig][it]->addParticle(P, destinationIsUnprocessedCell);
                }
            } else {
                bool destinationIsUnprocessedCell = (typeIndex >= ompLoopIt[omp_get_thread_num()]);
                size_t new_ig = box.ig(new_i);
                checkInitCell(new_ig, typeIndex);
                data[new_ig][typeIndex]->addParticle(P, destinationIsUnprocessedCell);
            }            
            particleCount[omp_get_thread_num()]++;
        }
    }
    void ompLoop_removeCurrentParticle() // sets the flag to remove the current particle after its processing is finished
    {
        if(ompLoopIsRunning)
            ompLoop_flagDeleteCurrent[omp_get_thread_num()] = true;
        else
            pipic_log.message("ensemble.ompLoop_removeCurrentParticle() error: this function is to be only used during ompLoop.", true);
    }
    template <typename particleHandler>
    void ompLoop(particleHandler *handler) // handler must have a function processParticle(particle &P, particleType &type)
    {
        if(!shuffle_doubleLoop){
            ompLoop_(handler);
        } else {
            shuffle_doubleLoopStage = 0;
            ompLoop_(handler);
            shuffle_doubleLoopStage = 1;
            ompLoop_(handler);
        }
    }
    template <typename particleHandler>
    void ompLoop_(particleHandler *handler) // wrapper to enable doubleLoop shuffle 
    {
        ompLoopIsRunning = true;
        for(int iTh = 0; iTh < omp_get_max_threads(); iTh++) particleCount[iTh] = 0;
        
        ompLoopOffset = {rand()%4, rand()%4, rand()%4};
        // non-parallel loop over shifts to make possible parallel loop with stride {4, 4, 4}
        for(int oz = 0; oz < 4; oz++)
        for(int oy = 0; oy < 4; oy++)
        for(int ox = 0; ox < 4; ox++)
        {
            ompLoopOg = ox + (oy + oz*4)*4;
            #pragma omp parallel for collapse(3)
            for(int ioz = oz; ioz < box.n.z; ioz += 4)
            for(int ioy = oy; ioy < box.n.y; ioy += 4)
            for(int iox = ox; iox < box.n.x; iox += 4)
            {
                int ix = (iox + ompLoopOffset.x)%box.n.x, iy = (ioy + ompLoopOffset.y)%box.n.y, iz = (ioz + ompLoopOffset.z)%box.n.z;
                ompLoopI[omp_get_thread_num()] = {ix, iy, iz};
                size_t ig = box.ig({ix, iy, iz});
                if(data[ig] != nullptr)
                for(int it = 0; it < type.size(); it++)
                if(data[ig][it] != nullptr)
                for(int ip = data[ig][it]->begin(shuffle_cellShuffleProbability, shuffle_doubleLoop, shuffle_doubleLoopStage); !data[ig][it]->end(ip); ip++)
                {
                    particle *Particle = &(data[ig][it]->P[ip]);
                    ompLoopIt[omp_get_thread_num()] = it;
                    ompLoop_flagDeleteCurrent[omp_get_thread_num()] = false;
                    handler->processParticle(*Particle, type[it]);
                    if(ompLoop_flagDeleteCurrent[omp_get_thread_num()]) // removing the current particle
                    {
                         data[ig][it]->removeParticle(ip);
                         particleCount[omp_get_thread_num()]--;
                    } else {
                        if(periodicCoordinates)
                        {
                            Particle->r.x += (box.max.x-box.min.x)*((Particle->r.x < box.min.x) - (Particle->r.x > box.max.x));
                            Particle->r.y += (box.max.y-box.min.y)*((Particle->r.y < box.min.y) - (Particle->r.y > box.max.y));
                            Particle->r.z += (box.max.z-box.min.z)*((Particle->r.z < box.min.z) - (Particle->r.z > box.max.z));
                        }
                        int3 new_i = int3_coord(Particle->r);
                        if((ix != new_i.x)||(iy != new_i.y)||(iz != new_i.z))
                        {
                            if((cyclicDist(ix, new_i.x, box.n.x) > 1)||(cyclicDist(iy, new_i.y, box.n.y) > 1)||(cyclicDist(iz, new_i.z, box.n.z) > 1))
                                violationList[omp_get_thread_num()].push_back({ig, it, ip});
                            else {
                                bool destinationIsUnprocessedCell = (new_i.x + 4 - ompLoopOffset.x)%box.n.x%4 + ((new_i.y + 4 - ompLoopOffset.y)%box.n.y%4 + ((new_i.z + 4 - ompLoopOffset.z)%box.n.z%4)*4)*4 >= ompLoopOg;
                                size_t new_ig = box.ig(new_i);
                                checkInitCell(new_ig, it);
                                data[ig][it]->sendParticle(ip, *data[new_ig][it], destinationIsUnprocessedCell);
                            }
                        }
                    }
                }
            }
        }
        for(int iTh = 0; iTh < omp_get_max_threads(); iTh++) totalNumberOfParticles += particleCount[iTh];
        ompLoopIsRunning = false;
        
        int multiCellRelocated = 0;
        for(int iTh = 0; iTh < omp_get_max_threads(); iTh++)
        {
            for(int il = violationList[iTh].size() - 1; il >= 0; il--)
            {
                int ig =  violationList[iTh][il].x;
                int it =  violationList[iTh][il].y;
                int ip =  violationList[iTh][il].z;
                addParticle((*data[ig][it])[ip], it); totalNumberOfParticles--;
                data[ig][it]->removeParticle(ip);
                multiCellRelocated++;
            }
            violationList[iTh].clear();
        }
        if(multiCellRelocated > 0) pipic_log.message("ensemble.ompLoop warning: relocation for more than a cell (" + to_string(multiCellRelocated) + " particle(s)); consider reducing time step.");
    }
    struct nonOmpIterator //forward iterator for use without OMP 
    {
        nonOmpIterator(int it, ensemble *Ensemble) : ig(0), ip(-1), it(it), Ensemble(Ensemble) { (*this)++; }
        particle& operator*() const { return *Particle; }
        friend bool operator < (nonOmpIterator const& lhs, int const& rhs)
        {
            return (lhs.Particle != nullptr);
        }
        void removeCurrentParticle()
        {
            if(Ensemble->ompLoopIsRunning){
                pipic_log.message("ensemble::nonOmpIterator error: nonOmpIterator cannot be used during ompLoop.", true);
            } else {
                Ensemble->data[ig][it]->removeParticle(ip);
                Ensemble->totalNumberOfParticles--;
                Particle = nullptr;
            }
        }
        bool operator++(int)
        {
            Particle = nullptr;
            if(Ensemble->ompLoopIsRunning){
                pipic_log.message("ensemble::nonOmpIterator error: nonOmpIterator cannot be used during ompLoop.", true);
                return false;
            }
            while(Particle == nullptr)
            {
                if(Ensemble->data[ig] == nullptr) {
                    if(ig < Ensemble->box.ng - 1) ig++;
                        else return false;
                } else {
                    if(Ensemble->data[ig][it] == nullptr){
                    if(ig < Ensemble->box.ng - 1) ig++;
                        else return false;
                    } else {
                        if(ip + 1 < Ensemble->data[ig][it]->P.size()){
                            ip++;
                            Particle = &(Ensemble->data[ig][it]->P[ip]);
                            return true;
                        } else {
                            ip = -1;
                            if(ig < Ensemble->box.ng - 1) ig++;
                            else return false;
                        }
                    }
                }
            }
            return false;
        }
        bool getLocation(int3 &indGrid, int &indParticle)
        {
            if(Particle != nullptr)
            {
                indGrid = Ensemble->box.ind(ig);
                indParticle = ip;
                return true;
            } else 
                return false;
        }
        size_t ig;
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
        for(int it = 0; it < type.size(); it++) if(typeName == type[it].name) typeIndex = it;
        if(typeIndex == -1){
            pipic_log.message("pipic error: unknown particle type name '" + typeName + "'.", true);
            exit(-1);
        }
        return typeIndex;
    }
};

#endif