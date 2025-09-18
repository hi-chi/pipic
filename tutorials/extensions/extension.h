// struct for managing the thread-specific data.
struct threadHandler{
    mt19937 rng;
    std::uniform_real_distribution<double> U1;
    std::normal_distribution<double> N1;
    threadHandler(): U1(0, 1.0), N1(0, 1.0) {}
    double random() {return U1(rng);} // returns a random number from [0, 1)
    double nrandom() {return N1(rng);} // returns a normal random number from [0, 1)
};

// Function for adding new particles to the cell interface.
void addParticle(cellInterface &CI, particle &P){
    if(CI.particleBufferSize < CI.particleBufferCapacity){ // checking if the buffer permits adding a particle 
        *CI.newParticle(CI.particleBufferSize) = P; // copy particle to a new particle (buffer)
        CI.particleBufferSize++;
    } else {
        pipic_log.message("pi-PIC error: particle buffer overflow.", true);
    };
};

// Function for removing particles according to some probability rate.
void removeParticles(cellContainer* C, threadHandler &cthread, double rate){ 
    int particles = int(C->P.size());
    // iterate over particles in cell. The C->endShift particles from the end of the list 
    // came from other cells during the paricle push and should thus not be processed again.
    for (int ip = 0; ip < (particles - C->endShift); ip++){
        if (cthread.random() > rate) {
            C->P[ip].w = 0.; // set particle weight to zero to mark it for removal
        }
    }
};