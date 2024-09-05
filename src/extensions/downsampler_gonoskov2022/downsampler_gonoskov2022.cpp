/*-------------------------------------------------------------------------------------------------------
downsampler_gonoskov2022, Copyright 2024 Arkady Gonoskov
---------------------------------------------------------------------------------------------------------
downsampler_gonoskov2022 is free software: you can redistribute it and/or modify it under the terms of the GNU General
Public License as published by the Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

downsampler_gonoskov2022 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with pi-PIC. If not, se
<https://www.gnu.org/licenses/>.
---------------------------------------------------------------------------------------------------------
Website: https://github.com/hi-chi/pipic
Contact: arkady.gonoskov@gu.se.
-------------------------------------------------------------------------------------------------------*/
// Description: An implementation of agnostic conservative downsampling based on the method described in
// [A. Gonoskov, Comput. Phys. Commun. 271 (2022)]. The extension considers each subset of particles
// that all give CIC contributions to a set of 2/4/8 nearby nodes, depending on dimensionality.
// The implementation is configured to always preserve total weight of particles in each subset.
// In addition, the energy, momentum and CIC contributions are preserved (can be switched off).
// The extension can act on several types of particles; this is to be configured with add_assignment().

#include "interfaces.h"
#include <pybind11/pybind11.h>
#include "pybind11/stl.h"
#include <pybind11/operators.h>

const string name = "downsampler_gonoskov2022";

//Warning: the implementation of this extension exploits accessing raw data from nearby cells,
// which is not provided by the cellInterface at the moment of development.

struct cellContainer
{
    vector<particle> P;
    int endShift; // the number of particles from the end of the vector that have been already processed (came from other cells) during the ongoing OMP loop
    cellContainer(): endShift(0) {}
};
static cellContainer ***cell;

struct Spec{ // scpecification of downsampling parameters for each type of interest
    int typeIndex;
    bool preserverEnergy, preserveMomentum, preserveCICWeight;
    int cap;
    double targetRatio;
    int numOfConstraints(int dim){
        return (!preserveCICWeight) + preserverEnergy + 3*preserveMomentum + (1 << dim)*preserveCICWeight;
    }
};
static vector<Spec> spec;

struct threadHandler{
    struct particleRef{
        particle* P;
        double e[12]; // coefficients of constraints
        particleRef(particle *P): P(P) {}
    };
    double v[13][13];
    vector<particleRef> PR;
    vector<particleRef*> PP;
    int dim;
    mt19937 rng;
    std::uniform_real_distribution<double> U1;
    threadHandler(): U1(0, 1.0){
        PR.reserve(100);
        PP.reserve(100);
    }
    double random() {return U1(rng);}
    void downsampleSingle(int is, cellInterface &C){
        int m = spec[is].numOfConstraints(C.dim);
        int n = m + 1; // apply algorithm to first m + 1 particles
	if(n > int(PP.size())) return;
        for(int in = 0; in < n; in++){
            for(int jn = 0; jn < n; jn++) v[in][jn] = 0;
            v[in][in] = 1;
        }
        int vnum = n;
        for(int im = 0; im < m; im++){
            for(int ivn = 0; ivn < vnum - 1; ivn++){
                double k1 = 0; for(int in = 0; in < n; in++) k1 += v[ivn+1][in]*PP[in]->e[im];
                double k2 = 0; for(int in = 0; in < n; in++) k2 -= v[ivn+0][in]*PP[in]->e[im];
                if ((k1 != 0) || (k2 != 0)) {
                    double s = 0;
                    for(int in = 0; in < n; in++){
                        v[ivn][in] = k1*v[ivn][in] + k2*v[ivn+1][in];
                        s += sqr(v[ivn][in]);
                    }
                    if(s == 0) return; // cancel downsampling for this case
                    s = sqrt(s);
                    for(int in = 0; in < n; in++) v[ivn][in] /= s;

                    //a check on limitations of numerical arithmetic:
                    s = 0; for(int in = 0; in < n; in++) s += sqr(v[ivn][in]);
                    if(abs(s - 1.0) > 1e-11) return;
                }
            }
            vnum--;
        }
        double ap = 0, an = 0; //positive and negrative values of a that are closest to 0;
		int iap = 0, ian = 0; // the corresponding indices;
		bool fn = false, fp = false; // found negative, found positive
		for(int in = 0; in < n; in++){
			if(v[0][in] != 0){
				double a = -PP[in]->P->w / v[0][in];
				if(a > 0){
					if(!fp){
						fp = true;
						ap = a;
						iap = in;
					}
					if(a < ap){
						ap = a;
						iap = in;
					}
				}
				if(a < 0){
					if(!fn){
						fn = true;
						an = a;
						ian = in;
					}
					if(a > an){
						an = a;
						ian = in;
					}
				}
			}
		}
		if((!fp)||(!fn)) return;
        double threshold = 1e-13;
        double w_max = 0; for(int in = 0; in < n; in++) if(w_max < PP[in]->P->w) w_max = PP[in]->P->w;
        if(random() < ap/(ap - an)){
            for(int in = 0; in < n; in++){
                PP[in]->P->w += an*v[0][in];
                if((in == ian)||(PP[in]->P->w < w_max*threshold)) PP[in]->P->w = 0;
            }
		}else{
            for(int in = 0; in < n; in++){
                PP[in]->P->w += ap*v[0][in];
                if((in == iap)||(PP[in]->P->w < w_max*threshold)) PP[in]->P->w = 0;
            }
		}

        for(int ip = int(PP.size()) - 1; ip >= 0; ip--)
        if(PP[ip]->P->w == 0){
            PP[ip] = PP[PP.size() - 1];
            PP.pop_back();
        }
    }
    void addParticles(int is, int3 i, int3 n, double3 min, double3 step){
        PR.clear();
        int it = spec[is].typeIndex;

        dim = 3;
        if(n.z == 1) dim = 2;
        if((n.y == 1)&&(n.z == 1)) dim = 1;

        bool exclProcessed = false; // if true, the cell::endShift is accounted for to exclude processed particles

        int ix = i.x, iy = i.y, iz = i.z;
        int ig;
        ig = ix + (iy + iz*n.y)*n.x;
        if(cell[ig] != nullptr)
            if(cell[ig][it] != nullptr)
                for(int ip = 0; ip < int(cell[ig][it]->P.size()) - exclProcessed*cell[ig][it]->endShift; ip++)if(cell[ig][it]->P[ip].w != 0){
                    particle *PC = &(cell[ig][it]->P[ip]); // particle candidate
                    bool add = true;
                    if(PC->r.x < min.x + step.x*(ix + 0.5)) add = false;
                    if(dim > 1) if(PC->r.y < min.y + step.y*(iy + 0.5)) add = false;
                    if(dim > 2) if(PC->r.z < min.z + step.z*(iz + 0.5)) add = false;
                    if(add) PR.emplace_back(PC);
                }
        ix++; if (ix == n.x) ix = 0;
        ig = ix + (iy + iz*n.y)*n.x;
        if(cell[ig] != nullptr)
            if(cell[ig][it] != nullptr)
                for(int ip = 0; ip < int(cell[ig][it]->P.size()) - exclProcessed*cell[ig][it]->endShift; ip++)if(cell[ig][it]->P[ip].w != 0){
                    particle *PC = &(cell[ig][it]->P[ip]); // particle candidate
                    bool add = true;
                    if(PC->r.x >= min.x + step.x*(ix + 0.5)) add = false;
                    if(dim > 1) if(PC->r.y < min.y + step.y*(iy + 0.5)) add = false;
                    if(dim > 2) if(PC->r.z < min.z + step.z*(iz + 0.5)) add = false;
                    if(add) PR.emplace_back(PC);
                }
        if(n.y == 1) return;
        ix = i.x;
        iy++; if (iy == n.y) iy = 0;
        ig = ix + (iy + iz*n.y)*n.x;
        if(cell[ig] != nullptr)
            if(cell[ig][it] != nullptr)
                for(int ip = 0; ip < int(cell[ig][it]->P.size()) - exclProcessed*cell[ig][it]->endShift; ip++)if(cell[ig][it]->P[ip].w != 0){
                    particle *PC = &(cell[ig][it]->P[ip]); // particle candidate
                    bool add = true;
                    if(PC->r.x < min.x + step.x*(ix + 0.5)) add = false;
                    if(dim > 1) if(PC->r.y >= min.y + step.y*(iy + 0.5)) add = false;
                    if(dim > 2) if(PC->r.z < min.z + step.z*(iz + 0.5)) add = false;
                    if(add) PR.emplace_back(PC);
                }
        ix++; if (ix == n.x) ix = 0;
        ig = ix + (iy + iz*n.y)*n.x;
        if(cell[ig] != nullptr)
            if(cell[ig][it] != nullptr)
                for(int ip = 0; ip < int(cell[ig][it]->P.size()) - exclProcessed*cell[ig][it]->endShift; ip++)if(cell[ig][it]->P[ip].w != 0){
                    particle *PC = &(cell[ig][it]->P[ip]); // particle candidate
                    bool add = true;
                    if(PC->r.x >= min.x + step.x*(ix + 0.5)) add = false;
                    if(dim > 1) if(PC->r.y >= min.y + step.y*(iy + 0.5)) add = false;
                    if(dim > 2) if(PC->r.z < min.z + step.z*(iz + 0.5)) add = false;
                    if(add) PR.emplace_back(PC);
                }
        if(n.z == 1) return;
        ix = i.x; iy = i.y;
        iz++; if (iz == n.z) iz = 0;
        ig = ix + (iy + iz*n.y)*n.x;
        if(cell[ig] != nullptr)
            if(cell[ig][it] != nullptr)
                for(int ip = 0; ip < int(cell[ig][it]->P.size()) - exclProcessed*cell[ig][it]->endShift; ip++)if(cell[ig][it]->P[ip].w != 0){
                    particle *PC = &(cell[ig][it]->P[ip]); // particle candidate
                    bool add = true;
                    if(PC->r.x < min.x + step.x*(ix + 0.5)) add = false;
                    if(dim > 1) if(PC->r.y < min.y + step.y*(iy + 0.5)) add = false;
                    if(dim > 2) if(PC->r.z >= min.z + step.z*(iz + 0.5)) add = false;
                    if(add) PR.emplace_back(PC);
                }
        ix++; if (ix == n.x) ix = 0;
        ig = ix + (iy + iz*n.y)*n.x;
        if(cell[ig] != nullptr)
            if(cell[ig][it] != nullptr)
                for(int ip = 0; ip < int(cell[ig][it]->P.size()) - exclProcessed*cell[ig][it]->endShift; ip++)if(cell[ig][it]->P[ip].w != 0){
                    particle *PC = &(cell[ig][it]->P[ip]); // particle candidate
                    bool add = true;
                    if(PC->r.x >= min.x + step.x*(ix + 0.5)) add = false;
                    if(dim > 1) if(PC->r.y < min.y + step.y*(iy + 0.5)) add = false;
                    if(dim > 2) if(PC->r.z >= min.z + step.z*(iz + 0.5)) add = false;
                    if(add) PR.emplace_back(PC);
                }
        ix = i.x;
        iy++; if (iy == n.y) iy = 0;
        ig = ix + (iy + iz*n.y)*n.x;
        if(cell[ig] != nullptr)
            if(cell[ig][it] != nullptr)
                for(int ip = 0; ip < int(cell[ig][it]->P.size()) - exclProcessed*cell[ig][it]->endShift; ip++)if(cell[ig][it]->P[ip].w != 0){
                    particle *PC = &(cell[ig][it]->P[ip]); // particle candidate
                    bool add = true;
                    if(PC->r.x < min.x + step.x*(ix + 0.5)) add = false;
                    if(dim > 1) if(PC->r.y >= min.y + step.y*(iy + 0.5)) add = false;
                    if(dim > 2) if(PC->r.z >= min.z + step.z*(iz + 0.5)) add = false;
                    if(add) PR.emplace_back(PC);
                }
        ix++; if (ix == n.x) ix = 0;
        ig = ix + (iy + iz*n.y)*n.x;
        if(cell[ig] != nullptr)
            if(cell[ig][it] != nullptr)
                for(int ip = 0; ip < int(cell[ig][it]->P.size()) - exclProcessed*cell[ig][it]->endShift; ip++)if(cell[ig][it]->P[ip].w != 0){
                    particle *PC = &(cell[ig][it]->P[ip]); // particle candidate
                    bool add = true;
                    if(PC->r.x >= min.x + step.x*(ix + 0.5)) add = false;
                    if(dim > 1) if(PC->r.y >= min.y + step.y*(iy + 0.5)) add = false;
                    if(dim > 2) if(PC->r.z >= min.z + step.z*(iz + 0.5)) add = false;
                    if(add) PR.emplace_back(PC);
                }
    }
    void set_constraints(int is, cellInterface &C){
        double m2c2_ = 1/sqr(C.particleMass*lightVelocity);
        double mc2 = C.particleMass*sqr(lightVelocity);
        for(int ip = 0; ip < int(PR.size()); ip++){
            particle &P(*PR[ip].P);
            int ie = 0;
            if(!spec[is].preserveCICWeight){
                PR[ip].e[ie] = 1; ie++;
            }
            if(spec[is].preserverEnergy){
                if(C.particleMass == 0){
                    PR[ip].e[ie] = lightVelocity*P.p.norm();
                } else {
                    double d = m2c2_*P.p.norm2();
                    if(d > 1e-7) PR[ip].e[ie] = mc2*(sqrt(1 + d) - 1);
                    else PR[ip].e[ie] = 0.5*mc2*d;
                }
                ie++;
            }
            if(spec[is].preserveMomentum){
                PR[ip].e[ie] = P.p.x; ie++;
                PR[ip].e[ie] = P.p.y; ie++;
                PR[ip].e[ie] = P.p.z; ie++;
            }
            if(spec[is].preserveCICWeight){
                if(dim == 1){
                    double xf = (C.globalMin.x + (C.i.x + 1.5)*C.step.x - P.r.x)*C.invStep.x;
                    PR[ip].e[ie] = xf; ie++;
                    PR[ip].e[ie] = (1 - xf); ie++;
                }
                if(dim == 2){
                    double xf = (C.globalMin.x + (C.i.x + 1.5)*C.step.x - P.r.x)*C.invStep.x;
                    double yf = (C.globalMin.y + (C.i.y + 1.5)*C.step.y - P.r.y)*C.invStep.y;
                    PR[ip].e[ie] = yf*xf; ie++;
                    PR[ip].e[ie] = yf*(1 - xf); ie++;
                    PR[ip].e[ie] = (1 - yf)*xf; ie++;
                    PR[ip].e[ie] = (1 - yf)*(1 - xf); ie++;
                }
                if(dim == 3){
                    double xf = (C.globalMin.x + (C.i.x + 1.5)*C.step.x - P.r.x)*C.invStep.x;
                    double yf = (C.globalMin.y + (C.i.y + 1.5)*C.step.y - P.r.y)*C.invStep.y;
                    double zf = (C.globalMin.z + (C.i.z + 1.5)*C.step.z - P.r.z)*C.invStep.z;
                    PR[ip].e[ie] = zf*yf*xf; ie++;
                    PR[ip].e[ie] = zf*yf*(1 - xf); ie++;
                    PR[ip].e[ie] = zf*(1 - yf)*xf; ie++;
                    PR[ip].e[ie] = zf*(1 - yf)*(1 - xf); ie++;
                    PR[ip].e[ie] = (1 - zf)*yf*xf; ie++;
                    PR[ip].e[ie] = (1 - zf)*yf*(1 - xf); ie++;
                    PR[ip].e[ie] = (1 - zf)*(1 - yf)*xf; ie++;
                    PR[ip].e[ie] = (1 - zf)*(1 - yf)*(1 - xf); ie++;
                }
            }
        }
    }
    int hypoteticalMax(cellInterface &C, int it){ // computes the sum of P.size() of involved cells
        int hMax = 0, ig, ix = C.i.x, iy = C.i.y, iz = C.i.z;
        ig = ix + (iy + iz*C.n.y)*C.n.x;
        if(cell[ig] != nullptr) if(cell[ig][it] != nullptr) hMax += cell[ig][it]->P.size();
        ix++; if(ix == C.n.x) ix = 0; ig = ix + (iy + iz*C.n.y)*C.n.x;
        if(cell[ig] != nullptr) if(cell[ig][it] != nullptr) hMax += cell[ig][it]->P.size();
        if(C.dim == 1) return hMax;
        iy++; if(iy == C.n.y) iy = 0; ig = ix + (iy + iz*C.n.y)*C.n.x;
        if(cell[ig] != nullptr) if(cell[ig][it] != nullptr) hMax += cell[ig][it]->P.size();
        ix--; if(ix == -1) ix = C.n.x - 1; ig = ix + (iy + iz*C.n.y)*C.n.x;
        if(cell[ig] != nullptr) if(cell[ig][it] != nullptr) hMax += cell[ig][it]->P.size();
        if(C.dim == 2) return hMax;
        iz++; if(iz == C.n.z) iz = 0; ig = ix + (iy + iz*C.n.y)*C.n.x;
        if(cell[ig] != nullptr) if(cell[ig][it] != nullptr) hMax += cell[ig][it]->P.size();
        ix++; if(ix == C.n.x) ix = 0; ig = ix + (iy + iz*C.n.y)*C.n.x;
        if(cell[ig] != nullptr) if(cell[ig][it] != nullptr) hMax += cell[ig][it]->P.size();
        iy++; if (iy == C.n.y) iy = 0; ig = ix + (iy + iz*C.n.y)*C.n.x;
        if(cell[ig] != nullptr) if(cell[ig][it] != nullptr) hMax += cell[ig][it]->P.size();
        ix--; if (ix == -1) ix = C.n.x - 1; ig = ix + (iy + iz*C.n.y)*C.n.x;
        if(cell[ig] != nullptr) if(cell[ig][it] != nullptr) hMax += cell[ig][it]->P.size();
        return hMax;
    }
    void downsample(int is, cellInterface &C){
        if((spec[is].numOfConstraints(C.dim) >= spec[is].cap)||(spec[is].numOfConstraints(C.dim) >= int(spec[is].targetRatio*spec[is].cap))){
            if(C.threadNum == 0){cout << name << ": Error: The number of constraints must be less than cap and int(target_ratio*cap)." << endl;}
            return;
        }
        if(hypoteticalMax(C, spec[is].typeIndex) <= spec[is].cap) return;
        addParticles(is, C.i, C.n, C.globalMin, C.step);
        if(int(PR.size()) <= spec[is].cap) return;
        set_constraints(is, C);
        PP.clear();
        for(int ip = 0; ip < int(PR.size()); ip++) PP.push_back(&PR[ip]);
        rng.seed(C.rngSeed);
        for(int ia = 0; ia < int(PR.size()) - int(spec[is].targetRatio*spec[is].cap); ia++){
            std::shuffle(begin(PP), next(begin(PP), PP.size()), rng);
            downsampleSingle(is, C);
        }
    }
};
static vector<threadHandler> Thread;

// function that is called to process particles in a given cell
void Handler(int *I, double *D, double *F, double *P, double *NP, double *dataDouble, int *dataInt){
    cellInterface C(I, D, F, P, NP);
    if(C.particleTypeIndex != -1){cout << name << ": Error: Inconsistent subject, set it to 'cells'." << endl; exit(0);}
    for(int is = 0; is < int(spec.size()); is++) Thread[C.threadNum].downsample(is, C);
};

// extension initialization
int64_t handler(int64_t ensembleData, int typeIndex, bool preserveEnergy = true,
                bool preserveMomentum = true, bool preserveCICWeight = true,
                int cap = 15, double targetRatio = 1.0){
    double _targetRatio = targetRatio;
    if(_targetRatio > 1.0){cout << name << ": Warning: target_ratio must be <= 1; setting to 1.0." << endl; _targetRatio = 1.0;}
    spec.push_back({typeIndex, preserveEnergy, preserveMomentum, preserveCICWeight, cap, _targetRatio});
    Thread.resize(omp_get_max_threads());
    cell = (cellContainer***)ensembleData;
    return (int64_t)Handler;
};

void addAssignment(int typeIndex, bool preserveEnergy = true,
                    bool preserveMomentum = true, bool preserveCICWeight = true,
                    int cap = 15, double targetRatio = 1.0){
    double _targetRatio = targetRatio;
    if(_targetRatio > 1.0){cout << name << ": Warning: target_ratio must be <= 1; setting to 1.0." << endl; _targetRatio = 1.0;}
    spec.push_back({typeIndex, preserveEnergy, preserveMomentum, preserveCICWeight, cap, _targetRatio});
};

namespace py = pybind11;
PYBIND11_MODULE(_downsampler_gonoskov2022, object) {
    object.attr("name") = name;
    object.def("handler", &handler, py::arg("ensemble_data"), py::arg("type_index"), py::arg("preserve_energy") = true,
               py::arg("preserve_momentum") = true, py::arg("preserve_cic_weight") = true,
               py::arg("cap") = 15, py::arg("target_ratio") = 1.0);
    object.def("add_assignment", &addAssignment, py::arg("type_index"), py::arg("preserve_energy") = true,
               py::arg("preserve_momentum") = true, py::arg("preserve_cic_weight") = true,
               py::arg("cap") = 15, py::arg("target_ratio") = 1.0);
}
