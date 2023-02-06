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
// Description of the file: File introduces functions for CIC weighting in 1D/2D/3D geometry.

#ifndef CIC_WEIGHTING_H
#define CIC_WEIGHTING_H

#include "emfield.h"

inline void CIC1d(double &Lw, int &Li, double coord, double min, double invSize, int n) // 1D CIC: returns weight (Lw) and index (Li) of the low node; coord is coordinate, n is the mnumber of cells
{
    double u = n*(coord - min)*invSize - 0.5;
    Lw = 1 - (u - floor(u));
    Li = ((int)floor(u) % n + n) % n;
}

void CIC(double3 &r, simulationBox &box, double *c, size_t *cig) // for a given coordinate r determines contributing nodes and compute weights
{
    int Lx = 0, Ly = 0, Lz = 0, Ux = 0, Uy = 0, Uz = 0; // indeces of lower (L) and upper (U) bound nodes
    
    double wx;
    CIC1d(wx, Lx, r.x, box.min.x, box.invSize.x, box.n.x);
    Ux = (Lx + 1) % box.n.x;
    c[0] = wx;
    c[1] = (1 - wx);
    
    if(box.n.y != 1){ // 2D or 3D
        double wy;
        CIC1d(wy, Ly, r.y, box.min.y, box.invSize.y, box.n.y);
        Uy = (Ly + 1) % box.n.y;
        for(int i = 0; i < 2; i++) c[i+2] = c[i]*(1 - wy);
        for(int i = 0; i < 2; i++) c[i] *= wy;
    }

    if(box.n.z != 1){ // 3D 
        double wz;
        CIC1d(wz, Lz, r.z, box.min.z, box.invSize.z, box.n.z);
        Uz = (Lz + 1) % box.n.z;
        for(int i = 0; i < 4; i++) c[i+4] = c[i]*(1 - wz);
        for(int i = 0; i < 4; i++) c[i] *= wz;
    }

    cig[0] = box.ig({Lx, Ly, Lz});
    cig[1] = box.ig({Ux, Ly, Lz});
    if(box.n.y != 1){ // 2D or 3D
        cig[2] = box.ig({Lx, Uy, Lz});
        cig[3] = box.ig({Ux, Uy, Lz});
        if(box.n.z != 1){ // 3D
            cig[4] = box.ig({Lx, Ly, Uz});
            cig[5] = box.ig({Ux, Ly, Uz});
            cig[6] = box.ig({Lx, Uy, Uz});
            cig[7] = box.ig({Ux, Uy, Uz});
        }
    }
}
#endif