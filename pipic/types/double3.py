"""
-------------------------------------------------------------------------------------------------------
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
"""

import numba as nb
import numpy as np
from numba.experimental import jitclass

# __all__ is used to not pollute the namespace when using 'from xyz import *'
__all__ = ["Double3"]


@jitclass(
    [
        ("x", nb.types.double),
        ("y", nb.types.double),
        ("z", nb.types.double),
    ]
)
class Double3:
    def __init__(self, x=0, y=0, z=0):
        self.x = x
        self.y = y
        self.z = z

    def normalize(self):
        r = np.sqrt(self.x * self.x + self.y * self.y + self.z * self.z)
        if r > 0:
            inv_r = 1 / r
            self.x *= inv_r
            self.y *= inv_r
            self.z *= inv_r

    def norm(self) -> float:
        return np.sqrt(self.x * self.x + self.y * self.y + self.z * self.z)

    def norm2(self) -> float:
        return self.x * self.x + self.y * self.y + self.z * self.z

    def __add__(self, other):
        return Double3(self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__(self, other):
        return Double3(self.x - other.x, self.y - other.y, self.z - other.z)

    def __iadd__(self, other):
        self.x += other.x
        self.y += other.y
        self.z += other.z
        return self

    def __mul__(self, other):
        return Double3(self.x * other, self.y * other, self.z * other)

    def __imul__(self, other):
        self.x *= other
        self.y *= other
        self.z *= other
        return self

    def cross(self, b):
        return Double3(
            self.y * b.z - self.z * b.y,
            self.z * b.x - self.x * b.z,
            self.x * b.y - self.y * b.x,
        )

    def dot(self, b) -> float:
        return self.x * b.x + self.y * b.y + self.z * b.z
