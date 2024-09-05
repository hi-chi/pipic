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

import numpy as np

from . import _pipic

# __all__ is used to not pollute the namespace when using 'from xyz import *'
__all__ = [
    "light_velocity",
    "electron_mass",
    "electron_charge",
    "elementary_charge",
    "plancks_constant",
    "hbar",
    "boltzmann_constant",
    "avogadro_constant",
]


light_velocity = _pipic.light_velocity
electron_mass = _pipic.electron_mass
electron_charge = _pipic.electron_charge

elementary_charge = np.abs(electron_charge)

plancks_constant = 6.62607015e-27  # 2019 definition
hbar = plancks_constant / (2 * np.pi)

boltzmann_constant = 1.380649e-16  # 2019 definition
avogadro_constant = 6.02214076e23  # 2019 definition
