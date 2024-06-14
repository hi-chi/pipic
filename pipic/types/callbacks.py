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

# __all__ is used to not pollute the namespace when using 'from xyz import *'
__all__ = [
    "add_particles_callback",
    "particle_loop_callback",
    "field_loop_callback",
    "it2r_callback",
    "field2data_callback",
    "handler_callback",
]


add_particles_callback = nb.types.double(
    nb.types.CPointer(nb.types.double),
    nb.types.CPointer(nb.types.double),
    nb.types.CPointer(nb.types.int32),
)

particle_loop_callback = nb.types.void(
    nb.types.CPointer(nb.types.double),
    nb.types.CPointer(nb.types.double),
    nb.types.CPointer(nb.types.double),
    nb.types.CPointer(nb.types.uint64),
    nb.types.CPointer(nb.types.double),
    nb.types.CPointer(nb.types.int32),
)

field_loop_callback = nb.types.void(
    nb.types.CPointer(nb.types.int32),
    nb.types.CPointer(nb.types.double),
    nb.types.CPointer(nb.types.double),
    nb.types.CPointer(nb.types.double),
    nb.types.CPointer(nb.types.double),
    nb.types.CPointer(nb.types.int32),
)

it2r_callback = nb.types.void(
    nb.types.CPointer(nb.types.int32),
    nb.types.CPointer(nb.types.double),
    nb.types.CPointer(nb.types.double),
    nb.types.CPointer(nb.types.int32),
)

field2data_callback = nb.types.void(
    nb.types.CPointer(nb.types.int32),
    nb.types.CPointer(nb.types.double),
    nb.types.CPointer(nb.types.double),
    nb.types.CPointer(nb.types.double),
    nb.types.CPointer(nb.types.double),
    nb.types.CPointer(nb.types.int32),
)

handler_callback = nb.types.void(
    nb.types.CPointer(nb.types.int32),
    nb.types.CPointer(nb.types.double),
    nb.types.CPointer(nb.types.double),
    nb.types.CPointer(nb.types.double),
    nb.types.CPointer(nb.types.double),
    nb.types.CPointer(nb.types.double),
    nb.types.CPointer(nb.types.int32),
)
