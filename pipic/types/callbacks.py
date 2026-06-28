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

import warnings

import numba as nb

# __all__ is used to not pollute the namespace when using 'from xyz import *'
__all__ = [
    "add_particles",
    "particle_loop",
    "field_loop",
    "it2r",
    "custom_field_loop",
    "add_particles_callback",
    "particle_loop_callback",
    "field_loop_callback",
    "it2r_callback",
    "custom_field_loop_callback",
    "handler_callback",
]

_LEGACY_CALLBACK_NAMES = {
    "add_particles_callback": "add_particles",
    "particle_loop_callback": "particle_loop",
    "field_loop_callback": "field_loop",
    "it2r_callback": "it2r",
    "custom_field_loop_callback": "custom_field_loop",
}


def _warn_legacy_callback_name(name, replacement, stacklevel=2):
    warnings.warn(
        f"`{name}` is deprecated and will be removed in a future release. "
        f"Use `{replacement}` instead.",
        FutureWarning,
        stacklevel=stacklevel,
    )


add_particles = nb.types.double(
    nb.types.CPointer(nb.types.double),
    nb.types.CPointer(nb.types.double),
    nb.types.CPointer(nb.types.int32),
)

particle_loop = nb.types.void(
    nb.types.CPointer(nb.types.double),
    nb.types.CPointer(nb.types.double),
    nb.types.CPointer(nb.types.double),
    nb.types.CPointer(nb.types.uint64),
    nb.types.CPointer(nb.types.double),
    nb.types.CPointer(nb.types.int32),
)

field_loop = nb.types.void(
    nb.types.CPointer(nb.types.int32),
    nb.types.CPointer(nb.types.double),
    nb.types.CPointer(nb.types.double),
    nb.types.CPointer(nb.types.double),
    nb.types.CPointer(nb.types.double),
    nb.types.CPointer(nb.types.int32),
)

it2r = nb.types.void(
    nb.types.CPointer(nb.types.int32),
    nb.types.CPointer(nb.types.double),
    nb.types.CPointer(nb.types.double),
    nb.types.CPointer(nb.types.int32),
)

custom_field_loop = nb.types.void(
    nb.types.CPointer(nb.types.int32),
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


def __getattr__(name):
    if name in _LEGACY_CALLBACK_NAMES:
        replacement = _LEGACY_CALLBACK_NAMES[name]
        _warn_legacy_callback_name(name, replacement, stacklevel=3)
        return globals()[replacement]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
