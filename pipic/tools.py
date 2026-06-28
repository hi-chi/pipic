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

# .tools is used as shorthand for pulling everything into the namespace
from . import consts as _consts
from . import ctypes as _ctypes
from . import interfaces as _interfaces
from . import types as _types

__all__ = []

for _module in (_consts, _ctypes, _interfaces, _types):
    for _name in _module.__all__:
        if _name in _types._LEGACY_CALLBACK_NAMES:
            continue
        globals()[_name] = getattr(_module, _name)
        __all__.append(_name)

__all__.extend(_types._LEGACY_CALLBACK_NAMES)


def __getattr__(name):
    if name in _types._LEGACY_CALLBACK_NAMES:
        replacement = _types._LEGACY_CALLBACK_NAMES[name]
        _types._warn_legacy_callback_name(name, replacement, stacklevel=3)
        return getattr(_types, replacement)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
