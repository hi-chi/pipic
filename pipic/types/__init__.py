# -*- coding: utf-8 -*-

from .callbacks import (
    _LEGACY_CALLBACK_NAMES,
    _warn_legacy_callback_name,
    add_particles,
    custom_field_loop,
    field_loop,
    handler_callback,
    it2r,
    particle_loop,
)
from .double3 import Double3
from .int3 import Int3

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
    "Int3",
    "Double3",
]


def __getattr__(name):
    if name in _LEGACY_CALLBACK_NAMES:
        replacement = _LEGACY_CALLBACK_NAMES[name]
        _warn_legacy_callback_name(name, replacement, stacklevel=3)
        return globals()[replacement]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
