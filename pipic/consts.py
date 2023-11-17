import _pipic
import numpy as np

# __all__ is used to not pollute the namespace when using 'from xyz import *'
__all__ = ['light_velocity',
           'electron_mass',
           'electron_charge',
           'elementary_charge',
           'plancks_constant',
           'hbar',
           'boltzmann_constant',
           'avogadro_constant']


light_velocity = _pipic.light_velocity
electron_mass = _pipic.electron_mass
electron_charge = _pipic.electron_charge

elementary_charge = np.abs(electron_charge)

plancks_constant = 6.62607015e-27  # 2019 definition
hbar = plancks_constant/(2*np.pi)

boltzmann_constant = 1.380649e-16  # 2019 definition
avogadro_constant = 6.02214076e23  # 2019 definition
