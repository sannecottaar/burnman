from __future__ import absolute_import

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1, os.path.abspath('..'))

import burnman
from burnman import minerals

# P_0 = 20 GPa and T_0 = 1500 K
P0 = 0.
T0 = 0.
Pref = 20.e9
Tref = 1500.

E0_13 = -100885.
E0_23 = -49.e3
V0_13 = 8.8e-8
V0_23 = 1.9e-7
S0_13 = 4.5
S0_23 = 0.
dSdT_13 = 17.e-3
dSdT_23 = 0.

print([E0_13 - Tref*S0_13 - Tref*Tref*dSdT_13/2. + Pref*V0_13,
       E0_23 - Tref*S0_23 - Tref*Tref*dSdT_23/2. + Pref*V0_23])
print([S0_13 + Tref*dSdT_13, S0_23 + Tref*dSdT_23])

melt = burnman.SolidSolution(name = 'FMS melt',
                             solution_type = 'boukare',
                             endmembers = [[minerals.DKS_2013_liquids.MgO_liquid(), '[Fe]O'],
                                           [minerals.DKS_2013_liquids.MgO_liquid(), '[Mg]O'],
                                           [minerals.DKS_2013_liquids.SiO2_liquid(), '[Si]O2']],
                             energy_interaction = [[0., E0_23],
                                                   [E0_13]],
                             volume_interaction = [[0., V0_23],
                                                   [V0_13]],
                             entropy_interaction = [[0., S0_23],
                                                    [S0_13]],
                             dentropydT_interaction = [[0., dSdT_23],
                                                       [dSdT_13]],
                             alphas = [1., 1., 1.],
                             molar_fractions = [0.5, 0.5, 0.0])


print(melt)
