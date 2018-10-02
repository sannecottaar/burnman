# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2017 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
Solid and liquid endmembers and solid solution from Boukare et al., 2015
^^^^^^^^^^^^^^

"""
from __future__ import absolute_import

from .. import mineral_helpers as helpers
from ..mineral import Mineral
from ..solidsolution import SolidSolution
from ..solutionmodel import *
from ..processchemistry import dictionarize_formula, formula_mass

from .DKS_2013_liquids import MgO_liquid, SiO2_liquid


class mg_bridgmanite(Mineral):
    """
    MgSiO3 bridgmanite
    Boukare et al. (2015)
    """
    def __init__(self):
        formula = 'MgSiO3'
        formula = dictionarize_formula(formula)
        self.params = {
            'name': 'Mg bridgmanite',
            'formula': formula,
            'T_0': 300.,
            'V_0': 24.45e-6,
            'K_0': 251.e9,
            'Kprime_0': 4.14,
            'grueneisen_0': 1.57,
            'q_0': 1.1,
            'Debye_0': 905.,
            'F_0': -1408.e3,
            'equation_of_state': 'mgd3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula)}
        Mineral.__init__(self)

class fe_bridgmanite(Mineral):
    """
    FeSiO3 bridgmanite
    Boukare et al. (2015)
    """
    def __init__(self):
        formula = 'FeSiO3'
        formula = dictionarize_formula(formula)
        self.params = {
            'name': 'Fe bridgmanite',
            'formula': formula,
            'T_0': 300.,
            'V_0': 25.49e-6,
            'K_0': 272.e9,
            'Kprime_0': 4.1,
            'grueneisen_0': 1.57,
            'q_0': 1.1,
            'Debye_0': 871.,
            'F_0': -1048.e3,
            'equation_of_state': 'mgd3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula)}
        Mineral.__init__(self)

        
class periclase(Mineral):
    """
    Periclase
    Boukare et al. (2015)
    """
    def __init__(self):
        formula = 'MgO'
        formula = dictionarize_formula(formula)
        self.params = {
            'name': 'Periclase',
            'formula': formula,
            'T_0': 300.,
            'V_0': 11.24e-6,
            'K_0': 161.e9,
            'Kprime_0': 4.8,
            'grueneisen_0': 1.3,
            'q_0': 1.7,
            'Debye_0': 767.,
            'F_0': -569.e3,
            'equation_of_state': 'mgd3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula)}
        Mineral.__init__(self)


class wuestite(Mineral):
    """
    Wuestite
    Boukare et al. (2015)
    """
    def __init__(self):
        formula = 'FeO'
        formula = dictionarize_formula(formula)
        self.params = {
            'name': 'Wuestite',
            'formula': formula,
            'T_0': 300.,
            'V_0': 12.256e-6,
            'K_0': 149.e9,
            'Kprime_0': 3.6,
            'grueneisen_0': 1.41,
            'q_0': 0.5,
            'Debye_0': 417.,
            'F_0': -165.e3,
            'equation_of_state': 'mgd3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula)}
        Mineral.__init__(self)


class stishovite(Mineral):
    """
    Stishovite
    Boukare et al. (2015)
    """
    def __init__(self):
        formula = 'SiO2'
        formula = dictionarize_formula(formula)
        self.params = {
            'name': 'Stishovite',
            'formula': formula,
            'T_0': 300.,
            'V_0': 14.02e-6,
            'K_0': 314.e9,
            'Kprime_0': 3.8,
            'grueneisen_0': 1.37,
            'q_0': 2.8,
            'Debye_0': 1108.,
            'F_0': -819.e3,
            'equation_of_state': 'mgd3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula)}
        Mineral.__init__(self)


class ferropericlase(SolidSolution):
    def __init__(self, molar_fractions=None):
        self.name = 'magnesiowustite/ferropericlase'
        self.solution_type = 'ideal'
        self.endmembers = [[periclase(), '[Mg]O'], [wuestite(), '[Fe]O']]

        SolidSolution.__init__(self, molar_fractions=molar_fractions)

        
class bridgmanite(SolidSolution):
    def __init__(self, molar_fractions=None):
        self.name = 'Mg-Fe bridgmanite'
        self.solution_type = 'ideal'
        self.endmembers = [[mg_bridgmanite(), '[Mg]SiO3'], [fe_bridgmanite(), '[Fe]SiO3']]
        
        SolidSolution.__init__(self, molar_fractions=molar_fractions)



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

#print([E0_13 - Tref*S0_13 - Tref*Tref*dSdT_13/2. + Pref*V0_13,
#       E0_23 - Tref*S0_23 - Tref*Tref*dSdT_23/2. + Pref*V0_23])
#print([S0_13 + Tref*dSdT_13, S0_23 + Tref*dSdT_23])

class melt(SolidSolution):
    def __init__(self, molar_fractions=None):
        self.name = 'FMS melt'
        self.solution_type = 'boukare'
        self.endmembers = [[MgO_liquid(), '[Fe]O'],
                           [MgO_liquid(), '[Mg]O'],
                           [SiO2_liquid(), '[Si]O2']]
        self.energy_interaction = [[0., E0_23],
                                   [E0_13]]
        self.volume_interaction = [[0., V0_23],
                                   [V0_13]]
        self.entropy_interaction = [[0., S0_23],
                                    [S0_13]]
        self.dentropydT_interaction = [[0., dSdT_23],
                                       [dSdT_13]]
        self.alphas = [1., 1., 1.]
        SolidSolution.__init__(self, molar_fractions=molar_fractions)



