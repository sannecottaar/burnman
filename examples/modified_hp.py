# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.

'''
example_equilibrate
--------------------

This example demonstrates how burnman may be used to calculate the
equilibrium phase proportions and compositions for an assemblage
of a fixed bulk composition.

*Uses:*

* :doc:`mineral_database`
* :class:`burnman.composite.Composite`
* :func:`burnman.equilibrate.equilibrate`
'''
from __future__ import absolute_import
from __future__ import print_function

import os
import sys
# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1, os.path.abspath('..'))

import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ggplot')

import burnman
from burnman import Mineral



def convert_to_hpht(m):
    from scipy.optimize import curve_fit
    temperatures = np.linspace(300., 2000., 101)
    Pref =  10.e9
    Cps = m.evaluate(['C_p'], temperatures*0. + Pref, temperatures)[0]
    m.set_state(Pref, m.params['T_0'])

    func_Cp = lambda T, a, b, c, d: a + b*T + c/T/T + d/np.sqrt(T)
    popt, pcov = curve_fit(func_Cp, temperatures, Cps)

    params = dict(m.params)
    params['equation_of_state'] = 'mod_hp_tmt'
    params['Pref'] = Pref
    params['H_Pref'] = m.H
    params['S_Pref'] = m.S
    params['Cp_Pref'] = popt

    del params['P_0']
    del params['H_0']
    del params['S_0']
    del params['Cp']

    return Mineral(params)
    


per1 = burnman.minerals.HP_2011_ds62.per()
per2 = convert_to_hpht(per1)


per1.set_state(1.e10, 6000.)
per2.set_state(1.e10, 6000.)

print(per1.gibbs, per2.gibbs)
        
