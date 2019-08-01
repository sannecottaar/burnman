from __future__ import absolute_import

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1, os.path.abspath('..'))

import burnman
from burnman import Mineral
from burnman.processchemistry import dictionarize_formula, formula_mass




fo = burnman.minerals.HP_2011_ds62.fo()
hp_fo = burnman.minerals.HP_2011_ds62.fo()

Pref = 20.e9

hp_fo.params['Pref'] = Pref
hp_fo.params['Cp_Pref'] = [0., 0., 0., 0.]
del(hp_fo.params['H_0'])
del(hp_fo.params['S_0'])
del(hp_fo.params['Cp'])

fo.set_state(Pref, 298.15)
hp_fo.params['H_Pref'] = fo.H
hp_fo.params['S_Pref']: fo.S

hp_fo.set_method('mod_hp_tmt')
#hp_fo.params['Cp_Pref']: [0., 0., 0., 0.]


from scipy.optimize import curve_fit
func_Cp = lambda T, a, b, c, d: a + b*T + c/T/T + d/np.sqrt(T)
temperatures = np.linspace(298.15, 400., 101)
pressures = Pref + 0.*temperatures
print('New heat capacities at a reference pressure of {0} GPa:'.format(Pref/1.e9))
for m, m_mod in [[fo, hp_fo]]:
    heat_capacities = m.evaluate(['C_p'], pressures, temperatures)[0]

    guesses = m.params['Cp']
    popt, _ = curve_fit(func_Cp, temperatures, heat_capacities, guesses)
    m_mod.params['Cp_Pref'] = popt
    print('{0} {1}'.format(m.name, popt))
print('(Calibrated between {0} and {1} K)'.format(temperatures[0], temperatures[-1]))

temperatures = np.linspace(298.15, 2000, 101)
plt.plot(temperatures, fo.evaluate(['C_p'], pressures, temperatures)[0], linestyle=':')
plt.plot(temperatures, hp_fo.evaluate(['C_p'], pressures, temperatures)[0])

pressures = 0.*temperatures
plt.plot(temperatures, fo.evaluate(['C_p'], pressures, temperatures)[0], linestyle=':')
plt.plot(temperatures, hp_fo.evaluate(['C_p'], pressures, temperatures)[0])

plt.show()
