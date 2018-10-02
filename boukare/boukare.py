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
from burnman.equilibrate import equilibrate
from burnman.minerals.boukare import bridgmanite, ferropericlase, stishovite, melt



fper = ferropericlase()
bdg = bridgmanite()
stv = stishovite()
liq = melt()

fper.set_composition([0.9, 0.1])
fper.set_state(100.e9, 2000)


#bdg = burnman.minerals.SLB_2011.mg_fe_perovskite()
#fper = burnman.minerals.SLB_2011.ferropericlase()
# M_FeO = 0.071844 kg/mol

#bdg.guess = np.array([0.985, 0.015, 0.0])
bdg.guess = np.array([0.985, 0.015])
fper.guess = np.array([0.90, 0.10])


pressures = np.linspace(40.e9, 140.e9, 6)
for press in pressures:
    temperatures = np.linspace(2000., 4000., 21)
    composition = {'Fe': 0.2, 'Mg': 1.8, 'Si': 1.5, 'O': 5.}
    assemblage = burnman.Composite([bdg, fper])
    equality_constraints = [('P', press), ('T', temperatures)]
    sol, prm = equilibrate(composition, assemblage, equality_constraints, store_iterates=False)
    P, T, x_bdg, p_fbdg, x_per, p_wus = np.array([s.x for s in sol]).T

    plt.plot(T, p_fbdg*(1. - p_wus)/((1. - p_fbdg)*p_wus), label='{0} GPa'.format(press/1.e9))

plt.legend(loc='best')
plt.xlabel('Temperature (K)')
plt.ylabel('[p$_{FeSiO_3}$p$_{MgO}$]/[p$_{MgSiO_3}$p$_{FeO}$]')
plt.savefig('Fe_Mg_partitioning_Boukare_solids.pdf')
plt.show()

exit()
    


params = {'T_0': 2000.,
          'P_0': 20.e9, # ?
          'molar_mass': 0.071844,
          'rho_0': 3826.41, # apparently at 0 GPa, in order to fit Fig 1. from notes
          #'K_0': 30.961e9, # from notes
          'K_0': 30.196e9, # from paper, typo in notes? This value fits Fig 1 better
          #'rho_0': 0.071844/18.48e-6, # from paper
          'a_0': 9.54e-5, # +ve/-ve in Boukare paper/notes. +ve fits the density in Fig 1
          'q_0': -0.97, # from paper and notes
          'n': 3.26, # from paper and notes
          'S_0': 175.617,
          'C_0': 3.e-3,
          'C_1': 73.7532, # Cp = C1 + C0*T
          'G_0': -115.352e3} # at 20 GPa and 2000 K - with wustite, this is near melting point (Fig 2)


def gibbs(pressure, temperature, params):
    return (params['G_0'] - 
            params['S_0']*(temperature - params['T_0']) -
            params['C_1']*(temperature*np.log(temperature/params['T_0']) -
                           (temperature - params['T_0'])) -
            params['C_0']/2.*(temperature - params['T_0'])*(temperature - params['T_0']) -
            params['molar_mass']*params['a_0']*params['K_0']/(params['rho_0']*params['q_0']) *
            (np.power((1. + params['n']*(pressure/params['K_0'])),
                      -params['q_0']/params['n']) -
            np.power((1. + params['n']*params['P_0']/params['K_0']),
                     -params['q_0']/params['n'])) *
            (temperature - params['T_0']) +
            params['molar_mass']*params['K_0']/(params['rho_0']*(params['n'] - 1.)) *
            (np.power((1. + params['n']*pressure/params['K_0']),
                      1. - 1./params['n']) -
             np.power((1. + params['n']*params['P_0']/params['K_0']),
                      1. - 1./params['n'])))

def entropy(pressure, temperature, params):
    return (params['S_0'] +
            params['C_0']*(temperature - params['T_0']) + 
            params['C_1']*np.log(temperature/params['T_0']) +
            params['molar_mass']*params['a_0']*params['K_0']/(params['rho_0']*params['q_0']) *
            (np.power((1. + params['n']*pressure/params['K_0']), -params['q_0']/params['n']) -
             np.power((1. + params['n']*params['P_0']/params['K_0']), -params['q_0']/params['n'])))

def density(pressure, temperature, params):
    return params['rho_0']/(np.power(1. + params['n']*pressure/params['K_0'], -1./params['n']) +
                            np.power(1. + params['n']*pressure/params['K_0'], -(params['n'] + params['q_0'])/params['n'])*params['a_0']*(temperature - params['T_0']))


def volume(pressure, temperature, params):
    return ( (params['molar_mass'] / params['rho_0']) *
             (np.power(1. + params['n']*pressure/params['K_0'],
                       -1./params['n']) +
              np.power(1. + params['n']*pressure/params['K_0'],
                       -(params['n'] + params['q_0'])/params['n'])*params['a_0'] *
              (temperature - params['T_0'])) )


per = fper.endmembers[0][0]
wus = fper.endmembers[1][0]

FeO_melting_curve = np.loadtxt('boukare_melting_curves_FeO.dat')
MgO_melting_curve = np.loadtxt('boukare_melting_curves_MgO.dat')
SiO2_melting_curve = np.loadtxt('boukare_melting_curves_SiO2.dat')

boukare_densities = mpimg.imread('boukare_densities.png')

fig = plt.figure()
ax = [fig.add_subplot(1, 2, i) for i in range(1, 3)]
ax[0].imshow(boukare_densities, extent=[20., 140., 3.025, 9.0], aspect='auto')

pressures, temperatures = FeO_melting_curve.T
pressures *= 1.e9
ax[0].plot(pressures/1.e9, density(pressures, temperatures, params)/1000.,
           linestyle=':', linewidth=10., label='FeO liq')
ax[0].plot(pressures/1.e9, wus.evaluate(['density'], pressures, temperatures)[0]/1000.,
           linestyle=':', linewidth=10., label='FeO sol')

ax[1].plot(pressures, entropy(pressures, temperatures, params), label='liq')
ax[1].plot(pressures, wus.evaluate(['S'], pressures, temperatures)[0], label='sol')

pressures, temperatures = MgO_melting_curve.T
pressures *= 1.e9
ax[0].plot(pressures/1.e9, liq.endmembers[1][0].evaluate(['density'],
                                                         pressures, temperatures)[0]/1000.,
           linestyle=':', linewidth=10., label='MgO liq')
ax[0].plot(pressures/1.e9, per.evaluate(['density'], pressures, temperatures)[0]/1000.,
           linestyle=':', linewidth=10., label='MgO sol')

pressures, temperatures = SiO2_melting_curve.T
pressures *= 1.e9
ax[0].plot(pressures/1.e9, liq.endmembers[2][0].evaluate(['density'],
                                                         pressures, temperatures)[0]/1000.,
           linestyle=':', linewidth=10., label='SiO2 liq')
ax[0].plot(pressures/1.e9, stv.evaluate(['density'], pressures, temperatures)[0]/1000.,
           linestyle=':', linewidth=10., label='SiO2 sol')


ax[0].legend(loc='best')

plt.show()





from scipy.optimize import fsolve


boukare_melting = mpimg.imread('boukare_melting_crop.png')
plt.imshow(boukare_melting, extent=[20., 140., 2000, 8000], aspect='auto')

pressures = np.linspace(20.e9, 140.e9, 25)
temperatures = np.empty_like(pressures)


# FeO
def del_gibbs_T(T, P):
    return gibbs(P, T[0], params) - wus.evaluate(['gibbs'], P, T[0])[0]
guess = 2000.
for i, P in enumerate(pressures):
    temperatures[i] = fsolve(del_gibbs_T, [guess], args=(P))[0]
    guess = temperatures[i]

plt.plot(pressures/1.e9, temperatures, linestyle=':', linewidth=10.)



# SiO2
liq.set_composition([0.0, 0.0, 1.0])
print(liq.evaluate(['gibbs'], 20.e9, 3000.)[0])
print(stv.evaluate(['gibbs'], 20.e9, 3000.)[0])
per2 = burnman.minerals.DKS_2013_solids.periclase()
stv2 = burnman.minerals.DKS_2013_solids.stishovite()
stv3 = burnman.minerals.SLB_2011.stishovite()
stv2.set_state(20.e9, 3000)
stv3.set_state(20.e9, 3000)
print(stv2.gibbs)
print(stv3.gibbs)
def del_gibbs_MgO_T(T, P):
    liq.set_composition([0.0, 1.0, 0.0])
    return liq.evaluate(['gibbs'], P, T[0])[0] - per2.evaluate(['gibbs'], P, T[0])[0]

def del_gibbs_SiO2_T(T, P):
    liq.set_composition([0.0, 0.0, 1.0])
    return liq.evaluate(['gibbs'], P, T[0])[0] - stv2.evaluate(['gibbs'], P, T[0])[0]

guess = 3000.
for i, P in enumerate(pressures):
    temperatures[i] = fsolve(del_gibbs_MgO_T, [guess], args=(P))[0]
    guess = temperatures[i]

plt.plot(pressures/1.e9, temperatures, linestyle=':', linewidth=10.)

guess = 3000.
for i, P in enumerate(pressures):
    temperatures[i] = fsolve(del_gibbs_SiO2_T, [guess], args=(P))[0]
    guess = temperatures[i]

plt.plot(pressures/1.e9, temperatures, linestyle=':', linewidth=10.)


plt.show()
