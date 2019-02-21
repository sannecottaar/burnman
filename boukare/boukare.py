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
from burnman.minerals.boukare import bridgmanite_boukare, ferropericlase_boukare, stishovite_boukare, melt_boukare


fper = ferropericlase_boukare()
bdg = bridgmanite_boukare()
stv = stishovite_boukare()
liq = melt_boukare()

mpv = bdg.endmembers[0][0]
fpv = bdg.endmembers[1][0]

per = fper.endmembers[0][0]
wus = fper.endmembers[1][0]

liq_FeO = liq.endmembers[0][0]
liq_MgO = liq.endmembers[1][0]
liq_SiO2 = liq.endmembers[2][0]

for m in [liq_MgO, liq_SiO2]:
    if not burnman.tools.check_eos_consistency(m, 10.e9, 3000.):
        burnman.tools.check_eos_consistency(m, 10.e9, 3000., verbose=True)
        raise Exception('{0} EoS invalid'.format(m.name))


"""
liq_SiO2.params['S_0'] -= 15.
T = 3650.
liq_SiO2.set_state(20.e9, T)
stv.set_state(20.e9, T)
print(stv.gibbs, liq_SiO2.gibbs)

liq_SiO2.params['F_0'] += stv.gibbs - liq_SiO2.gibbs
print(stv.gibbs - liq_SiO2.gibbs)
"""

props = []
for m in [liq_MgO, liq_SiO2, liq_FeO]:
    m.set_state(1.e5, 3000.)
    props.append([m.V*1.e6, m.K_T/1.e9, m.alpha*1.e5, m.C_p/m.params['n']/burnman.constants.gas_constant])
props = np.array(props)
n = ['V:  ', 'K_T: ', 'alpha:', 'Cp/nR:']
for i, p in enumerate(props.T):
    print('{0} {1:0.2f}, {2:0.2f}, {3:0.2f}'.format(n[i], *p))

p = [1447.2, -0.24865, 10.27e6, -1.1258]
P_st = lambda p, x: p[0]*np.exp(p[1]*x) + p[2]*np.exp(p[3]*x)
print(P_st(p, 27.4))


fig = plt.figure()
ax = [fig.add_subplot(2, 2, i) for i in range(1, 5)]

boukare_sio2_volumes = mpimg.imread('boukare_sio2_volumes.png')
boukare_sio2_energies = mpimg.imread('boukare_sio2_energies.png')
boukare_sio2_akts = mpimg.imread('boukare_sio2_akts.png')
boukare_sio2_cvs = mpimg.imread('boukare_sio2_cvs.png')

ax[0].imshow(boukare_sio2_volumes, extent=[10., 30., -0.3, 300], aspect='auto')
ax[1].imshow(boukare_sio2_energies, extent=[10., 30., -2400., -500.], aspect='auto')
ax[2].imshow(boukare_sio2_cvs, extent=[7., 33., 120., 170.], aspect='auto')
ax[3].imshow(boukare_sio2_akts, extent=[10., 30., 0., 0.02], aspect='auto')

volumes = np.linspace(10., 30., 101)
ax[0].plot(volumes, P_st(p, volumes), linestyle=':', linewidth=8.)

stv2 = burnman.minerals.SLB_2011.stishovite()

pressures = np.linspace(5.e9, 200.e9, 101)
for T in np.linspace(2000., 7000., 6):
    temperatures = T + pressures*0.

    
    G, E, V, K_T, alpha, C_v = liq_SiO2.evaluate(['gibbs', 'molar_internal_energy', 'V', 'K_T', 'alpha', 'C_v'], pressures, temperatures)
    ax[0].plot(V*1.e6, pressures/1.e9, linestyle='--', linewidth=2.)
    ax[1].plot(V*1.e6, E/1000. - 1640., linestyle='--', linewidth=2.)
    ax[2].plot(V*1.e6, C_v, linestyle='--', linewidth=2.)
    ax[2].plot([0., 30.], [104.51, 104.51+30.e3*0.1353e-2], linestyle=':')
    ax[2].plot([0., 30.], [104.51+0.6e-4*10000.,
                           104.51+30.e3*0.1353e-2 + 0.6e-4*10000.], linestyle=':')
    ax[3].plot(V*1.e6, alpha*K_T/1.e9, linestyle='--', linewidth=2.)

ax[0].scatter(liq_SiO2.params['V_0']*1.e6, 1.e-4)

plt.show()





fper.set_composition([0.9, 0.1])
fper.set_state(100.e9, 2000)

"""
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
"""

# FeO-SiO2 phase diagram
P = 130.e9
temperatures = np.linspace(3400., 3700., 7)
composition = {'Fe': 0.99, 'Si': 0.01, 'O': 1.01}
liq.guess = np.array([0.99, 0., 0.01])

assemblage = burnman.Composite([wus, liq])
assemblage.set_state(P, 3000.)
equality_constraints = [('P', P), ('T', temperatures)]
sols, prm = equilibrate(composition, assemblage,
                       equality_constraints,
                       initial_state_from_assemblage=True,
                       store_iterates=False)
T, Si_melt = np.array([[sol.x[1], sol.x[-1]] for sol in sols]).T

plt.plot(Si_melt, T, label='{0} GPa'.format(P/1.e9))

composition = {'Fe': 0.01, 'Si': 0.99, 'O': 1.99}
liq.guess = np.array([0.01, 0., 0.99])

temperatures = np.linspace(3400., 5200., 7)
assemblage = burnman.Composite([stv, liq])
assemblage.set_state(P, 3000.)
equality_constraints = [('P', P), ('T', temperatures)]
sols, prm = equilibrate(composition, assemblage,
                       equality_constraints,
                       initial_state_from_assemblage=True,
                       store_iterates=False)
T, Si_melt = np.array([[sol.x[1], sol.x[-1]] for sol in sols]).T

plt.plot(Si_melt, T, label='{0} GPa'.format(P/1.e9))
plt.xlim(0., 1.)
plt.show()



# MgO-SiO2 phase diagram
P = 130.e9
temperatures = np.linspace(4800., 7400., 7)
composition = {'Mg': 0.99, 'Si': 0.01, 'O': 1.01}
liq.guess = np.array([0.0, 0.99, 0.01])

assemblage = burnman.Composite([per, liq])
assemblage.set_state(P, 3000.)
equality_constraints = [('P', P), ('T', temperatures)]
sols, prm = equilibrate(composition, assemblage,
                       equality_constraints,
                       initial_state_from_assemblage=True,
                       store_iterates=False)
T, Si_melt = np.array([[sol.x[1], sol.x[-1]] for sol in sols]).T

plt.plot(Si_melt, T, label='{0} GPa'.format(P/1.e9))

composition = {'Mg': 0.01, 'Si': 0.99, 'O': 1.99}
liq.guess = np.array([0., 0.01, 0.99])

temperatures = np.linspace(4400., 5200., 7)
assemblage = burnman.Composite([stv, liq])
assemblage.set_state(P, 3000.)
equality_constraints = [('P', P), ('T', temperatures)]
sols, prm = equilibrate(composition, assemblage,
                       equality_constraints,
                       initial_state_from_assemblage=True,
                       store_iterates=False)
T, Si_melt = np.array([[sol.x[1], sol.x[-1]] for sol in sols]).T

plt.plot(Si_melt, T, label='{0} GPa'.format(P/1.e9))


composition = {'Mg': 0.49, 'Si': 0.51, 'O': 1.51}
liq.guess = np.array([0., 0.49, 0.51])

temperatures = np.linspace(4400., 4900., 7)
assemblage = burnman.Composite([mpv, liq])
assemblage.set_state(P, 3000.)
equality_constraints = [('P', P), ('T', temperatures)]
sols, prm = equilibrate(composition, assemblage,
                       equality_constraints,
                       initial_state_from_assemblage=True,
                       store_iterates=False)
T, Si_melt = np.array([[sol.x[1], sol.x[-1]] for sol in sols]).T

plt.plot(Si_melt, T, label='{0} GPa'.format(P/1.e9))


composition = {'Mg': 0.51, 'Si': 0.49, 'O': 1.49}
liq.guess = np.array([0., 0.51, 0.49])

temperatures = np.linspace(4400., 4900., 7)
assemblage = burnman.Composite([mpv, liq])
assemblage.set_state(P, 3000.)
equality_constraints = [('P', P), ('T', temperatures)]
sols, prm = equilibrate(composition, assemblage,
                       equality_constraints,
                       initial_state_from_assemblage=True,
                       store_iterates=False)
T, Si_melt = np.array([[sol.x[1], sol.x[-1]] for sol in sols]).T

plt.plot(Si_melt, T, label='{0} GPa'.format(P/1.e9))
plt.xlim(0., 1.)
plt.show()


exit()

# pressures, XFeOs
pressures = np.linspace(100.e9, 140.e9, 6)
x_Fe = 0.4
for press in pressures:
    composition = {'Fe': 2.*x_Fe, 'Mg': 2.*(1. - x_Fe), 'Si': 1.5, 'O': 5.}
    
    bdg.guess = np.array([1. - x_Fe, x_Fe])
    fper.guess = np.array([1. - x_Fe, x_Fe])
    x_Si = 0.5
    liq.guess = np.array([x_Fe*2.*(1 - x_Si), (1. - x_Fe/2.)*(1 - x_Si), x_Si])
    assemblage = burnman.Composite([bdg, fper, liq])
    assemblage.set_state(press, 3500.)
    equality_constraints = [('P', press), ('phase_proportion', (liq, np.array([0.])))]
    sol, prm = equilibrate(composition, assemblage, equality_constraints, initial_state_from_assemblage=True, store_iterates=False)
    #P, T, x_bdg, p_fbdg, x_per, p_wus = np.array([s.x for s in sol]).T
    print(assemblage)
    #plt.plot(T, p_fbdg*(1. - p_wus)/((1. - p_fbdg)*p_wus), label='{0} GPa'.format(press/1.e9))

exit()



FeO_melting_curve = np.loadtxt('boukare_melting_curves_FeO.dat')
MgO_melting_curve = np.loadtxt('boukare_melting_curves_MgO.dat')
SiO2_melting_curve = np.loadtxt('boukare_melting_curves_SiO2.dat')

boukare_densities = mpimg.imread('boukare_densities.png')

fig = plt.figure()
ax = [fig.add_subplot(1, 2, i) for i in range(1, 3)]
ax[0].imshow(boukare_densities, extent=[20., 140., 3.025, 9.0], aspect='auto')

pressures, temperatures = FeO_melting_curve.T
pressures *= 1.e9
ax[0].plot(pressures/1.e9, liq.endmembers[0][0].evaluate(['density'],
                                                         pressures, temperatures)[0]/1000.,
           linestyle=':', linewidth=10., label='FeO liq')
ax[0].plot(pressures/1.e9, wus.evaluate(['density'], pressures, temperatures)[0]/1000.,
           linestyle=':', linewidth=10., label='FeO sol')

ax[1].plot(pressures, liq.endmembers[0][0].evaluate(['S'],
                                                    pressures, temperatures)[0], label=' FeO liq')
ax[1].plot(pressures, wus.evaluate(['S'], pressures, temperatures)[0], label='FeO sol')

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


def affinity_2_mins(T, P, m1, m2):
    return m1.evaluate(['gibbs'], P, T[0])[0] - m2.evaluate(['gibbs'], P, T[0])[0]

guess = 2000.
for i, P in enumerate(pressures):
    temperatures[i] = fsolve(affinity_2_mins, [guess], args=(P, liq_FeO, wus))[0]
    guess = temperatures[i]

plt.plot(pressures/1.e9, temperatures, linestyle=':', linewidth=10.)


guess = 3000.
for i, P in enumerate(pressures):
    temperatures[i] = fsolve(affinity_2_mins, [guess], args=(P, liq_MgO, per))[0]
    guess = temperatures[i]

plt.plot(pressures/1.e9, temperatures, linestyle=':', linewidth=10.)


guess = 3000.
for i, P in enumerate(pressures):
    temperatures[i] = fsolve(affinity_2_mins, [guess], args=(P, liq_SiO2, stv))[0]
    guess = temperatures[i]

plt.plot(pressures/1.e9, temperatures, linestyle=':', linewidth=10.)


plt.show()
