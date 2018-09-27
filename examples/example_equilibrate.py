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
from burnman.minerals import HP_2011_ds62, SLB_2011, JH_2015
from burnman import equilibrate

ordering = False # This example plots the state of order of the Jennings and Holland orthopyroxene in the simple en-fs binary at 1 bar.
aluminosilicates = False # The example plots the andalusite-sillimanite-kyanite phase diagram
gt_solvus = False # This example demonstrates the shape of the pyrope-grossular solvus
lower_mantle = True # This example calculates temperatures and assemblage properties along an isentrope in the lower mantle
upper_mantle = False # This example produces a 2D grid of the ol-opx-gt field
olivine_polymorphs = True # This example produces a P-T pseudosection for a fo90 composition 

gt = SLB_2011.garnet()
ol = SLB_2011.mg_fe_olivine()
wad = SLB_2011.mg_fe_wadsleyite()
rw = SLB_2011.mg_fe_ringwoodite()
bdg = SLB_2011.mg_fe_bridgmanite()
ppv = SLB_2011.post_perovskite()
per = SLB_2011.ferropericlase()
opx = SLB_2011.orthopyroxene()
stv = SLB_2011.stishovite()
coe = SLB_2011.coesite()
cpv = SLB_2011.ca_perovskite()

ol.guess = np.array([0.93, 0.07])
wad.guess = np.array([0.91, 0.09]) # 0.91 0.09 works for olivine polymorphs...
rw.guess = np.array([0.93, 0.07])
opx.guess = np.array([0.68, 0.08, 0.15, 0.09])
gt.guess = np.array([0.42, 0.12, 0.46, 0.0, 0.00])
bdg.guess = np.array([0.86, 0.1, 0.04]) # 
ppv.guess = np.array([0.86, 0.1, 0.01]) # bdg-in works if guess[2] = 0.
per.guess = np.array([0.9, 0.1])


if ordering:
    orthopyroxene = JH_2015.orthopyroxene()
    orthopyroxene.guess = np.array([1./3., 1./3., 1./3., 0., 0., 0.])
    temperatures = np.linspace(300., 2000., 41)

    Mg_numbers = np.linspace(10., 50., 5)
    assemblage = burnman.Composite([orthopyroxene])
    equality_constraints = [('P', 1.e5), ('T', temperatures)]
    
    for Mg_number in Mg_numbers:
        composition = {'Mg': Mg_number/100.*2., 'Fe': (1.-Mg_number/100.)*2., 'Si': 2., 'O': 6.}
        sols, prm = equilibrate(composition, assemblage, equality_constraints, store_iterates=False)
        Ts = np.array([sol.x[1] for sol in sols if sol.success])
        p_fms = np.array([sol.x[-1] for sol in sols if sol.success])
        plt.plot(Ts, p_fms, label='Mg# = {0}'.format(Mg_number))
    plt.xlabel("Temperature (K)")
    plt.ylabel("Proportion of ordered orthopyroxene")
    plt.legend(loc='best')
    plt.show()


    
    
if aluminosilicates:
    sillimanite = HP_2011_ds62.sill()
    andalusite = HP_2011_ds62.andalusite()
    kyanite = HP_2011_ds62.ky()
    
    composition = sillimanite.formula
    assemblage = burnman.Composite([sillimanite, andalusite, kyanite])
    equality_constraints = [('phase_proportion', (kyanite, np.array([0.0]))), ('phase_proportion', (sillimanite, np.array([0.0])))]
    sol, prm = equilibrate(composition, assemblage, equality_constraints, store_iterates=False)
    P_inv, T_inv = sol.x[0:2]
    print('invariant point found at {0} GPa, {1} K'.format(P_inv/1e9, T_inv))
    
    pressures = np.linspace(1.e5, P_inv, 21)
    assemblage = burnman.Composite([andalusite, kyanite])
    equality_constraints = [('P', pressures), ('phase_proportion', (andalusite, np.array([0.0])))]
    sols, prm = equilibrate(composition, assemblage, equality_constraints, store_iterates=False)
    Ps = np.array([sol.x[0] for sol in sols if sol.success])
    Ts = np.array([sol.x[1] for sol in sols if sol.success])
    plt.plot(Ts, Ps/1.e9)
    
    assemblage = burnman.Composite([sillimanite, andalusite])
    equality_constraints = [('P', pressures), ('phase_proportion', (sillimanite, np.array([0.0])))]
    sols, prm = equilibrate(composition, assemblage, equality_constraints, store_iterates=False)
    Ps = np.array([sol.x[0] for sol in sols if sol.success])
    Ts = np.array([sol.x[1] for sol in sols if sol.success])
    plt.plot(Ts, Ps/1.e9)

    
    pressures = np.linspace(P_inv, 1.e9, 21)
    assemblage = burnman.Composite([sillimanite, kyanite])
    equality_constraints = [('P', pressures), ('phase_proportion', (sillimanite, np.array([0.0])))]
    sols, prm = equilibrate(composition, assemblage, equality_constraints, store_iterates=False)
    Ps = np.array([sol.x[0] for sol in sols if sol.success])
    Ts = np.array([sol.x[1] for sol in sols if sol.success])
    plt.plot(Ts, Ps/1.e9)
    plt.xlabel('Temperature (K)')
    plt.ylabel('Pressure (GPa)')
    
    plt.show()

if gt_solvus:
    gt1 = SLB_2011.garnet()
    gt2 = SLB_2011.garnet()
    
    gt1.guess = np.array([0.99, 0.00, 0.01, 0.0, 0.00])
    gt2.guess = np.array([0.01, 0.00, 0.99, 0.0, 0.00])
    
    
    pressure = 1.e5
    
    composition = {'Mg': 1.5, 'Ca': 1.5, 'Al': 2., 'Si': 3.0, 'O':12.0}
    
    
    
    assemblage = burnman.Composite([gt1, gt2])
    equality_constraints = [('P', pressure), ('phase_proportion', (gt1, 0.))]
    sol, prm = equilibrate(composition, assemblage, equality_constraints, store_iterates=False)
    T_closure = assemblage.temperature
    
    temperatures = np.linspace(200., T_closure - 0.01, 21)
    
    assemblage = burnman.Composite([gt1, gt2])
    equality_constraints = [('P', pressure), ('T', temperatures)]
    sols, prm = equilibrate(composition, assemblage, equality_constraints, store_iterates=False)
    x1s = np.array([sol.x[3] for sol in sols if sol.code==0])
    x2s = np.array([sol.x[6] for sol in sols if sol.code==0])
    Ts = np.array([sol.x[1] for sol in sols if sol.code==0])
    plt.plot(x1s, Ts)
    plt.plot(x2s, Ts)
    plt.xlabel('Molar proportion of pyrope')
    plt.ylabel('Temperature (K)')
    plt.show()

if lower_mantle:
    P0 = 25.e9
    T0 = 1600.
    
    composition = {'Fe': 0.2, 'Mg': 2.0, 'Si': 1.9, 'Ca': 0.2, 'Al': 0.4, 'O':6.8}
    assemblage = burnman.Composite([bdg, per, cpv])
    equality_constraints = [('P', P0), ('T', T0)]
    sol, prm = equilibrate(composition, assemblage, equality_constraints, store_iterates=False)
    
    S = np.array([assemblage.molar_entropy*assemblage.n_moles])
    
    
    assemblage = burnman.Composite([bdg, per, ppv, cpv])
    equality_constraints = [('S', S), ('phase_proportion', (bdg, np.array([0.])))]
    sol, prm = equilibrate(composition, assemblage, equality_constraints,
                           store_iterates=False,
                           initial_state = sol.x[0:2])
    P_bdg_in = assemblage.pressure
    
    
    equality_constraints = [('S', S), ('phase_proportion', (ppv, np.array([0.])))]
    sol, prm = equilibrate(composition, assemblage, equality_constraints,
                           store_iterates=False,
                           initial_state = sol.x[0:2])
    P_ppv_in = assemblage.pressure
    T_ppv_in = assemblage.temperature

    pressures = np.linspace(P_ppv_in, P_bdg_in, 21)
    equality_constraints = [('P', pressures), ('S', S)]
    sols1, prm1 = equilibrate(composition, assemblage, equality_constraints,
                              initial_state = [P_ppv_in, T_ppv_in],
                              initial_composition_from_assemblage = True,
                              store_iterates=False)
    p1 = np.array([sol.x for sol in sols1 if sol.success]).T


    assemblage = burnman.Composite([bdg, per, cpv])
    pressures = np.linspace(25.e9, P_ppv_in, 21)
    equality_constraints = [('P', pressures), ('S', S)]
    sols2, prm2 = equilibrate(composition, assemblage, equality_constraints, store_iterates=False)
    
    p2 = np.array([sol.x for sol in sols2 if sol.success]).T
    
    assemblage = burnman.Composite([ppv, per, cpv])
    pressures = np.linspace(P_bdg_in, 140.e9, 21)
    equality_constraints = [('P', pressures), ('S', S)]
    sols3, prm3 = equilibrate(composition, assemblage, equality_constraints, store_iterates=False)
    
    p3 = np.array([sol.x for sol in sols3 if sol.success]).T

    fig = plt.figure()
    ax = [fig.add_subplot(2, 3, i) for i in range(1, 7)]
    phases = [bdg, per, cpv, ppv]
    colors = ['red', 'green', 'blue', 'orange', 'purple']
    for p, prm in [(p1, prm1), (p2, prm2), (p3, prm3)]:
        
        pressures, temperatures = p[[0, 1],:]
        ax[0].plot(pressures/1.e9, temperatures, color='black')
        ax[0].set_xlabel('Pressure (GPa)')
        ax[0].set_ylabel('Temperature (K)')
        for i, phase in enumerate(phases):
            try:
                idx = prm.parameter_names.index('x({0})'.format(phase.name))
                x_phase = p[idx,:]
                ax[i+1].plot(pressures/1.e9, x_phase, color=colors[i])
                ax[i+1].set_ylim(0,2)
                ax[i+1].set_xlim(0,140)
                ax[i+1].set_xlabel('Pressure (GPa)')
                ax[i+1].set_ylabel('x({0})'.format(phase.name))
            except:
                pass
        
        try:
            pv_idx = prm.parameter_names.index('x({0})'.format(phases[0].name))
            per_idx = prm.parameter_names.index('x({0})'.format(phases[1].name))
            KD = p[pv_idx+1,:]*(1. - p[per_idx+1,:])/((1. - p[pv_idx+1,:] - p[pv_idx+2,:])*p[per_idx+1,:])
            ax[5].plot(pressures/1.e9, KD, color='red', label='pv KD')
        except:
            pass
        
        try:
            ppv_idx = prm.parameter_names.index('x({0})'.format(phases[3].name))
            per_idx = prm.parameter_names.index('x({0})'.format(phases[1].name))
            KD = p[ppv_idx+1,:]*(1. - p[per_idx+1,:])/((1. - p[ppv_idx+1,:] - p[ppv_idx+2,:])*p[per_idx+1,:])
            ax[5].plot(pressures/1.e9, KD, color='blue', label='ppv KD')
        except:
            pass

        ax[5].set_ylim(0., 1.)
        ax[5].set_xlabel('Pressure (GPa)')
        ax[5].set_ylabel('[FeSiO3/MgSiO3]/[FeO/MgO]')
        
        
    plt.show()

if upper_mantle:
    temperatures = np.linspace(800., 1500., 8)
    pressures = np.linspace(1.e9, 14.e9, 11)
    
    composition = {'Fe': 0.2, 'Mg': 2.0, 'Si': 1.9, 'Ca': 0.2, 'Al': 0.4, 'O':6.8}
    assemblage = burnman.Composite([ol, opx, gt])
    #equality_constraints = [('phase_proportion', (opx, np.array([0.0]))), ('T', temperatures)]
    equality_constraints = [('P', pressures), ('T', temperatures)]
    sol, prm = equilibrate(composition, assemblage, equality_constraints,
                           store_iterates=False)
    

# This example produces a P-T pseudosection for a fo90 composition 
if olivine_polymorphs:
    
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    
    composition = {'Fe': 0.2, 'Mg': 1.8, 'Si': 1.0, 'O':4.0}
    assemblage = burnman.Composite([ol, wad, rw])
    equality_constraints = [('phase_proportion', (ol, 0.0)),
                            ('phase_proportion', (rw, 0.0))]
    
    sol, prm = equilibrate(composition, assemblage, equality_constraints, store_iterates=False)
    Pinv1, Tinv1 = sol.x[0:2]

    if sol.code != 0:
        raise Exception("Couldn't find ol-wad-rw invariant")
          
    assemblage = burnman.Composite([ol, wad, rw])
    equality_constraints = [('phase_proportion', (wad, 0.0)),
                            ('phase_proportion', (rw, 0.0))]
    sol, prm = equilibrate(composition, assemblage, equality_constraints, store_iterates=False)
    Pinv2, Tinv2 = sol.x[0:2]
    
    temperatures = np.linspace(573., Tinv1, 8)
    assemblage = burnman.Composite([ol, wad, rw])
    
    equality_constraints = [('T', temperatures), ('phase_proportion', (ol, 0.0))]
    sols, prm = equilibrate(composition, assemblage, equality_constraints, store_iterates=False)
    Ps = np.array([sol.x[0] for sol in sols if sol.success])
    Ts = np.array([sol.x[1] for sol in sols if sol.success])
    ax.plot(Ts, Ps/1.e9, color='k')

    temperatures = np.linspace(573., Tinv2, 8)
    equality_constraints = [('T', temperatures), ('phase_proportion', (wad, 0.0))]
    sols, prm = equilibrate(composition, assemblage, equality_constraints, store_iterates=False)
    Ps = np.array([sol.x[0] for sol in sols if sol.success])
    Ts = np.array([sol.x[1] for sol in sols if sol.success])
    ax.plot(Ts, Ps/1.e9, color='k')    
    
    temperatures = np.linspace(Tinv2, Tinv1, 8)
    equality_constraints = [('T', temperatures), ('phase_proportion', (rw, 0.0))]
    sols, prm = equilibrate(composition, assemblage, equality_constraints, store_iterates=False)
    Ps = np.array([sol.x[0] for sol in sols if sol.success])
    Ts = np.array([sol.x[1] for sol in sols if sol.success])
    ax.plot(Ts, Ps/1.e9, color='k')    

    temperatures = np.linspace(573., Tinv2, 8)
    assemblage = burnman.Composite([ol, rw])
    equality_constraints = [('T', temperatures), ('phase_proportion', (rw, 0.0))]
    sols, prm = equilibrate(composition, assemblage, equality_constraints, store_iterates=False)
    Ps = np.array([sol.x[0] for sol in sols if sol.success])
    Ts = np.array([sol.x[1] for sol in sols if sol.success])
    ax.plot(Ts, Ps/1.e9, color='k')
    
    temperatures = np.linspace(Tinv2, 1773., 8)
    assemblage = burnman.Composite([ol, wad])
    equality_constraints = [('T', temperatures), ('phase_proportion', (wad, 0.0))]
    sols, prm = equilibrate(composition, assemblage, equality_constraints, store_iterates=False)
    Ps = np.array([sol.x[0] for sol in sols if sol.success])
    Ts = np.array([sol.x[1] for sol in sols if sol.success])
    ax.plot(Ts, Ps/1.e9, color='k')

    
    temperatures = np.linspace(Tinv1, 1773., 8)
    equality_constraints = [('T', temperatures), ('phase_proportion', (wad, 1.0))]
    sols, prm = equilibrate(composition, assemblage, equality_constraints, store_iterates=False)
    Ps = np.array([sol.x[0] for sol in sols if sol.success])
    Ts = np.array([sol.x[1] for sol in sols if sol.success])
    ax.plot(Ts, Ps/1.e9, color='k')
    
    
    temperatures = np.linspace(Tinv1, 1773., 8)
    assemblage = burnman.Composite([wad, rw])
    equality_constraints = [('T', temperatures), ('phase_proportion', (rw, 0.0))]
    sols, prm = equilibrate(composition, assemblage, equality_constraints, store_iterates=False)
    Ps = np.array([sol.x[0] for sol in sols if sol.success])
    Ts = np.array([sol.x[1] for sol in sols if sol.success])
    ax.plot(Ts, Ps/1.e9, color='k')
    
    temperatures = np.linspace(573., 1773., 8)
    equality_constraints = [('T', temperatures), ('phase_proportion', (rw, 1.0))]
    sols, prm = equilibrate(composition, assemblage, equality_constraints, store_iterates=False)
    Ps = np.array([sol.x[0] for sol in sols if sol.success])
    Ts = np.array([sol.x[1] for sol in sols if sol.success])
    ax.plot(Ts, Ps/1.e9, color='k')
    
    
    assemblage = burnman.Composite([rw, bdg, per])
    equality_constraints = [('T', temperatures), ('phase_proportion', (rw, 0.0))]
    sols, prm = equilibrate(composition, assemblage, equality_constraints, store_iterates=False)
    Ps = np.array([sol.x[0] for sol in sols if sol.success])
    Ts = np.array([sol.x[1] for sol in sols if sol.success])
    ax.plot(Ts, Ps/1.e9, color='k')
    
    equality_constraints = [('T', temperatures), ('phase_proportion', (rw, 1.0))]
    sols, prm = equilibrate(composition, assemblage, equality_constraints, store_iterates=False)
    Ps = np.array([sol.x[0] for sol in sols if sol.success])
    Ts = np.array([sol.x[1] for sol in sols if sol.success])
    ax.plot(Ts, Ps/1.e9, color='k')

    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel('Pressure (GPa)')
    plt.show()

    
