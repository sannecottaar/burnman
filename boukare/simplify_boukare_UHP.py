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

from burnman.minerals.boukare import wuestite_boukare, periclase_boukare, stishovite_boukare, mg_bridgmanite_boukare, fe_bridgmanite_boukare
from burnman.minerals.boukare import FeO_liquid_boukare, MgO_liquid_boukare, SiO2_liquid_boukare
from burnman.minerals.DKS_2013_solids import periclase, stishovite
from burnman.minerals.DKS_2013_liquids import MgO_liquid, SiO2_liquid
from burnman.minerals.SLB_2011 import seifertite


# Boukare liquids and solids
wus = wuestite_boukare()
per = periclase_boukare()
stv = stishovite_boukare()
mpv = mg_bridgmanite_boukare()
fpv = fe_bridgmanite_boukare()

liq_FeO = FeO_liquid_boukare()
liq_MgO = MgO_liquid_boukare()
liq_SiO2 = SiO2_liquid_boukare()

# BOUKARE MODEL
fig = plt.figure()
ax = [fig.add_subplot(2, 2, i) for i in range(1, 5)]


P_list = np.linspace(50.e9, 150.e9, 11)
T_list = np.linspace(2500., 4500., 5)

for i, (fper, stv, bdg) in enumerate([(per, stv, mpv),
                                      (wus, stv, fpv)]):
    print(i)
    
    
    for T in T_list:
        pressures = P_list
        temperatures = T + pressures*0.
        Vfper, Sfper = fper.evaluate(['V', 'S'], pressures, temperatures)
        Vstv, Sstv = stv.evaluate(['V', 'S'], pressures, temperatures)
        Vbdg, Sbdg = bdg.evaluate(['V', 'S'], pressures, temperatures)

        
        ax[i].set_title(bdg.name)
        ax[i].set_ylabel('V$_{melt}$ (/f.u.)')
        ax[i].plot(pressures/1.e9, Vbdg-(Vfper + Vstv), label='{0} K'.format(T))
        
    for P in P_list:
        temperatures = T_list
        pressures = P + temperatures*0.
        Vfper, Sfper = fper.evaluate(['V', 'S'], pressures, temperatures)
        Vstv, Sstv = stv.evaluate(['V', 'S'], pressures, temperatures)
        Vbdg, Sbdg = bdg.evaluate(['V', 'S'], pressures, temperatures)
                                     
        ax[i + 2].set_title(bdg.name)
        ax[i + 2].set_ylabel('S$_{melt}$ (/atom)')
        ax[i + 2].plot(temperatures, (Sbdg-(Sfper + Sstv))/bdg.params['n'], label='{0} GPa'.format(P/1.e9))
                                     

for i in range(0, 4):
    ax[i].legend()

plt.show()

        
    

# BOUKARE MODEL
fig = plt.figure()
ax = [fig.add_subplot(2, 3, i) for i in range(1, 7)]


P_list = np.linspace(50.e9, 150.e9, 11)
T_list = np.linspace(2500., 4500., 5)

for i, (solid, liquid) in enumerate([(per, liq_MgO),
                                     (stv, liq_SiO2),
                                     (wus, liq_FeO)]):
    print(i)
    
    
    for T in T_list:
        pressures = P_list
        temperatures = T + pressures*0.
        Vs, Ss = solid.evaluate(['V', 'S'], pressures, temperatures)
        Vl, Sl = liquid.evaluate(['V', 'S'], pressures, temperatures)

        
        ax[i].set_title(solid.name)
        ax[i].set_ylim(-5e-7, 5e-7)
        ax[i].set_ylabel('V$_{melt}$ (/f.u.)')
        ax[i].plot(pressures/1.e9, Vl-Vs, label='{0} K'.format(T))
        #ax[i + 3].plot(pressures/1.e9, Sl-Ss, label='{0} K'.format(T))
        
    for P in P_list:
        temperatures = T_list
        pressures = P + temperatures*0.
        Vs, Ss = solid.evaluate(['V', 'S'], pressures, temperatures)
        Vl, Sl = liquid.evaluate(['V', 'S'], pressures, temperatures)

        ax[i + 3].set_title(solid.name)
        ax[i + 3].set_ylim(0, 20)
        ax[i + 3].set_ylabel('S$_{melt}$ (/atom)')
        #ax[i].plot(temperatures, Vl-Vs, label='{0} GPa'.format(P/1.e9))
        ax[i + 3].plot(temperatures, (Sl-Ss)/solid.params['n'], label='{0} GPa'.format(P/1.e9))

for i in range(0, 6):
    ax[i].legend()




# DKS LIQUIDS AND SOLIDS
per = periclase()
stv = stishovite()
seif = seifertite()

liq_MgO = MgO_liquid()
liq_SiO2 = SiO2_liquid()


fig2 = plt.figure()
ax = [fig2.add_subplot(2, 3, i) for i in range(1, 7)]

for i, (solid, liquid) in enumerate([(per, liq_MgO),
                                     (stv, liq_SiO2),
                                     (wus, liq_FeO)]):
    print(i)
    
    for T in T_list:
        pressures = P_list
        temperatures = T + pressures*0.
        Vs, Ss = solid.evaluate(['V', 'S'], pressures, temperatures)
        Vl, Sl = liquid.evaluate(['V', 'S'], pressures, temperatures)

        
        ax[i].set_title(solid.name)
        ax[i].set_ylim(-5e-7, 5e-7)
        ax[i].set_ylabel('V$_{melt}$ (/f.u.)')
        ax[i].plot(pressures/1.e9, Vl-Vs, label='{0} K'.format(T))
        #ax[i + 3].plot(pressures/1.e9, Sl-Ss, label='{0} K'.format(T))
        
    for P in P_list:
        temperatures = T_list
        pressures = P + temperatures*0.
        Vs, Ss = solid.evaluate(['V', 'S'], pressures, temperatures)
        Vl, Sl = liquid.evaluate(['V', 'S'], pressures, temperatures)

        ax[i + 3].set_title(solid.name)
        ax[i + 3].set_ylim(0, 20)
        ax[i + 3].set_ylabel('S$_{melt}$ (/atom)')
        #ax[i].plot(temperatures, Vl-Vs, label='{0} GPa'.format(P/1.e9))
        ax[i + 3].plot(temperatures, (Sl-Ss)/solid.params['n'], label='{0} GPa'.format(P/1.e9))

for i in range(0, 6):
    ax[i].legend()

plt.show()

