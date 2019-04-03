from __future__ import absolute_import
import numpy as np
import matplotlib.pyplot as plt
from eos_and_averaging import thermodynamic_properties
from model_parameters import *

def liq_sol_molar_compositions(P, T,
                               melting_reference_pressure,
                               melting_temperatures,
                               melting_entropies,
                               melting_volumes,
                               n_mole_mix):
    
    dG = ((melting_temperatures - T) * melting_entropies +
          (P - melting_reference_pressure) * melting_volumes) # this is a linearised version (we don't use the melt and solid endmember equations of state directly).
    
    Xls = 1.0 - ((1.0 - np.exp(dG[1]/(n_mole_mix[1]*gas_constant*T))) /
                 (np.exp(dG[0]/(n_mole_mix[0]*gas_constant*T)) -
                  np.exp(dG[1]/(n_mole_mix[1]*gas_constant*T))))
           
    Xss = Xls * np.exp(dG[1]/(n_mole_mix[1]*gas_constant*T))

    return Xls, Xss

def calculate_xfe_in_fper_and_pv(pressure, temperature,
                                 x_Fe_in_solid, c_mantle):
    """
    Calculates the proportion of wus in fper, fpv in pv, and 
    the MASS fraction of pv in the solid given
    the Fe number in the solid.
    """
    x_MgO = (1. - x_Fe_in_solid)*c_mantle[0]['MgO']
    x_FeO = x_Fe_in_solid*c_mantle[1]['FeO']
    x_SiO2 = (1. - x_Fe_in_solid)*c_mantle[0]['SiO2'] + x_Fe_in_solid*c_mantle[1]['SiO2']
    
    norm = x_MgO + x_FeO
    f_FeO = x_FeO/norm
    f_SiO2 = x_SiO2/norm

    mpv_gibbs = thermodynamic_properties(pressure, temperature, mpv_params)['gibbs']
    fpv_gibbs = thermodynamic_properties(pressure, temperature, fpv_params)['gibbs']
    per_gibbs = thermodynamic_properties(pressure, temperature, per_params)['gibbs']
    wus_gibbs = thermodynamic_properties(pressure, temperature, wus_params)['gibbs']

    gibbs_rxn = mpv_gibbs + wus_gibbs - fpv_gibbs - per_gibbs
    KD = np.exp(gibbs_rxn/(gas_constant*temperature))

    # Solving equation 6 in Nakajima et al., 2012 for X_Fe_fp and X_Fe_pv
    # Solved using the definition of the distribution coefficient to define X_Fe_fp as a function of X_Fe_pv

    num_to_sqrt = ((-4. * f_FeO * (KD - 1.) * KD * f_SiO2) +
                   (np.power(1. + (f_FeO * (KD - 1)) + ((KD - 1.) * f_SiO2), 2.)))

    p_fpv = ((-1. + f_FeO - (f_FeO * KD) + f_SiO2 - (f_SiO2 * KD) + np.sqrt(num_to_sqrt)) /
               (2. * f_SiO2 * (1. - KD)))

    p_wus = p_fpv / (((1. - p_fpv) * KD) + p_fpv)
    
    f_pv = f_SiO2 # MOLAR mass of bridgmanite in the solid

    # finally, convert to mass fraction of bridgmanite in the solid
    molar_mass_fper = p_wus*wus_params['molar_mass'] + (1. - p_wus)*per_params['molar_mass']
    molar_mass_pv = p_fpv*fpv_params['molar_mass'] + (1. - p_fpv)*mpv_params['molar_mass']

    mass_fraction_pv = f_pv*molar_mass_pv/(f_pv*molar_mass_pv + (1. - f_pv)*molar_mass_fper)
    
    return p_wus, p_fpv, mass_fraction_pv # f_pv is the molar fraction of pv


def calculate_xfe_in_solid_from_molar_proportions(p_wus, p_fpv, mass_fraction_pv):
    """
    Calculates the Fe number in the solid given 
    the proportion of wus in fper, fpv in pv, and 
    the molar fraction of pv in the solid.
    """

    molar_mass_fper = p_wus*wus_params['molar_mass'] + (1. - p_wus)*per_params['molar_mass']
    molar_mass_pv = p_fpv*fpv_params['molar_mass'] + (1. - p_fpv)*mpv_params['molar_mass']

    f_pv = ((mass_fraction_pv/molar_mass_pv)/
            ((mass_fraction_pv/molar_mass_pv) + (1. - mass_fraction_pv)/molar_mass_fper)) # molar fraction of bridgmanite
    
    return f_pv*p_fpv + (1. - f_pv)*p_wus




################# BEGIN PLOTS ######################
plt.rc('font', family='DejaVu sans', size=15.)


# 1) Plot phase proportions as a function of the 1D compositional parameter x_Fe
x_Fes = np.linspace(1.e-12, 1.-1.e-12, 101)
p_wus = np.empty_like(x_Fes)
p_fpv = np.empty_like(x_Fes)
mass_f_pv = np.empty_like(x_Fes)

pressure = 100.e9
temperature = 3000.
for i, x_Fe_in_solid in enumerate(x_Fes):

    p_wus[i], p_fpv[i], mass_f_pv[i] = calculate_xfe_in_fper_and_pv(pressure, temperature,
                                                                    x_Fe_in_solid, c_mantle)

plt.plot(x_Fes, p_wus, label='molar proportion wus in fper')
plt.plot(x_Fes, p_fpv, label='molar proportion fpv in pv')
plt.plot(x_Fes, mass_f_pv, label='mass fraction pv', linestyle=':')

plt.xlabel('$x_{Fe}$')
plt.ylabel('Proportions (fractional)')
plt.legend()
plt.show()


# 2) Plot melt curves and melting entropies/volumes as a function of the 1D compositional parameter x_Fe
fig = plt.figure()
fig.set_size_inches(18.0, 12.0)

ax = [fig.add_subplot(2, 2, i) for i in range(1, 5)]
for P in [120.e9, 130.e9, 140.e9]:
    dT = melting_volumes/melting_entropies*(P - melting_reference_pressure) # dT/dP = DV/DS
    T_melt = melting_temperatures + dT
    temperatures = np.linspace(T_melt[0], T_melt[1], 101)

    
    Xls = np.empty_like(temperatures)
    Xss = np.empty_like(temperatures)
    for i, T in enumerate(temperatures):
        Xls[i], Xss[i] = liq_sol_molar_compositions(P, T,
                                                    melting_reference_pressure,
                                                    melting_temperatures,
                                                    melting_entropies,
                                                    melting_volumes,
                                                    n_mole_mix)


    ax[0].plot(Xls, temperatures, label='liquid composition, {0} GPa'.format(P/1.e9))
    ax[0].plot(Xss, temperatures, label='solid composition, {0} GPa'.format(P/1.e9))
ax[0].legend()

xs = np.linspace(0., 1., 101)
ax[1].plot(xs, melting_entropies.dot(np.array([1 - xs, xs])))
ax[2].plot(xs, melting_volumes.dot(np.array([1 - xs, xs]))*1.e6)

for i in range(4):
    ax[i].set_xlim(0., 1.)
    ax[i].set_xlabel('$x_{Fe}$')

ax[0].set_ylabel('Temperature (K)')
ax[1].set_ylabel('$\Delta S_{melting}$ (J/K/mol_cations)')
ax[2].set_ylabel('$\Delta V_{melting}$ (cm$^3$/mol_cations)')
    
plt.show()

