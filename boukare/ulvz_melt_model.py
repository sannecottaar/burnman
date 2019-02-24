from __future__ import absolute_import

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1, os.path.abspath('..'))


from modified_HP import thermodynamic_properties

gas_constant = 8.31446

#################### BEGIN EOS PARAMETERS #######################
mpv_params={'name': 'mpv',
            'a_0': 1.87e-05,
            'K_0': 251.0e9,
            'Pref': 100.0e9,
            'Cp_Pref': np.array([ 1.61546581e+02, -3.31714290e-03, -3.57533814e+06, -1.11254791e+03]),
            'H_Pref': -1442310.0,
            'Kprime_0': 4.14,
            'T_0': 298.15,
            'T_einstein': 561.0,
            'Kdprime_0': -1.6e-11,
            'V_0': 2.445e-05,
            'molar_mass': 0.1003887,
            'S_Pref': 62.6,
            'equation_of_state': 'mod_hp_tmt',
            'n': 5.0,
            'formula': {'Mg': 1.0, 'Si': 1.0, 'O': 3.0}}


fpv_params={'a_0': 1.87e-05,
            'K_0': 281.0e9,
            'Pref': 100.0e9,
            'Cp_Pref': np.array([ 1.39546209e+02,  6.36191292e-03, -4.13886524e+06, -4.64775577e+02]),
            'H_Pref': -1082910.0,
            'Kprime_0': 4.14,
            'T_0': 298.15,
            'T_einstein': 418.1,
            'Kdprime_0': -1.6e-11,
            'V_0': 2.534e-05,
            'name': 'fpv',
            'molar_mass': 0.1319287,
            'S_Pref': 95.0,
            'equation_of_state': 'mod_hp_tmt',
            'n': 5.0,
            'formula': {'Si': 1.0, 'Fe': 1.0, 'O': 3.0}}


per_params={'a_0': 3.11e-05,
            'K_0': 161.6e9,
            'Pref': 100.0e9,
            'Cp_Pref': np.array([ 7.31147154e+01, -6.35318887e-03, -7.33679285e+05, -5.92994207e+02]),
            'H_Pref': -601570.0,
            'Kprime_0': 3.95,
            'T_0': 298.15,
            'T_einstein': 540.2,
            'Kdprime_0': -2.4e-11,
            'V_0': 1.125e-05,
            'name': 'per',
            'molar_mass': 0.0403044,
            'S_Pref': 26.5,
            'equation_of_state': 'mod_hp_tmt',
            'n': 2.0,
            'formula': {'Mg': 1.0, 'O': 1.0}}


wus_params={'a_0': 3.22e-05,
            'K_0': 152.0e9,
            'Pref': 100.0e9,
            'Cp_Pref': np.array([ 5.20016403e+01, 3.36163516e-03, -1.19540964e+06,  2.55067110e+01]),
            'H_Pref': -262240.0,
            'Kprime_0': 4.9,
            'T_0': 298.15,
            'T_einstein': 297.6,
            'Kdprime_0': -3.2e-11,
            'V_0': 1.206e-05,
            'name': 'fper',
            'molar_mass': 0.0718444,
            'S_Pref': 58.6,
            'equation_of_state': 'mod_hp_tmt',
            'n': 2.0,
            'formula': {'Fe': 1.0, 'O': 1.0}}

fe_mantle_melt_params={'a_0': 2.9614332469401705e-05,
                       'K_0': 166652774642.11273,
                       'Pref': 100.e9,
                       'Cp_Pref': np.array([ 7.95326013e+01, -2.41909947e-03, -1.61692272e+06, -5.62222634e+02]),
                       'H_Pref': -195245.49100022088,
                       'Kprime_0': 5.0802472229003905,
                       'T_0': 298.15,
                       'T_einstein': 505.75,
                       'Kdprime_0': -3.9742163085937504e-11,
                       'V_0': 1.2325484447664221e-05,
                       'name': 'fe mantle',
                       'molar_mass': 0.0707624708,
                       'S_Pref': 95.0299295525918,
                       'equation_of_state': 'mod_hp_tmt',
                       'n': 2.092,
                       'formula': {'O': 1.092, 'Fe': 0.908, 'Si': 0.092}}

mg_mantle_melt_params={'a_0': 2.0572748142847914e-05,
                       'K_0': 231645314972.72287,
                       'Pref': 100.e9,
                       'Cp_Pref': np.array([ 7.95326013e+01, -2.41909947e-03, -1.61692272e+06, -5.62222634e+02]),
                       'H_Pref': -538009.8593335259,
                       'Kprime_0': 4.252832366943359,
                       'T_0': 298.15,
                       'T_einstein': 558.0924045503805,
                       'Kdprime_0': -2.1381292724609374e-11,
                       'V_0': 1.2180438865657191e-05,
                       'name': 'mg mantle',
                       'molar_mass': 0.048592178,
                       'S_Pref': 64.88469713598576,
                       'equation_of_state': 'mod_hp_tmt',
                       'n': 2.419,
                       'formula': {'O': 1.419, 'Mg': 0.581, 'Si': 0.419}}
##################### END EOS PARAMETERS ########################

################# BEGIN MELT MODEL PARAMETERS ###################
# All the lists are Mg endmember first, then Fe endmember
c_mantle = [{'MgO': 0.581, 'SiO2': 0.419},
            {'FeO': 0.908, 'SiO2': 0.092}] # molar_proportions

melting_reference_pressure = 120.e9 # Pa
melting_temperatures=np.array([4821.2, 3470.0]) # K 
melting_entropies=np.array([34.33, 33.77]) # J/K/mol-cations
melting_volumes=np.array([9.29e-08, 1.51e-07]) # m^3/mol-cations
n_mole_mix = np.array([0.52, 0.56]) # having different "n"s would be incorrect for a solid solution, but this does a slightly better job than assuming the same number of moles mixing for each "endmember"... this model is not strictly thermodynamically correct anyway.


################# END MELT MODEL PARAMETERS ###################

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

