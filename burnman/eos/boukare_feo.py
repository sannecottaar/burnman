from __future__ import absolute_import
# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2017 by the BurnMan team, released under the GNU
# GPL v2 or later.


import numpy as np
import warnings

from . import modified_tait as mt
from . import equation_of_state as eos

from . import einstein


class BOUKARE_FEO(eos.EquationOfState):

    """
    Class for the Boukare FeO liquid equation of state
    """

    def volume(self, pressure, temperature, params):
        """
        Returns volume [m^3] as a function of pressure [Pa] and temperature [K]
        """
        return ( (params['molar_mass'] / params['rho_0']) *
                 (np.power(1. + params['en']*pressure/params['K_0'],
                           -1./params['en']) +
                  np.power(1. + params['en']*pressure/params['K_0'],
                           -(params['en'] + params['q_0'])/params['en'])*params['a_0'] *
                  (temperature - params['T_0'])) )

    def pressure(self, temperature, volume, params):
        """
        Returns pressure [Pa] as a function of temperature [K] and volume[m^3]
        """
        _delta_volume = lambda p, args: args[0] - self.volume(p, args[1], args[2])
        return opt.brentq(_delta_volume, 0., 500.e9, args=(volume, temperature, params))

    def isothermal_bulk_modulus(self, pressure, temperature, volume, params):
        """
        Returns isothermal bulk modulus [Pa] as a function of pressure [Pa],
        temperature [K], and volume [m^3].  EQ 13+2
        """
        # Value found by numerical differentiation (fine for now)
        # TODO: Use analytical solution
        dP = 1000.
        V0 = self.volume(pressure-dP/2., temperature, params)
        V1 = self.volume(pressure+dP/2., temperature, params)
        return -(((V1 + V0) / 2.) * (dP/(V1 - V0)))

    # calculate the shear modulus as a function of P, V, and T
    def shear_modulus(self, pressure, temperature, volume, params):
        """
        Shear modulus. Exactly zero for a liquid
        """
        return 0.

    def thermal_expansivity(self, pressure, temperature, volume, params):
        """
        Returns thermal expansivity at the pressure, temperature, and volume [1/K]
        """
        # Value found by numerical differentiation (fine for now)
        # TODO: Use analytical solution
        dT = 0.1
        V0 = self.volume(pressure, temperature-dT/2., params)
        V1 = self.volume(pressure, temperature+dT/2., params)
        return ((2. / (V1 + V0)) * ((V1 - V0)/dT))


    def gibbs_free_energy(self, pressure, temperature, volume, params):
        """
        Returns the gibbs free energy [J/mol] as a function of pressure [Pa]
        and temperature [K].
        """
        return (params['G_0'] - 
                params['S_0']*(temperature - params['T_0']) -
                params['C_1']*(temperature*np.log(temperature/params['T_0']) -
                               (temperature - params['T_0'])) -
                params['C_0']/2.*(temperature - params['T_0'])*(temperature - params['T_0']) -
                params['molar_mass']*params['a_0']*params['K_0']/(params['rho_0']*params['q_0']) *
                (np.power((1. + params['en']*(pressure/params['K_0'])),
                          -params['q_0']/params['en']) -
                 np.power((1. + params['en']*params['P_0']/params['K_0']),
                          -params['q_0']/params['en'])) *
                (temperature - params['T_0']) +
                params['molar_mass']*params['K_0']/(params['rho_0']*(params['en'] - 1.)) *
                (np.power((1. + params['en']*pressure/params['K_0']),
                          1. - 1./params['en']) -
                 np.power((1. + params['en']*params['P_0']/params['K_0']),
                          1. - 1./params['en'])))
        

    def entropy(self, pressure, temperature, volume, params):
        """
        Returns the entropy [J/K/mol] as a function of pressure [Pa]
        and temperature [K].
        """
        return (params['S_0'] +
                params['C_0']*(temperature - params['T_0']) + 
                params['C_1']*np.log(temperature/params['T_0']) +
                params['molar_mass']*params['a_0']*params['K_0']/(params['rho_0']*params['q_0']) *
                (np.power((1. + params['en']*pressure/params['K_0']), -params['q_0']/params['en']) -
                 np.power((1. + params['en']*params['P_0']/params['K_0']), -params['q_0']/params['en'])))
    

    def enthalpy(self, pressure, temperature, volume, params):
        """
        Returns the enthalpy [J/mol] as a function of pressure [Pa]
        and temperature [K].
        """
        gibbs = self.gibbs_free_energy(pressure, temperature, volume, params)
        entropy = self.entropy(pressure, temperature, volume, params)
        return gibbs + temperature * entropy

    def molar_heat_capacity_p(self, pressure, temperature, volume, params):
        """
        Returns the heat capacity [J/K/mol] as a function of pressure [Pa]
        and temperature [K].
        """
        # Value found by numerical differentiation (fine for now)
        # TODO: Use analytical solution
        dT = 0.1
        S0 = self.entropy(pressure, temperature-dT/2., volume, params)
        S1 = self.entropy(pressure, temperature+dT/2., volume, params)
        return temperature*(S1 - S0)/dT

    def validate_parameters(self, params):
        """
        Check for existence and validity of the parameters
        """

        # Check all the required keys are in the dictionary
        expected_keys = ['T_0', 'P_0', 'molar_mass',
                         'rho_0', 'K_0', 'a_0', 'q_0', 'en',
                         'S_0', 'C_0', 'C_1', 'G_0']
        
        for k in expected_keys:
            if k not in params:
                raise KeyError('params object missing parameter : ' + k)

        # Check that the values are reasonable.
        if params['T_0'] < 0.:
            warnings.warn('Unusual value for T_0', stacklevel=2)
