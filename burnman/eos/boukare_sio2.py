# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2017 by the BurnMan team, released under the GNU
# GPL v2 or later.

from __future__ import absolute_import

import numpy as np
import scipy.optimize as opt
import warnings

# Try to import the jit from numba.  If it is
# not available, just go with the standard
# python interpreter
try:
    from numba import jit
except ImportError:
    def jit(fn):
        return fn


from . import birch_murnaghan as bm
from . import debye
from . import equation_of_state as eos
from ..tools import bracket


@jit
def _delta_pressure(x, pressure, temperature, T_0, p, f, Cv1):
    # x = V
    P_static = p[0]*np.exp(p[1]*x) + p[2]*np.exp(p[3]*x) # N.B. THIS DOES NOT EVALUATE TO ZERO AT V_O!!!
    P_th = ((temperature - T_0)*f[0]*pow(x, f[1]) +
            Cv1*(temperature*np.log(temperature / T_0) -
                 (temperature - T_0))) 

    return P_static + P_th - pressure  # EQ 21


class BOUKARE_SIO2(eos.EquationOfState):

    """
    Class for the Boukare SiO2 liquid
    """

    # calculate isotropic thermal pressure
    def _thermal_pressure(self, T, V, params):
        P_th = ((T - params['T_0'])*params['f'][0]*pow(V, params['f'][1]) +
                params['Cv'][1]*(T*np.log(T / params['T_0']) -
                                 (T - params['T_0'])))
        return P_th

    def volume(self, pressure, temperature, params):
        """
        Returns molar volume. :math:`[m^3]`
        """
        # we need to have a sign change in [a,b] to find a zero. Let us start with a
        # conservative guess:
        args = (pressure, temperature, params['T_0'],
                params['p'], params['f'], params['Cv'][1])
        try:
            sol = bracket(_delta_pressure, params[
                          'V_0'], 1.e-2 * params['V_0'], args)
        except ValueError:
            raise Exception(
                'Cannot find a volume, perhaps you are outside of the range of validity for the equation of state?')
        return opt.brentq(_delta_pressure, sol[0], sol[1], args=args)

    def pressure(self, temperature, volume, params):
        """
        Returns the pressure of the mineral at a given temperature and volume [Pa]
        """
        p = params['p']
        P_static = p[0]*np.exp(p[1]*volume) + p[2]*np.exp(p[3]*volume) # N.B. THIS EVALUATES TO ABOUT 1.591 GPa AT V_0
        P_th = self._thermal_pressure(temperature, volume, params)
        return P_static + P_th


    def isothermal_bulk_modulus(self, pressure, temperature, volume, params):
        """
        Returns isothermal bulk modulus :math:`[Pa]`
        """
        # Value found by numerical differentiation (fine for now)
        # TODO: Use analytical solution
        dP = 1000.
        V0 = self.volume(pressure-dP/2., temperature, params)
        V1 = self.volume(pressure+dP/2., temperature, params)
        return -(((V1 + V0) / 2.) * (dP/(V1 - V0)))


    def shear_modulus(self, pressure, temperature, volume, params):
        """
        Shear modulus. Exactly zero for a liquid
        """
        return 0.

    def molar_heat_capacity_v(self, pressure, temperature, volume, params):
        """
        Returns heat capacity at constant volume. :math:`[J/K/mol]`
        """
        return params['Cv'][0] + params['Cv'][1]*volume + params['Cv'][2]*temperature

    def grueneisen_parameter(self, pressure, temperature, volume, params):
        """
        Returns heat capacity at constant pressure. :math:`[J/K/mol]`
        """
        alpha = self.thermal_expansivity(pressure, temperature, volume, params)
        K = self.isothermal_bulk_modulus(pressure, temperature, volume, params)
        C_v = self.molar_heat_capacity_v(pressure, temperature, volume, params)
        
        return alpha * K * volume / C_v
    
    def molar_heat_capacity_p(self, pressure, temperature, volume, params):
        """
        Returns heat capacity at constant pressure. :math:`[J/K/mol]`
        """
        alpha = self.thermal_expansivity(pressure, temperature, volume, params)
        K = self.isothermal_bulk_modulus(pressure, temperature, volume, params)
        C_v = self.molar_heat_capacity_v(pressure, temperature, volume, params)
        C_p = C_v + alpha * alpha * volume * K * temperature
        return C_p

    def thermal_expansivity(self, pressure, temperature, volume, params):
        """
        Returns thermal expansivity. :math:`[1/K]`
        """
        K = self.isothermal_bulk_modulus(pressure, temperature, volume, params)
        alpha = (params['f'][0]*pow(volume, params['f'][1]) +
                 params['Cv'][1]*np.log(temperature/params['T_0'])) / K


        # Numerical differentiation
        
        dT = 0.1
        V0 = self.volume(pressure, temperature-dT/2., params)
        V1 = self.volume(pressure, temperature+dT/2., params)
        alpha = ((2. / (V1 + V0)) * ((V1 - V0)/dT))
        return alpha

    def gibbs_free_energy(self, pressure, temperature, volume, params):
        """
        Returns the Gibbs free energy at the pressure and temperature of the mineral [J/mol]
        """
        G = self.helmholtz_free_energy(
            pressure, temperature, volume, params) + pressure * volume
        return G

    def molar_internal_energy(self, pressure, temperature, volume, params):
        """
        Returns the internal energy at the pressure and temperature of the mineral [J/mol]
        """
        return self.helmholtz_free_energy(pressure, temperature, volume, params) + \
            temperature * \
            self.entropy(pressure, temperature, volume, params)

    def entropy(self, pressure, temperature, volume, params):
        """
        Returns the entropy at the pressure and temperature of the mineral [J/K/mol]
        """
        # Note that equation 21 in Boukare notes is incorrect
        # dS = a*K_T dV + Cv/T dT
        # Cv = Cv0 + Cv1*V + Cv2*T
        # aKT = f1*V^f2 + CV1 ln(T/T0)
        
        # S = S_0 + f1/(f2+1)*(V^(f2+1) - V0^(f2+1)) + (Cv0 + Cv1*V)*ln(T/T0) + Cv2*(T - T0)
        f = params['f']
        S = (params['S_0'] +
             f[0]/(f[1] + 1.)*(pow(volume, f[1]+1.) - pow(params['V_0'], f[1]+1.)) + 
             (params['Cv'][0] + params['Cv'][1]*volume)*np.log(temperature/
                                                               params['T_0']) +
             params['Cv'][2]*(temperature - params['T_0']))

        return S

    def enthalpy(self, pressure, temperature, volume, params):
        """
        Returns the enthalpy at the pressure and temperature of the mineral [J/mol]
        """

        return self.helmholtz_free_energy(pressure, temperature, volume, params) + \
            temperature * self.entropy(pressure, temperature, volume, params) + \
            pressure * volume

    def _thermal_helmholtz(self, T, V, params):
        f = params['f']
        a = (params['S_0'] + f[0]/(f[1] + 1.) *
             (pow(V, f[1]+1.) - pow(params['V_0'], f[1]+1.)))
        b = (params['Cv'][0] + params['Cv'][1]*V)
        c = params['Cv'][2]

        F_th = -((T - params['T_0']) *
                 (a - b + 0.5*c*(T - params['T_0'])) +
                 b*T*np.log(T/params['T_0']))

        return F_th
    
    def helmholtz_free_energy(self, pressure, temperature, volume, params):
        """
        Returns the Helmholtz free energy at the pressure and temperature of the mineral [J/mol]
        """
        p = params['p']
        F_static = -(p[0]/p[1]*(np.exp(p[1]*volume) - np.exp(p[1]*params['V_0'])) +
                     p[2]/p[3]*(np.exp(p[3]*volume) - np.exp(p[3]*params['V_0'])))
        return params['F_0'] + F_static + self._thermal_helmholtz(temperature, volume, params)

    def validate_parameters(self, params):
        """
        Check for existence and validity of the parameters
        """

        # Now check all the required keys for the EoS are in the dictionary
        expected_keys = ['molar_mass', 'T_0', 'V_0', 'f', 'p',
                         'Cv', 'S_0', 'F_0']
        for k in expected_keys:
            if k not in params:
                raise KeyError('params object missing parameter : ' + k)


