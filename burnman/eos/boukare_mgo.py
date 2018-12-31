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
def _grueneisen_parameter_fast(V_0, volume, gruen_0, q_0):
    """global function with plain parameters so jit will work"""
    return gruen_0*np.power(volume/V_0, q_0)


@jit
def _delta_pressure(x, pressure, temperature, V_0, T_0, b_iikk, b_iikkmm, Cv_0, gruen_0, q_0):
    # x = V
    # Remember that dF/dV = -P
    f = 0.5 * (pow(V_0 / x, 2. / 3.) - 1.)
    P_static = (1. / 3.) * (pow(1. + 2. * f, 5. / 2.)) * ((b_iikk * f) + (0.5 * b_iikkmm * f * f))
    P_th = Cv_0*gruen_0/x*(temperature - T_0)*pow(x/V_0, q_0)

    return P_static + P_th - pressure  # EQ 21


class BOUKARE_MGO(eos.EquationOfState):

    """
    Class for the Boukare MgO liquid
    """

    # calculate isotropic thermal pressure
    def _thermal_pressure(self, T, V, params):
        P_th = (params['Cv_0'] * params['grueneisen_0'] / V *
                (T - params['T_0'])*pow(V/params['V_0'], params['q_0']))
        return P_th

    def volume(self, pressure, temperature, params):
        """
        Returns molar volume. :math:`[m^3]`
        """
        T_0 = params['T_0']
        V_0 = params['V_0']

        b_iikk = 9. * params['K_0']  # EQ 28
        b_iikkmm = 27. * params['K_0'] * (params['Kprime_0'] - 4.)  # EQ 29z

        Cv_0 = params['Cv_0']
        gruen_0 = params['grueneisen_0']
        q_0 = params['q_0']
        
        # we need to have a sign change in [a,b] to find a zero. Let us start with a
        # conservative guess:
        args = (pressure, temperature, V_0, T_0,
                b_iikk, b_iikkmm, Cv_0, gruen_0, q_0)
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
        
        P_th = self._thermal_pressure(temperature, volume, params)
        
        b_iikk = 9. * params['K_0']  # EQ 28
        b_iikkmm = 27. * params['K_0'] * (params['Kprime_0'] - 4.)  # EQ 29
        f = 0.5 * (pow(params['V_0'] / volume, 2. / 3.) - 1.)  # EQ 24
        P = (1. / 3.) * (pow(1. + 2. * f, 5. / 2.)) \
            * ((b_iikk * f) + (0.5 * b_iikkmm * pow(f, 2.)))\
            + P_th
        
        return P

    def grueneisen_parameter(self, pressure, temperature, volume, params):
        """
        Returns grueneisen parameter :math:`[unitless]`
        """
        return _grueneisen_parameter_fast(params['V_0'], volume, params['grueneisen_0'], params['q_0'])

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
        return params['Cv_0']

    def molar_heat_capacity_p(self, pressure, temperature, volume, params):
        """
        Returns heat capacity at constant pressure. :math:`[J/K/mol]`
        """
        alpha = self.thermal_expansivity(pressure, temperature, volume, params)
        gr = self.grueneisen_parameter(pressure, temperature, volume, params)
        C_v = self.molar_heat_capacity_v(pressure, temperature, volume, params)
        C_p = C_v * (1. + gr * alpha * temperature)
        return C_p

    def thermal_expansivity(self, pressure, temperature, volume, params):
        """
        Returns thermal expansivity. :math:`[1/K]`
        """
        C_v = self.molar_heat_capacity_v(pressure, temperature, volume, params)
        gr = self.grueneisen_parameter(pressure, temperature, volume, params)
        K = self.isothermal_bulk_modulus(pressure, temperature, volume, params)
        alpha = gr * C_v / K / volume
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
        S = (params['S_0'] +
             (params['Cv_0']*params['grueneisen_0']/params['q_0'] *
              (pow(volume/params['V_0'], params['q_0']) - 1.)) +
             params['Cv_0']*np.log(temperature/params['T_0']))
        return S

    def enthalpy(self, pressure, temperature, volume, params):
        """
        Returns the enthalpy at the pressure and temperature of the mineral [J/mol]
        """

        return self.helmholtz_free_energy(pressure, temperature, volume, params) + \
            temperature * self.entropy(pressure, temperature, volume, params) + \
            pressure * volume

    def helmholtz_free_energy(self, pressure, temperature, volume, params):
        """
        Returns the Helmholtz free energy at the pressure and temperature of the mineral [J/mol]
        """
        x = params['V_0'] / volume
        f = 1. / 2. * (pow(x, 2. / 3.) - 1.)
        b_iikk = 9. * params['K_0']  # EQ 28
        b_iikkmm = 27. * params['K_0'] * (params['Kprime_0'] - 4.)  # EQ 29

        F_th = (-params['S_0']*(temperature - params['T_0']) -
                params['Cv_0']*(temperature*np.log(temperature/params['T_0']) -
                                (temperature - params['T_0'])) - 
                (params['Cv_0']*params['grueneisen_0']/params['q_0'] *
                 (temperature - params['T_0']) *
                 (pow(volume/params['V_0'], params['q_0']) - 1.)))
        F = (params['F_0'] + 
            (0.5 * b_iikk * f * f * params['V_0'] +
             (1. / 6.) * params['V_0'] * b_iikkmm * f * f * f) +
             F_th)

        return F

    def validate_parameters(self, params):
        """
        Check for existence and validity of the parameters
        """

        # Now check all the required keys for the EoS are in the dictionary
        expected_keys = ['molar_mass', 'T_0', 'V_0', 'K_0', 'Kprime_0',
                         'grueneisen_0', 'q_0', 'Cv_0', 'S_0', 'F_0']
        for k in expected_keys:
            if k not in params:
                raise KeyError('params object missing parameter : ' + k)

