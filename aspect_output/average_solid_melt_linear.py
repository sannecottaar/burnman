# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2017 by the BurnMan team, released under the GNU
# GPL v2 or later.

from __future__ import absolute_import
import numpy as np
import warnings
import sys,os
sys.path.insert(1, os.path.abspath('../burnman'))
from averaging_schemes import AveragingScheme


class equilibrium_geometry(AveragingScheme):

    """
    Class for computing the elastic properties for a partial molten case, assuming an equilibrium spheroid model, currently simplified to a linear relationship for a dihedral angle of 80.
    This derives from :class:`burnman.averaging_schemes.averaging_scheme`, and implements
    the :func:`burnman.averaging_schemes.averaging_scheme.average_bulk_moduli` and
    :func:`burnman.averaging_schemes.averaging_scheme.average_shear_moduli` functions.
    """

    def average_bulk_moduli(self, volumes, bulk_moduli, shear_moduli):
        """
        Average the bulk moduli of a composite with the Reuss (iso-stress)
        bound, given by:

        .. math::
            K_R = \\left(\\Sigma_i \\frac{V_i}{K_i} \\right)^{-1}

        Parameters
        ----------
        volumes : list of floats
            [volume solid, volume fluid] :math:`[m^3]`
        bulk_moduli : list of floats
            [K solid, K fluid]:math:`[Pa]`
        shear_moduli : list of floats
            [G solid, G fluid] :math:`G` of each phase in the composite.
            Not used in this average.
        Returns
        -------

        Kb : float
            The bulk modulus for partial molten case :math:`K_b`. :math:`[Pa]`
        """
        volume_fraction = volumes/sum(volumes)
        Kb = (1.-0.5*volume_fraction[1]/0.072)*bulk_moduli[0]
        print(volume_fraction, bulk_moduli, Kb)
        
        
        if Kb<0:
            Kb =0
        
        return Kb

    def average_shear_moduli(self, volumes, bulk_moduli, shear_moduli):
        """
        Average the shear moduli of a composite with the Reuss (iso-stress)
        bound, given by:

        .. math::
            G_R = \\left( \\Sigma_i \\frac{V_i}{G_i} \\right)^{-1}

        Parameters
        ----------
        volumes : list of floats
            [volume solid, volume fluid] :math:`[m^3]`
        bulk_moduli : list of floats
            [K solid, K fluid]:math:`[Pa]`
        shear_moduli : list of floats
            [G solid, G fluid] :math:`G` of each phase in the composite.

        Returns
        -------

        N : float
            The shear modulus for partial molten case :math:`G_R`. :math:`[Pa]`
        """
        volume_fraction = volumes/sum(volumes)
        N = (1.-0.5*volume_fraction[1]/0.102)*shear_moduli[0]
        if N<0:
            N=0
        return N

