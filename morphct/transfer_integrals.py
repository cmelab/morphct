import glob
import os
import pickle
import shutil
import sys
import numpy as np


def calculate_delta_E(chromo1, chromo2):
    """Calculate the energy difference between the chromophores.

    Parameters
    ----------
    chromo1 : Chromophore
        A chromophore.
    chromo2 : Chromophore
        Another chromophore.

    Returns
    -------
    float
        The energy difference between the frontier orbitals in eV
    """
    return chromo2.get_MO_energy() - chromo1.get_MO_energy()


def calculate_ti(orbital_splitting, delta_e):
    """Calculate the electronic transfer integral in eV.

    Use the energy splitting in dimer method to calculate the transfer integral.
    If the transfer integral would be imaginary, return zero.

    Parameters
    ----------
    orbital_splitting : float
        The energy level splitting induced by including the neighbor chromophore
        in the qcc calculation. Units are eV.
    delta_e : float
        The energy difference between the frontier orbitals of the chromophores
        in eV.

    Returns
    -------
    float
        The electronic transfer integral.
    """
    # Use the energy splitting in dimer method to calculate the electronic
    # transfer integral in eV
    if delta_e ** 2 > orbital_splitting ** 2:
        # Avoid an imaginary TI by returning zero.
        # (Could use KOOPMAN'S APPROXIMATION here if desired)
        return 0
    return 0.5 * np.sqrt((orbital_splitting ** 2) - (delta_e ** 2))
