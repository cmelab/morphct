import glob
import os
import pickle
import shutil
import sys
import numpy as np


def calculate_delta_E(chromo1, chromo2):
    return chromo2.get_MO_energy() - chromo1.get_MO_energy()


def calculate_ti(orbital_splitting, delta_e):
    # Use the energy splitting in dimer method to calculate the electronic
    # transfer integral in eV
    if delta_e ** 2 > orbital_splitting ** 2:
        # Avoid an imaginary TI by returning zero.
        # (Could use KOOPMAN'S APPROXIMATION here if desired)
        return 0
    return 0.5 * np.sqrt((orbital_splitting ** 2) - (delta_e ** 2))
