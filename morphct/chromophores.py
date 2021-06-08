import copy
import itertools
from collections import defaultdict
import os
import sys
from warnings import warn

import ele
import freud
from openbabel import openbabel
from openbabel import pybel
import numpy as np
from scipy.spatial import Delaunay

from morphct import execute_qcc as eqcc
from morphct import helper_functions as hf


class Chromophore:
    """An object for tracking infromation related to the chromophore.

    Chromphore refers to the electronically active part of the molecule where a
    charge is expected to reside.
    This class keeps track of what particles from the simulation compose each
    chromophore, the chromphore's neighbors, and the energy values of each
    chromophore and chromophore pair.

    Parameters
    ----------
    chromo_id : int
        Index of this chromophore in the system.
    snap : gsd.hoomd.Snapshot
        Atomistic simulation snapshot from a GSD file. It is expected that the
        lengths in this file have been converted to Angstroms.
    atom_ids : numpy.ndarray of int
        Snapshot indices of the particles which belong to this chromophore.
    species : str
        Chromophore species ('donor' or 'acceptor')
    conversion_dict : dictionary
        A dictionary that maps the atom type to its element. e.g., `{'c3': C}`
        An instance that maps AMBER types to their element can be found in
        `amber_dict`.
    reorganization_energy : float, default 0.3064
        Energy (in eV) required to "reorganize" the system structure from the
        final to initial coordinates. The default value is for P3HT and came
        from DOI:10.1016/j.synthmet.2004.03.005
    vrh_delocalization : float, default 2e-10
        Variable-range hopping modifier in meters. Hopping rates are scaled by
        exp(r/vrh_delocalization) when `use_vrh` is True.

    Attributes
    ----------
    id : int
        Index of this chromophore in the system.
    species : str
        Chromophore species ('donor' or 'acceptor')
    reorganization_energy : float
        Energy required to "reorganize" the system structure from the final to
        initial coordinates.
    vrh_delocalization : float
        Variable-range hopping modifier in meters. Hopping rates are scaled by
        exp(r/vrh_delocalization) when `use_vrh` is True.
    atom_ids : numpy.ndarray of int
        Snapshot indices of the particles which belong to this chromophore.
    n_atoms : int
        The number of atoms in the chromophore.
    qcc_input : str
        The input for the MINDO3 quantum chemical calculation run in pySCF. See
        https://pyscf.org/quickstart.html for more information.
    neighbors : list of (int, numpy.ndarray(size=3))
        Each list entry is the chromophore index of the neighbor followed by the
        relative image of that neighbor. On initialization this is an empty
        list, values are set after running `set_neighbors_voronoi`.
    homo : float
        HOMO energy value in eV for this chromophore. On initialization this
        is set to None, values are set after running `set_energyvalues`.
    homo_1 : float
        HOMO-1 energy value in eV for this chromophore. On initialization this
        is set to None, values are set after running `set_energyvalues`.
    lumo : float
        LUMO energy value in eV for this chromophore. On initialization this
        is set to None, values are set after running `set_energyvalues`.
    lumo_1 : float
        LUMO+1 energy value in eV for this chromophore. On initialization this
        is set to None, values are set after running `set_energyvalues`.
    neighbors_delta_e : list of floats
        List of the differences in the frontier orbital energy values in eV
        between this chromophore and its neighbors. If the chromphore is an
        acceptor, the LUMO is used as the frontier orbital; if a donor, the HOMO
        is used. The indices are the same as in `neighbors`. On initialization
        this is an empty list, values are set after running `set_energyvalues`.
    neighbors_ti : list of floats
        List of the transfer integral energy in eV between this chromophore and
        its neighbors. See `calculate_ti` for more information. The indices are
        the same as in `neighbors`. On initialization this is an empty list,
        values are set after running `set_energyvalues`.

    Methods
    -------
    get_MO_energy
        Get the frontier orbital energy for this chromophore.
    """

    def __init__(
        self,
        chromo_id,
        snap,
        atom_ids,
        species,
        conversion_dict,
        reorganization_energy=0.3064,
        vrh_delocalization=2e-10,
    ):
        self.id = chromo_id
        if species.lower() not in ["donor", "acceptor"]:
            raise TypeError("Species must be either donor or acceptor")
        self.species = species.lower()
        self.reorganization_energy = reorganization_energy
        self.vrh_delocalization = vrh_delocalization
        self.atom_ids = atom_ids
        self.n_atoms = len(atom_ids)

        # Sets unwrapped_center, center, and image attributes
        self._set_center(snap, atom_ids)

        self.qcc_input = eqcc.write_qcc_inp(snap, atom_ids, conversion_dict)

        # Now to create a load of placeholder parameters to update later when we
        # have the full list/energy levels.
        # The self.neighbors list contains one element for each chromophore
        # within `d_cut` passed to `set_neighbors_voronoi` of this one
        # (considering periodic boundary conditions). Its format is
        # [[neighbor1_ID, relative_image_of_neighbor1],...]
        self.neighbors = []

        # The molecular orbitals of this chromophore have not yet been
        # calculated, but they will simply be floats.
        self.homo = None
        self.homo_1 = None
        self.lumo = None
        self.lumo_1 = None

        # The neighbor_delta_e and neighbor_ti are lists where each element
        # describes the difference in important molecular orbitals or transfer
        # integral between this chromophore and each neighbor. The list indices
        # here are the same as in self.neighbors for coherence.
        self.neighbors_delta_e = []
        self.neighbors_ti = []

    def __repr__(self):
        """Return the Chromophore representation."""
        return "Chromophore {} ({}): {} atoms at {:.3f} {:.3f} {:.3f}".format(
            self.id, self.species, self.n_atoms, *self.center
        )

    def _set_center(self, snap, atom_ids):
        box = snap.configuration.box[:3]
        unwrapped_pos = snap.particles.position + snap.particles.image * box
        center = np.mean(unwrapped_pos[atom_ids], axis=0)
        img = np.zeros(3)
        while (center + img * box < -box / 2).any() or (
            center + img * box > box / 2
        ).any():
            img[np.where(center < -box / 2)] += 1
            img[np.where(center > box / 2)] -= 1
        self.unwrapped_center = center
        self.center = center + img * box
        self.image = img

    def get_MO_energy(self):
        """Get the frontier molecular orbital energy for this chromophore.

        If the chromphore is an acceptor, the LUMO energy is returned. If the
        chromphore is a donor, the HOMO energy is returned.

        Returns
        -------
        float
            MO energy in eV
        """
        if self.species == "acceptor":
            return self.lumo
        return self.homo


def get_chromo_ids_smiles(snap, smarts_str, conv_dict):
    """Get the atom indices in a snapshot associated with a SMARTS string.

    This function can be used to determine the atom indices for each
    chromophore. SMARTS matching depends on the molecular structures making
    chemical sense (e.g., aromatic structures are planar, etc). Often snapshots
    from molecular simulations based on classical methods (e.g., MC, MD) may
    have distortions that are chemically unphysical, in which case this function
    may not find all chromophores. A solution is to use this function on a
    snapshot of the initial frame of the trajectory, and then apply these
    indices to a later frame.

    Parameters
    ----------
    snap : gsd.hoomd.Snapshot
        Atomistic simulation snapshot from a GSD file. It is expected that the
        lengths in this file have been converted to Angstroms.
    smarts_str : str
        SMARTS string used to find the atom indices.
    conversion_dict : dictionary
        A dictionary that maps the atom type to its element. e.g., `{'c3': C}`
        An instance that maps AMBER types to their element can be found in
        `amber_dict`.

    Returns
    -------
    list of numpy.ndarray of int
        atom indices of each SMARTS match

    Note
    ----
    If no matches are found, a warning is raised and the pybel.Molecule object
    is returned for debugging.
    """
    box = snap.configuration.box[:3]
    unwrapped_positions = snap.particles.position + snap.particles.image * box

    mol = openbabel.OBMol()
    for i, typeid in enumerate(snap.particles.typeid):
        a = mol.NewAtom()
        element = conv_dict[snap.particles.types[typeid]]
        a.SetAtomicNum(element.atomic_number)
        a.SetVector(*[float(x) for x in unwrapped_positions[i]])

    for i, j in snap.bonds.group:
        # openbabel indexes atoms from 1
        # AddBond(i_index, j_index, bond_order)
        mol.AddBond(int(i + 1), int(j + 1), 1)

    # This will correctly set the bond order
    # (necessary for smarts matching)
    mol.PerceiveBondOrders()
    mol.SetAromaticPerceived()

    pybelmol = pybel.Molecule(mol)

    smarts = pybel.Smarts(smarts_str)
    # shift indices by 1
    atom_ids = [np.array(i) - 1 for i in smarts.findall(pybelmol)]
    if not atom_ids:
        warn(
            f"No matches found for smarts string {smarts_str}. "
            + "Please check the returned pybel.Molecule for errors.\n"
        )
        return pybelmol
    print(f"Found {len(atom_ids)} chromophores.")
    return atom_ids


def set_neighbors_voronoi(chromo_list, snap, conversion_dict, d_cut=10):
    """Set the chromophore neighbors using voronoi analysis.

    See https://freud.readthedocs.io/en/latest/modules/locality.html#freud.locality.Voronoi
    for more information.

    Parameters
    ----------
    chromo_list : list of Chromphore
        List of the Chromophore objects between which to calculate neighbors.
    snap : gsd.hoomd.Snapshot
        Atomistic simulation snapshot from a GSD file. It is expected that the
        lengths in this file have been converted to Angstroms.
    conversion_dict : dictionary
        A dictionary that maps the atom type to its element. e.g., `{'c3': C}`
        An instance that maps AMBER types to their element can be found in
        `amber_dict`.
    d_cut : float, default 10
        The distance cutoff for neighbors.

    Returns
    -------
    qcc_pairs : list of ((int,int), str)
        The information needed for calculating the pair energies. The first part
        of each entry is the pair indices followed by the input for the MINDO3
        quantum chemical calculation run in pySCF for the chromophore pairs.
        See https://pyscf.org/quickstart.html for more information.
    """
    voronoi = freud.locality.Voronoi()
    freudbox = freud.box.Box(*snap.configuration.box)
    centers = [chromo.center for chromo in chromo_list]
    voronoi.compute((freudbox, centers))

    box = snap.configuration.box[:3]
    qcc_pairs = []
    neighbors = []
    for (i, j) in voronoi.nlist:
        if i == j:
            pass
        elif (i, j) not in neighbors and (j, i) not in neighbors:
            chromo_i = chromo_list[i]
            chromo_j = chromo_list[j]
            centers = []
            distances = []
            images = []
            # calculate which of the periodic image is closest
            # shift chromophore j, hold chromophore i in place
            for xyz_image in itertools.product(range(-1, 2), repeat=3):
                xyz_image = np.array(xyz_image)
                sc_center = chromo_j.center + xyz_image * box
                images.append(xyz_image)
                centers.append(sc_center)
                distances.append(np.linalg.norm(sc_center - chromo_i.center))
            imin = distances.index(min(distances))
            if distances[imin] > d_cut:
                continue

            rel_image = images[imin]
            j_shift = centers[imin] - chromo_j.unwrapped_center
            chromo_i.neighbors.append([j, rel_image])
            chromo_i.neighbors_delta_e.append(None)
            chromo_i.neighbors_ti.append(None)
            chromo_j.neighbors.append([i, -rel_image])
            chromo_j.neighbors_delta_e.append(None)
            chromo_j.neighbors_ti.append(None)
            neighbors.append((i, j))
            qcc_input = eqcc.write_qcc_pair_input(
                snap, chromo_i, chromo_j, j_shift, conversion_dict
            )
            qcc_pairs.append(((i, j), qcc_input))
    return qcc_pairs


conversion_dict = {
    "S1": ele.element_from_symbol("S"),
    "H1": ele.element_from_symbol("H"),
    "C5": ele.element_from_symbol("C"),
    "C1": ele.element_from_symbol("C"),
    "C4": ele.element_from_symbol("C"),
    "C6": ele.element_from_symbol("C"),
    "C8": ele.element_from_symbol("C"),
    "C9": ele.element_from_symbol("C"),
    "C3": ele.element_from_symbol("C"),
    "C7": ele.element_from_symbol("C"),
    "C2": ele.element_from_symbol("C"),
    "C10": ele.element_from_symbol("C"),
}


amber_dict = {
    "c": ele.element_from_symbol("C"),
    "c1": ele.element_from_symbol("C"),
    "c2": ele.element_from_symbol("C"),
    "c3": ele.element_from_symbol("C"),
    "ca": ele.element_from_symbol("C"),
    "cp": ele.element_from_symbol("C"),
    "cq": ele.element_from_symbol("C"),
    "cc": ele.element_from_symbol("C"),
    "cd": ele.element_from_symbol("C"),
    "ce": ele.element_from_symbol("C"),
    "cf": ele.element_from_symbol("C"),
    "cg": ele.element_from_symbol("C"),
    "ch": ele.element_from_symbol("C"),
    "cx": ele.element_from_symbol("C"),
    "cy": ele.element_from_symbol("C"),
    "cu": ele.element_from_symbol("C"),
    "cv": ele.element_from_symbol("C"),
    "h1": ele.element_from_symbol("H"),
    "h2": ele.element_from_symbol("H"),
    "h3": ele.element_from_symbol("H"),
    "h4": ele.element_from_symbol("H"),
    "h5": ele.element_from_symbol("H"),
    "ha": ele.element_from_symbol("H"),
    "hc": ele.element_from_symbol("H"),
    "hn": ele.element_from_symbol("H"),
    "ho": ele.element_from_symbol("H"),
    "hp": ele.element_from_symbol("H"),
    "hs": ele.element_from_symbol("H"),
    "hw": ele.element_from_symbol("H"),
    "hx": ele.element_from_symbol("H"),
    "f": ele.element_from_symbol("F"),
    "cl": ele.element_from_symbol("Cl"),
    "br": ele.element_from_symbol("Br"),
    "i": ele.element_from_symbol("I"),
    "n": ele.element_from_symbol("N"),
    "n1": ele.element_from_symbol("N"),
    "n2": ele.element_from_symbol("N"),
    "n3": ele.element_from_symbol("N"),
    "n4": ele.element_from_symbol("N"),
    "na": ele.element_from_symbol("N"),
    "nb": ele.element_from_symbol("N"),
    "nc": ele.element_from_symbol("N"),
    "nd": ele.element_from_symbol("N"),
    "ne": ele.element_from_symbol("N"),
    "nf": ele.element_from_symbol("N"),
    "nh": ele.element_from_symbol("N"),
    "no": ele.element_from_symbol("N"),
    "o": ele.element_from_symbol("O"),
    "oh": ele.element_from_symbol("O"),
    "os": ele.element_from_symbol("O"),
    "ow": ele.element_from_symbol("O"),
    "p2": ele.element_from_symbol("P"),
    "p3": ele.element_from_symbol("P"),
    "p4": ele.element_from_symbol("P"),
    "p5": ele.element_from_symbol("P"),
    "pb": ele.element_from_symbol("P"),
    "pc": ele.element_from_symbol("P"),
    "pd": ele.element_from_symbol("P"),
    "pe": ele.element_from_symbol("P"),
    "pf": ele.element_from_symbol("P"),
    "px": ele.element_from_symbol("P"),
    "py": ele.element_from_symbol("P"),
    "s": ele.element_from_symbol("S"),
    "s2": ele.element_from_symbol("S"),
    "s4": ele.element_from_symbol("S"),
    "s6": ele.element_from_symbol("S"),
    "sh": ele.element_from_symbol("S"),
    "ss": ele.element_from_symbol("S"),
    "sx": ele.element_from_symbol("S"),
    "sy": ele.element_from_symbol("S"),
}
