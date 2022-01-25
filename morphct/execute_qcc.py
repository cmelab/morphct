import multiprocessing as mp
from multiprocessing import get_context

import ele
import numpy as np
import pyscf
from pyscf.semiempirical import MINDO3

from morphct import helper_functions as hf
from morphct import transfer_integrals as ti


def get_homolumo(molstr, charge=0, verbose=0, tol=1e-6):
    """Get the HOMO-1, HOMO, LUMO, LUMO+1 energies in eV using MINDO3.

    See https://pyscf.org/quickstart.html for more information.

    Parameters
    ----------
    molstr : str
        Input string for pySCF containing elements and positions in Angstroms
        (e.g., "C 0.0 0.0 0.0; H 1.54 0.0 0.0")
    charge : int, default 0
        If the molecule which we are calculating the energies of has a charge,
        it can be specified.
    verbose : int, default 0
        Verbosity level of the MINDO calculation output. 0 will silence output,
        4 will show convergence.
    tol : float, default 1e-6
        Tolerance of the MINDO convergence.

    Returns
    -------
    numpy.ndarray
        Array containing HOMO-1, HOMO, LUMO, LUMO+1 energies in eV
    """
    mol = pyscf.M(atom=molstr, charge=charge)
    mf = MINDO3(mol).run(verbose=verbose, conv_tol=tol)
    occ = mf.get_occ()
    i_lumo = np.argmax(occ < 1)
    energies = mf.mo_energy[i_lumo - 2 : i_lumo + 2]
    energies *= 27.2114  # convert Eh to eV
    return energies


def singles_homolumo(chromo_list, filename=None, nprocs=None):
    """Get the HOMO-1, HOMO, LUMO, LUMO+1 energies for all single chromophores.

    Parameters
    ----------
    chromo_list : list of Chromophore
        Chromophores to calculate energies of. Each Chromophore must have
        qcc_input attribute set.
    filename : str, default None
        Path to file where singles energies will be saved. If None, energies
        will not be saved.
    nprocs : int, default None
        Number of processes passed to multiprocessing.Pool. If None,
        multiprocessing.cpu_count will used to determine optimal number.

    Returns
    -------
    data : numpy.ndarray
        Array of energies where each row corresponds to the MO energies of each
        chromophore in the list.
    """
    if nprocs is None:
        nprocs = mp.cpu_count()
    with get_context("spawn").Pool(processes=nprocs) as p:
        data = p.map(
            _worker_wrapper, [(i.qcc_input, i.charge) for i in chromo_list]
        )

    data = np.stack(data)
    if filename is not None:
        np.savetxt(filename, data)
    return data


def dimer_homolumo(qcc_pairs, chromo_list, filename=None, nprocs=None):
    """Get the HOMO-1, HOMO, LUMO, LUMO+1 energies for all chromophore pairs.

    Parameters
    ----------
    qcc_pairs : list of ((int, int), str)
        Each list item contains a tuple with the indices of the pair and the
        qcc input string.
        qcc_pairs is returned by `morphct.chromophores.set_neighbors_voronoi`
    chromo_list : list of Chromophore
        List of chromphores to calculate dimer energies.
    filename : str, default None
        Path to file where the pair energies will be saved. If None, energies
        will not be saved.
    nprocs : int, default None
        Number of processes passed to multiprocessing.Pool. If None,
        multiprocessing.cpu_count will used to determine optimal number.

    Returns
    -------
    dimer_data : list of ((int, int), numpy.ndarray)
        Each list item contains the indices of the pair and an array of its MO
        energies.
    """
    if nprocs is None:
        nprocs = mp.cpu_count()

    with get_context("spawn").Pool(processes=nprocs) as p:
        args = [
            (qcc_input, chromo_list[i].charge + chromo_list[j].charge)
            for (i,j), qcc_input in qcc_pairs
        ]
        data = p.map(_worker_wrapper, args)

    dimer_data = [i for i in zip([pair for pair, qcc_input in qcc_pairs], data)]
    if filename is not None:
        with open(filename, "w") as f:
            f.writelines(
                f"{pair[0]} {pair[1]} {en[0]} {en[1]} {en[2]} {en[3]}\n"
                for pair, en in dimer_data
            )
    return dimer_data


def get_dimerdata(filename):
    """Read in the saved data created by `dimer_homolumo`.

    Parameters
    ----------
    filename : str
        Path to file where the dimer energies were saved.

    Returns
    -------
    dimer_data : list of ((int, int), numpy.ndarray)
        Each list item contains the indices of the pair and an array of its MO
        energies.
    """
    dimer_data = []
    with open(filename, "r") as f:
        for i in f.readlines():
            a, b, c, d, e, f = i.split()
            dimer_data.append(
                ((int(a), int(b)), (float(c), float(d), float(e), float(f)))
            )
    return dimer_data


def get_singlesdata(filename):
    """Read in the saved data created by `singles_homolumo`.

    Parameters
    ----------
    filename : str
        Path to file where the singles energies were saved.

    Returns
    -------
    numpy.ndarray
        Array of energies where each row corresponds to the MO energies of each
        chromophore in the list.
    """
    return np.loadtxt(filename)


def set_energyvalues(chromo_list, s_filename, d_filename):
    """Set the energy attributes of the Chromophore objects in chromo_list.

    Run singles_homolumo and dimer_homolumo first to get the energy files.
    Energy values set by this function:
        homo_1, homo, lumo, lumo_1, neighbors_delta_e, neighbors_ti

    Parameters
    ----------
    chromo_list : list of Chromophore
        Set the energy values of the chromphores in this list.
    s_filename : str
        Path to file where the singles energies were saved.
    d_filename : str
        Path to file where the pair energies were saved.
    """
    s_data = get_singlesdata(s_filename)
    d_data = get_dimerdata(d_filename)

    for i, chromo in enumerate(chromo_list):
        chromo.homo_1, chromo.homo, chromo.lumo, chromo.lumo_1 = s_data[i]

    for (i, j), (homo_1, homo, lumo, lumo_1) in d_data:
        ichromo = chromo_list[i]
        jchromo = chromo_list[j]
        ineighborind = [i for i, img in ichromo.neighbors].index(j)
        jneighborind = [i for i, img in jchromo.neighbors].index(i)
        deltaE = ti.calculate_delta_E(ichromo, jchromo)
        ichromo.neighbors_delta_e[ineighborind] = deltaE
        jchromo.neighbors_delta_e[jneighborind] = -deltaE

        assert ichromo.species == jchromo.species
        if ichromo.species == "donor":
            transint = ti.calculate_ti(homo - homo_1, deltaE)
        else:
            transint = ti.calculate_ti(lumo - lumo_1, deltaE)
        ichromo.neighbors_ti[ineighborind] = transint
        jchromo.neighbors_ti[jneighborind] = transint


def write_qcc_inp(snap, atom_ids, conversion_dict=None):
    """Write a quantum chemical input string.

    Input string for pySCF containing elements and positions in Angstroms
    (e.g., "C 0.0 0.0 0.0; H 1.54 0.0 0.0")
    See https://pyscf.org/quickstart.html for more information.

    Parameters
    ----------
    snap : gsd.hoomd.Snapshot
        Atomistic simulation snapshot from a GSD file. It is expected that the
        lengths in this file have been converted to Angstroms.
    atom_ids : numpy.ndarray of int
        Snapshot indices of the particles to include in the input string.
    conversion_dict : dictionary, default None
        A dictionary that maps the atom type to its element. e.g., `{'c3': C}`.
        An instance that maps AMBER types to their element can be found in
        `amber_dict`. If None is given, assume the particles already have
        element names.

    Returns
    -------
    str
        The input for the MINDO3 quantum chemical calculation run in pySCF.
    """
    atoms = []
    positions = []

    box = snap.configuration.box[:3]
    unwrapped_pos = snap.particles.position + snap.particles.image * box

    for i in atom_ids:
        if conversion_dict is not None:
            element = conversion_dict[
                snap.particles.types[snap.particles.typeid[i]]
            ]
        else:
            element = ele.element_from_symbol(
                snap.particles.types[snap.particles.typeid[i]]
            )
        atoms.append(element.symbol)
        positions.append(unwrapped_pos[i])

    # To determine where to add hydrogens, check the bonds that go to
    # particles outside of the ids provided
    for i, j in snap.bonds.group:
        if i in atom_ids and j not in atom_ids:
            if conversion_dict is not None:
                element = conversion_dict[
                    snap.particles.types[snap.particles.typeid[j]]
                ]
            else:
                element = ele.element_from_symbol(
                    snap.particles.types[snap.particles.typeid[j]]
                )
            # If it's already a Hydrogen, just add it
            if element.atomic_number == 1:
                atoms.append(element.symbol)
                positions.append(unwrapped_pos[j])
            # If it's not a hydrogen, use the existing bond vector to
            # determine the direction and scale it to a more reasonable
            # length for C-H bond
            else:
                # Average sp3 C-H bond is 1.094 Angstrom
                v = unwrapped_pos[j] - unwrapped_pos[i]
                unit_vec = v / np.linalg.norm(v)
                new_pos = unit_vec * 1.094 + unwrapped_pos[i]
                atoms.append("H")
                positions.append(new_pos)

        # Same as above but j->i instead of i->j
        elif j in atom_ids and i not in atom_ids:
            if conversion_dict is not None:
                element = conversion_dict[
                    snap.particles.types[snap.particles.typeid[i]]
                ]
            else:
                element = ele.element_from_symbol(
                    snap.particles.types[snap.particles.typeid[i]]
                )
            if element.atomic_number == 1:
                atoms.append(element.symbol)
                positions.append(unwrapped_pos[i])

            else:
                v = unwrapped_pos[i] - unwrapped_pos[j]
                unit_vec = v / np.linalg.norm(v)
                new_pos = unit_vec * 1.094 + unwrapped_pos[j]
                atoms.append("H")
                positions.append(new_pos)

    # Shift center to origin
    positions = np.stack(positions)
    positions -= np.mean(positions, axis=0)
    qcc_input = " ".join(
        [f"{atom} {x} {y} {z};" for atom, (x, y, z) in zip(atoms, positions)]
    )
    return qcc_input


def write_qcc_pair_input(
    snap, chromo_i, chromo_j, j_shift, conversion_dict=None
    ):
    """Write a quantum chemical input string for chromophore pairs.

    Pair input requires taking periodic images into account.
    Input string for pySCF containing elements and positions in Angstroms
    (e.g., "C 0.0 0.0 0.0; H 1.54 0.0 0.0")
    See https://pyscf.org/quickstart.html for more information.

    Parameters
    ----------
    snap : gsd.hoomd.Snapshot
        Atomistic simulation snapshot from a GSD file. It is expected that the
        lengths in this file have been converted to Angstroms.
    chromo_i : Chromophore
        One of the chromophores to be written.
    chromo_j : Chromophore
        One of the chromophores to be written.
    j_shift : numpy.ndarray(3)
        Vector to shift chromo_j.
        (chromo_j minimum image center - unwrapped center)
    conversion_dict : dictionary, default None
        A dictionary that maps the atom type to its element. e.g., `{'c3': C}`.
        An instance that maps AMBER types to their element can be found in
        `amber_dict`. If None is given, assume the particles already have
        element names.

    Returns
    -------
    str
        The input for the MINDO3 quantum chemical calculation run in pySCF.
    """
    box = snap.configuration.box[:3]
    unwrapped_pos = snap.particles.position + snap.particles.image * box

    # chromophore i is shifted into 0,0,0 image
    positions = [
        i + chromo_i.image * box for i in unwrapped_pos[chromo_i.atom_ids]
    ]
    # shift chromophore j's unwrapped positions
    positions += [i for i in unwrapped_pos[chromo_j.atom_ids] + j_shift]

    atom_ids = np.concatenate((chromo_i.atom_ids, chromo_j.atom_ids))
    typeids = snap.particles.typeid[atom_ids]
    if conversion_dict is not None:
        atoms = [
            conversion_dict[snap.particles.types[i]].symbol for i in typeids
        ]
    else:
        atoms = [
            ele.element_from_symbol(snap.particles.types[i]) for i in typeids
        ]

    # To determine where to add hydrogens, check the bonds that go to
    # particles outside of the ids provided
    for i, j in snap.bonds.group:
        if i in atom_ids and j not in atom_ids:
            # If bond is to chromophore j, additional shifting might be needed
            if i in chromo_j.atom_ids:
                shift = j_shift
            else:
                shift = chromo_i.image * box
            if conversion_dict is not None:
                element = conversion_dict[
                    snap.particles.types[snap.particles.typeid[j]]
                ]
            else:
                element = ele.element_from_symbol(
                    snap.particles.types[snap.particles.typeid[j]]
                )
            # If it's already a Hydrogen, just add it
            if element.atomic_number == 1:
                atoms.append(element.symbol)
                positions.append(unwrapped_pos[j] + shift)
            # If it's not a hydrogen, use the existing bond vector to
            # determine the direction and scale it to a more reasonable
            # length for C-H bond
            else:
                # Average sp3 C-H bond is 1.094 Angstrom
                v = unwrapped_pos[j] - unwrapped_pos[i]
                unit_vec = v / np.linalg.norm(v)
                new_pos = unit_vec * 1.094 + unwrapped_pos[i] + shift
                atoms.append("H")
                positions.append(new_pos)

        # Same as above but j->i instead of i->j
        elif j in atom_ids and i not in atom_ids:
            if j in chromo_j.atom_ids:
                shift = j_shift
            else:
                shift = chromo_i.image * box
            if conversion_dict is not None:
                element = conversion_dict[
                    snap.particles.types[snap.particles.typeid[i]]
                ]
            else:
                element = ele.element_from_symbol(
                    snap.particles.types[snap.particles.typeid[i]]
                )

            if element.atomic_number == 1:
                atoms.append(element.symbol)
                positions.append(unwrapped_pos[i] + shift)
            else:
                v = unwrapped_pos[i] - unwrapped_pos[j]
                unit_vec = v / np.linalg.norm(v)
                new_pos = unit_vec * 1.094 + unwrapped_pos[j] + shift
                atoms.append("H")
                positions.append(new_pos)

    # Shift center to origin
    positions = np.stack(positions)
    positions -= np.mean(positions, axis=0)

    qcc_input = " ".join(
        [f"{atom} {x} {y} {z};" for atom, (x, y, z) in zip(atoms, positions)]
    )
    return qcc_input


def _worker_wrapper(arg):
    qcc_input, charge = arg
    return get_homolumo(qcc_input, charge=charge)
