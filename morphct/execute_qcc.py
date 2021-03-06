import multiprocessing as mp
from multiprocessing import get_context
import numpy as np

import pyscf
from pyscf.semiempirical import MINDO3

from morphct import helper_functions as hf
from morphct import transfer_integrals as ti


def get_homolumo(molstr, verbose=0, tol=1e-6):
    """
    Get the HOMO-1, HOMO, LUMO, LUMO+1 energies for the given input string
    using the MINDO3 method in pySCF.

    Parameters
    ----------
    molstr : str
        Input string for pySCF containing elements and positions in Angstroms
        (e.g., "C 0.0 0.0 0.0; H 1.54 0.0 0.0")
    verbose : int
        Verbosity level of the MINDO calculation output. 0 will silence output,
        4 will show convergence. (default 0)
    tol : float
        Tolerance of the MINDO convergence (default 1e-6)

    Returns
    -------
    numpy.array
        Array containing HOMO-1, HOMO, LUMO, LUMO+1 energies in eV
    """
    mol = pyscf.M(atom=molstr)
    mf = MINDO3(mol).run(verbose=verbose, conv_tol=tol)
    occ = mf.get_occ()
    i_lumo = np.argmax(occ<1)
    energies = mf.mo_energy[i_lumo-2:i_lumo+2]
    energies *= 27.2114 # convert Eh to eV
    return energies


def singles_homolumo(chromo_list, filename=None, nprocs=None):
    """
    Obtain the HOMO-1, HOMO, LUMO, LUMO+1 energies in eV for all single
    chromophores

    Parameters
    ----------
    chromo_list : list of morphct.obtain_chromophores.Chromophore
        Chromophores to calculate energies of. Each Chromophore must have
        qcc_input attribute set.
    filename : str
        Path to file where singles energies will be saved. If None, energies
        will not be saved. (default None)
    nprocs : int
        Number of processes passed to multiprocessing.Pool. If None,
        multiprocessing.cpu_count will used to determine optimal number.
        (default None)

    Returns
    -------
    data : numpy.array
        Array of energies where each row corresponds to the MO energies of each
        chromophore in the list.
    """
    if nprocs is None:
        nprocs = mp.cpu_count()
    with get_context("spawn").Pool(processes=nprocs) as p:
        data = p.map(get_homolumo, [i.qcc_input for i in chromo_list])

    data = np.stack(data)
    if filename is not None:
        np.savetxt(filename, data)
    return data


def dimer_homolumo(qcc_pairs, filename=None, nprocs=None):
    """
    Obtain the HOMO-1, HOMO, LUMO, LUMO+1 energies in eV for all chromophore
    pairs

    Parameters
    ----------
    qcc_pairs : list of (tuple and string)
        Each list item contains a tuple with the indices of the pair and the
        qcc input string.
        qcc_pairs is returned by morphct.chromophores.set_neighbors_voronoi()
    filename : str
        Path to file where the pair energies will be saved. If None, energies
        will not be saved. (default None)
    nprocs : int
        Number of processes passed to multiprocessing.Pool. If None,
        multiprocessing.cpu_count will used to determine optimal number.
        (default None)

    Returns
    -------
    dimer_data : list of (tuple or ints, numpy.array)
        Each list item contains the indices of the pair and an array of its MO
        energies.
    """
    if nprocs is None:
        nprocs = mp.cpu_count()

    with get_context("spawn").Pool(processes=nprocs) as p:
        data = p.map(get_homolumo, [qcc_input for pair,qcc_input in qcc_pairs])

    dimer_data = [i for i in zip([pair for pair,qcc_input in qcc_pairs],data)]
    if filename is not None:
        with open(filename, "w") as f:
            f.writelines(
                f"{pair[0]} {pair[1]} {en[0]} {en[1]} {en[2]} {en[3]}\n"
                for pair,en in dimer_data
            )
    return dimer_data


def get_dimerdata(filename):
    """
    Read in the saved data created by dimer_homolumo().

    Parameters
    ----------
    filename : str
        Path to file where the dimer energies were saved.

    Returns
    -------
    dimer_data : list of (tuple or ints, numpy.array)
        Each list item contains the indices of the pair and an array of its MO
        energies.
    """
    dimer_data = []
    with open(filename, "r") as f:
        for i in f.readlines():
            a,b,c,d,e,f = i.split()
            dimer_data.append(
                ((int(a),int(b)),(float(c),float(d),float(e),float(f)))
            )
    return dimer_data


def get_singlesdata(filename):
    """
    Read in the saved data created by singles_homolumo().

    Parameters
    ----------
    filename : str
        Path to file where the singles energies were saved.

    Returns
    -------
    numpy.array
        Array of energies where each row corresponds to the MO energies of each
        chromophore in the list.
    """
    return np.loadtxt(filename)


def set_energyvalues(chromo_list, s_filename, d_filename):
    """
    Set the energy attributes of the Chromophore objects in chromo_list.
    Run singles_homolumo and dimer_homolumo first to get the energy files.
    Energy values set by this function:
        homo_1, homo, lumo, lumo_1, neighbors_delta_E, neighbors_TI

    Parameters
    ----------
    chromo_list : list of morphct.obtain_chromophores.Chromophore
        Set the energy values of the chromphores in this list.
    s_filename : str
        Path to file where the singles energies were saved.
    d_filename : str
        Path to file where the pair energies were saved.
    """
    s_data = get_singlesdata(s_filename)
    d_data = get_dimerdata(d_filename)

    for i,chromo in enumerate(chromo_list):
        chromo.homo_1, chromo.homo, chromo.lumo, chromo.lumo_1 = s_data[i]

    for (i,j), (homo_1, homo, lumo, lumo_1) in d_data:
        chromo1 = chromo_list[i]
        chromo2 = chromo_list[j]
        neighborind1 = [i[0] for i in chromo1.neighbors].index(j)
        neighborind2 = [i[0] for i in chromo2.neighbors].index(i)
        deltaE = ti.calculate_delta_E(chromo1,chromo2)
        chromo1.neighbors_delta_e[neighborind1] = deltaE
        chromo2.neighbors_delta_e[neighborind2] = -deltaE

        assert chromo1.species == chromo2.species
        if chromo1.species == "donor":
            TI = ti.calculate_TI(homo - homo_1, deltaE)
        else:
            TI = ti.calculate_TI(lumo - lumo_1, deltaE)
        chromo1.neighbors_ti[neighborind1] = TI
        chromo2.neighbors_ti[neighborind2] = TI


def write_qcc_inp(snap, atom_ids, conversion_dict):
    """

    Parameters
    ----------

    Returns
    -------

    """
    atoms = []
    positions = []

    box = snap.configuration.box[:3]
    unwrapped_pos = snap.particles.position + snap.particles.image * box

    for i in atom_ids:
        element = conversion_dict[
                snap.particles.types[snap.particles.typeid[i]]
                ]
        atoms.append(element.symbol)
        positions.append(unwrapped_pos[i])

    # To determine where to add hydrogens, check the bonds that go to
    # particles outside of the ids provided
    for i,j in snap.bonds.group:
        if i in atom_ids and j not in atom_ids:
            element = conversion_dict[
                    snap.particles.types[snap.particles.typeid[j]]
                    ]
            # If it's already a Hydrogen, just add it
            if element.atomic_number == 1:
                atoms.append(element.symbol)
                positions.append(unwrapped_pos[j])
            # If it's not a hydrogen, use the existing bond vector to
            # determine the direction and scale it to a more reasonable
            # length for C-H bond
            else:
                # Average sp3 C-H bond is 1.094 Angstrom
                v = unwrapped_pos[j]-unwrapped_pos[i]
                unit_vec = v/np.linalg.norm(v)
                new_pos = unit_vec * 1.094 + unwrapped_pos[i]
                atoms.append("H")
                positions.append(new_pos)

        # Same as above but j->i instead of i->j
        elif j in atom_ids and i not in atom_ids:
            element = conversion_dict[
                    snap.particles.types[snap.particles.typeid[i]]
                    ]
            if element.atomic_number == 1:
                atoms.append(element.symbol)
                positions.append(unwrapped_pos[i])

            else:
                v = unwrapped_pos[i]-unwrapped_pos[j]
                unit_vec = v/np.linalg.norm(v)
                new_pos = unit_vec * 1.094 + unwrapped_pos[j]
                atoms.append("H")
                positions.append(new_pos)

    # Shift center to origin
    positions = np.stack(positions)
    positions -= np.mean(positions,axis=0)
    qcc_input = " ".join(
            [f"{atom} {x} {y} {z};" for atom,(x,y,z) in zip(atoms,positions)]
            )
    return qcc_input


def write_qcc_pair_input(snap, chromo_i, chromo_j, j_shift, conversion_dict):
    """

    Parameters
    ----------
    j_shift : numpy.array(3)
        vector to shift chromophore j
        (chromophore j minimum image center - unwrapped center)

    Returns
    -------

    """
    box = snap.configuration.box[:3]
    unwrapped_pos = snap.particles.position + snap.particles.image * box

    # chromophore i is shifted into 0,0,0 image
    positions = [
            i+chromo_i.image*box for i in unwrapped_pos[chromo_i.atom_ids]
            ]
    # shift chromophore j's unwrapped positions
    positions += [i for i in unwrapped_pos[chromo_j.atom_ids] + j_shift]

    atom_ids = np.concatenate((chromo_i.atom_ids, chromo_j.atom_ids))
    typeids = snap.particles.typeid[atom_ids]
    atoms = [conversion_dict[snap.particles.types[i]].symbol for i in typeids]

    # To determine where to add hydrogens, check the bonds that go to
    # particles outside of the ids provided
    for i,j in snap.bonds.group:
        if i in atom_ids and j not in atom_ids:
            # If bond is to chromophore j, additional shifting might be needed
            if i in chromo_j.atom_ids:
                shift = j_shift
            else:
                shift = chromo_i.image*box
            element = conversion_dict[
                    snap.particles.types[snap.particles.typeid[j]]
                    ]
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
                unit_vec = v/np.linalg.norm(v)
                new_pos = unit_vec * 1.094 + unwrapped_pos[i] + shift
                atoms.append("H")
                positions.append(new_pos)

        # Same as above but j->i instead of i->j
        elif j in atom_ids and i not in atom_ids:
            if j in chromo_j.atom_ids:
                shift = j_shift
            else:
                shift = chromo_i.image*box
            element = conversion_dict[
                    snap.particles.types[snap.particles.typeid[i]]
                    ]

            if element.atomic_number == 1:
                atoms.append(element.symbol)
                positions.append(unwrapped_pos[i] + shift)
            else:
                v = unwrapped_pos[i] - unwrapped_pos[j]
                unit_vec = v/np.linalg.norm(v)
                new_pos = unit_vec * 1.094 + unwrapped_pos[j] + shift
                atoms.append("H")
                positions.append(new_pos)

    # Shift center to origin
    positions = np.stack(positions)
    positions -= np.mean(positions,axis=0)

    qcc_input = " ".join(
            [f"{atom} {x} {y} {z};" for atom,(x,y,z) in zip(atoms,positions)]
            )
    return qcc_input
