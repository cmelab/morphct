import os
import time

import gsd, gsd.hoomd
import numpy as np

from morphct.chromophores import Chromophore, set_neighbors_voronoi
from morphct.execute_qcc import (
    singles_homolumo, dimer_homolumo, set_energyvalues
)
from morphct.mobility_kmc import run_kmc
from morphct import kmc_analyze


class System():  # pragma: no cover
    """An object for managing all data for the KMC simulation.

    Parameters
    ----------
    gsdfile : path
        The path to a gsd file containing the structure for the KMC simulation.
    outpath : path
        The path to a directory where output files will be saved. If the path
        does not exist, it will be created.
    frame : int
        The frame number of the gsdfile to use.
    scale : float
        Scaling factor to convert the lengthscale in the gsdfile to Angstrom.
    conversion_dict : dict
        A dictionary to map atom types to an ele.element. An example for mapping
        the General AMBER forcefield (GAFF) types to elements can be found in
        `morphct.chromophores.amber_dict`.

    Attributes
    ----------
    snap : gsd.hoomd.Snapshot
        The system snapshot with lengths scaled to Angstroms.
    conversion_dict : dict
        A dictionary to map atom types to an ele.element.
    chromophores : list of Chromophore
        List of chromphores in the simulation.
    outpath : path
        The path to a directory where output files will be saved.
    qcc_pairs : list of ((int, int), str)
        QCC input for the pairs. Each list item contains a tuple of the pair
        indices and the QCC input string.

    Methods
    -------
    add_chromophores
        Add chromophore(s) to the system.
    compute_energies
        Compute the energies of the chromophores in the system.
    set_energies
        Set the computed energies.
    run_kmc
        Run the KMC simulation.
    visualize_qcc_input
        Visualize the input to QCC.
    visualize_system
        Visualize the entire system.
    visualize_chromophores
        Visualize the chromophores.

    Note
    ----
    The visualize functions require mbuild, a visualization backend (py3Dmol or
    nglview), and a jupyter notebook. They are designed only to help
    troubleshoot.
    """

    def __init__(
        self, gsdfile, outpath, frame=-1, scale=1.0, conversion_dict=None,
    ):
        with gsd.hoomd.open(name=gsdfile, mode="rb") as f:
            snap = f[frame]

        # It is expected that the snapshot is in the Angstrom length scale.
        # If not, the scaling factor is used to adjust
        snap.particles.position *= scale
        snap.configuration.box[:3] *= scale

        self.snap = snap
        self.conversion_dict = conversion_dict
        if not os.path.exists(outpath):
            os.makedirs(outpath)
        self.outpath = outpath
        self._chromophores = []
        self._comp = None
        self.qcc_pairs = None
        self._dinds = []
        self._ainds = []

    @property
    def chromophores(self):
        """Return the chromophores in the system."""
        return self._chromophores

    def add_chromophores(self, indices, species, chromophore_kwargs={}):
        """Add chromophore(s) to the system.

        Parameters
        ----------
        indices : list of numpy.ndarray,
            Atom indices in the snapshot for this chromophore type. Each array
            will be treated as a separate chromophore instance.
        species : str
            Chromophore species ("donor" or "acceptor").
        chromophore_kwargs : dict, default {}
            Additional keywrod arguments to be passed to the Chromophore class.
        """
        start = len(self.chromophores)
        for i, ind in enumerate(indices):
            self._chromophores.append(
                Chromophore(
                    i+start,
                    self.snap,
                    ind,
                    species,
                    self.conversion_dict,
                    **chromophore_kwargs
                )
            )
        if species == "acceptor":
            self._ainds += indices
        else:
            self._dinds += indices

    def compute_energies(self, dcut=None):
        """Compute the energies of the chromophores in the system.

        Parameters
        ----------
        dcut : float, default None
            The distance cutoff for chromophore neighbors. If None is provided,
            the cutoff will be set to half the smallest box length of the
            snapshot.
        """
        if dcut is None:
            dcut = min(self.snap.configuration.box[:3]/2)
        self.qcc_pairs = set_neighbors_voronoi(
            self.chromophores, self.snap, self.conversion_dict, d_cut=dcut
        )
        print(f"There are {len(self.qcc_pairs)} chromophore pairs")

        s_filename = os.path.join(self.outpath, "singles_energies.txt")
        d_filename = os.path.join(self.outpath, "dimer_energies.txt")

        t0 = time.perf_counter()
        print("Starting singles energy calculation...")
        data = singles_homolumo(self.chromophores, s_filename)
        t1 = time.perf_counter()
        print(f"Finished in {t1-t0:.2f} s. Output written to {s_filename}.")

        print("Starting dimer energy calculation...")
        dimer_data = dimer_homolumo(
            self.qcc_pairs, self.chromophores, d_filename
        )
        t2 = time.perf_counter()
        print(f"Finished in {t2-t1:.2f} s. Output written to {d_filename}.")

    def set_energies(self, dcut=None):
        """Set the computed energies.

        Parameters
        ----------
        dcut : float, default None
            The distance cutoff for chromophore neighbors. If None is provided,
            the cutoff will be set to half the smallest box length of the
            snapshot.
        """
        if self.qcc_pairs is None:
            if dcut is None:
                dcut = min(self.snap.configuration.box[:3]/2)
            self.qcc_pairs = set_neighbors_voronoi(
                self.chromophores, self.snap, self.conversion_dict, d_cut=dcut
            )

        s_filename = os.path.join(self.outpath, "singles_energies.txt")
        d_filename = os.path.join(self.outpath, "dimer_energies.txt")

        if not (os.path.isfile(s_filename) and os.path.isfile(d_filename)):
            raise FileNotFoundError(
                f"Expected to find {s_filename} and {d_filename}, but didn't."
            )

        set_energyvalues(self.chromophores, s_filename, d_filename)
        print("Energies set.")

    def run_kmc(
        self,
        lifetimes,
        temp,
        n_holes=0,
        n_elec=0,
        seed=42,
        combine=True,
        carrier_kwargs={},
        verbose=0
    ):
        """Run the KMC simulation.

        Parameters
        ----------
        lifetimes : list of float
            The potential lifetimes of the carriers. A value from these will be
            randomly assigned to each run.
        temp : float
            The simulation temperature in Kelvin.
        n_holes : int, default 0
            The number of holes to simulate.
        n_elec : int, default 0
            The number of electrons to simulate.
        seed : int, default 42
            A seed for the random processes.
        carrier_kwargs : dict, default {}
            Additional keyword arguments to be passed to the carrier instances.
        verbose : int, default 0
            The verbosity level of output.
        """
        kmc_dir = os.path.join(self.outpath, "kmc")
        if not os.path.exists(kmc_dir):
            os.makedirs(kmc_dir)
        data = run_kmc(
            lifetimes,
            kmc_dir,
            self.chromophores,
            self.snap,
            temp,
            n_holes=n_holes,
            n_elec=n_elec,
            seed=seed,
            combine=combine,
            carrier_kwargs=carrier_kwargs,
            verbose=verbose
        )

        kmc_analyze.main(data, temp, self.chromophores, self.snap, kmc_dir)

    def visualize_qcc_input(self, i, single=True):
        """Visualize the input for pyscf using mbuild.

        Parameters
        ----------
        i : int
            Index of the single or pair input to visualize.
        single : bool, default True
            Whether to visualize a single or pair input.

        Note
        ----
        The visualize functions require mbuild, a visualization backend (py3Dmol
        or nglview), and a jupyter notebook. They are designed only to help
        troubleshoot.
        """
        import mbuild as mb

        if single:
            qcc_input = self.chromophores[i].qcc_input
        else:
            qcc_input = self.qcc_pairs[i][1]
        comp = mb.Compound()
        for line in qcc_input.split(";")[:-1]:
            atom, x, y, z = line.split()
            xyz = np.array([x,y,z], dtype=float)
            # Angstrom -> nm
            xyz /= 10
            comp.add(mb.Particle(name=atom, pos=xyz))
        comp.visualize().show()

    def visualize_system(self):
        """Visualize the system snapshot.

        Note
        ----
        The visualize functions require mbuild, a visualization backend (py3Dmol
        or nglview), and a jupyter notebook. They are designed only to help
        troubleshoot.
        """
        if self._comp is None:
            self._comp = self._make_comp()
        self._comp.visualize().show()

    def visualize_chromophores(self):
        """Visualize the chromophores.

        Note
        ----
        The visualize functions require mbuild, a visualization backend (py3Dmol
        or nglview), and a jupyter notebook. They are designed only to help
        troubleshoot.
        """
        import mbuild as mb

        if self._comp is None:
            self._comp = self._make_comp()
        ccomp = mb.clone(self._comp)
        for i,p in enumerate(ccomp):
            if self._ainds:
                if i in np.hstack(self._ainds):
                    p.name = "Br"
            if self._dinds:
                if i in np.hstack(self._dinds):
                    p.name = "I"
        ccomp.visualize().show()

    def _make_comp(self):
        """Convert the system Snapshot to an mbuild Compound.

        Returns
        -------
        comp : mb.Compound
        """
        import mbuild as mb

        comp = mb.Compound()
        snap = self.snap
        bond_array = snap.bonds.group
        n_atoms = snap.particles.N

        box = snap.configuration.box[:3]
        # to_hoomdsnapshot shifts the coords, this will keep consistent
        shift = box/2

        unwrap_pos = snap.particles.position + snap.particles.image * box
        # angstrom -> nm
        unwrap_pos *= 0.1
        unwrap_pos += shift
        types = [
            self.conversion_dict[i].symbol for i in snap.particles.types
        ]

        # Add particles
        for i in range(n_atoms):
            name = types[snap.particles.typeid[i]]
            xyz = unwrap_pos[i]
            atom = mb.Particle(name=name, pos=xyz)
            comp.add(atom, label=str(i))

        # Add bonds
        for i in range(bond_array.shape[0]):
            atom1 = int(bond_array[i][0])
            atom2 = int(bond_array[i][1])
            comp.add_bond([comp[atom1], comp[atom2]])

        return comp
