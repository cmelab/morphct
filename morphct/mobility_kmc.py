import math
import multiprocessing as mp
import os
from sys import platform
import time
import warnings

import freud
import numpy as np
from scipy.sparse import lil_matrix

from morphct import helper_functions as hf
from morphct.helper_functions import v_print


try: # pragma: no cover
    if platform == "darwin":
        # OS X
        mp.set_start_method("fork")
except RuntimeError:
    pass


class Carrier:
    def __init__(
        self,
        chromo,
        lifetime,
        carrier_no,
        box,
        temp,
        n,
        hop_limit=None,
        record_history=True,
        mol_id_dict=None,
        use_avg_hoprates=False,
        avg_intra_rate=None,
        avg_inter_rate=None,
        use_koopmans=False,
        boltz=False,
        use_vrh=False,
        hopping_prefactor=1.0,
    ):
        both_rates = avg_inter_rate is None and avg_intra_rate is None
        any_rate = avg_inter_rate is None or avg_intra_rate is None
        if use_avg_hoprates and any_rate:
            err_msg = (
                    "If use_avg_hoprates is True, avg_inter_rate and "+
                    "avg_intra_rate must also be provided"
                    )
            raise ValueError(err_msg)

        elif not use_avg_hoprates and not both_rates:
            warn_msg = (
                    "You provided an average hop rate (avg_inter_rate or "+
                    "avg_intra_rate) but use_avg_hoprates is False. "+
                    "The provided hop rate will not be used"
                    )
            warnings.warn(warn_msg)

        elif use_avg_hoprates and mol_id_dict is None:
            err_msg = (
                    "If use_avg_hoprates is True, a molecule dictionary "+
                    "(mol_id_dict) must also be provided"
                    )
            raise ValueError(err_msg)

        self.id = carrier_no
        self.image = np.array([0, 0, 0])
        self.initial_chromo = chromo
        self.current_chromo = chromo
        self.lambda_ij = self.current_chromo.reorganization_energy
        self.hop_limit = hop_limit
        self.temp = temp
        self.lifetime = lifetime
        self.current_time = 0.0

        self.hole_history = None
        self.electron_history = None
        if self.current_chromo.species == "donor":
            self.c_type = "hole"
            if record_history:
                self.hole_history = lil_matrix((n, n), dtype=int)
        elif self.current_chromo.species == "acceptor":
            self.c_type = "electron"
            if record_history:
                self.electron_history = lil_matrix((n, n), dtype=int)

        self.n_hops = 0
        self.box = box
        self.displacement = 0
        self.mol_id_dict = mol_id_dict

        self.use_avg_hoprates = use_avg_hoprates
        self.avg_intra_rate = avg_intra_rate
        self.avg_inter_rate = avg_inter_rate

        # Set the use of Koopmans' approximation to false if the key does not
        # exist in the parameter dict
        self.use_koopmans = use_koopmans
        # Are we using a simple Boltzmann penalty?
        self.boltz = boltz
        # Are we applying a distance penalty beyond the transfer integral?
        self.use_vrh = use_vrh
        if self.use_vrh:
            self.vrh_delocalization = self.current_chromo.vrh_delocalization

        self.hopping_prefactor = hopping_prefactor

    def update_displacement(self):
        init_pos = self.initial_chromo.center
        final_pos = self.current_chromo.center
        displacement = final_pos - init_pos + self.image * self.box
        self.displacement = np.linalg.norm(displacement)

    def calculate_hop(self, chromo_list, verbose=0):
        # Terminate if the next hop would be more than the termination limit
        if self.hop_limit is not None:
            if self.n_hops + 1 > self.hop_limit:
                return False
        # Determine the hop times to all possible neighbors
        hop_times = []
        if self.use_avg_hoprates:
            # Use the average hop values given in the parameter dict to pick a
            # hop
            for i, img in self.current_chromo.neighbors:
                neighbor = chromo_list[i]
                assert neighbor.id == i

                current_mol = self.mol_id_dict[self.current_chromo.id]
                neighbor_mol = self.mol_id_dict[neighbor.id]
                if current_mol == neighbor_mol:
                    hop_rate = self.avg_intra_rate
                else:
                    hop_rate = self.avg_inter_rate
                hop_time = hf.get_event_tau(hop_rate)
                # Keep track of the chromophoreid and the corresponding tau
                hop_times.append([neighbor.id, hop_time, img])
        else:
            # Obtain the reorganization energy in J (from eV)
            for i_neighbor, ti in enumerate(self.current_chromo.neighbors_ti):
                # Ignore any hops with a NoneType transfer integral
                if ti is None:
                    continue
                # index of the neighbor in the main list
                n_ind = self.current_chromo.neighbors[i_neighbor][0]
                delta_E_ij = self.current_chromo.neighbors_delta_e[i_neighbor]
                # Load the specified hopping prefactor
                prefactor = self.hopping_prefactor
                # Get the relative image so we can update the carrier image
                # after the hop
                rel_img = self.current_chromo.neighbors[i_neighbor][1]
                # All of the energies are in eV currently, so convert them to J
                if self.use_vrh is True:
                    neighbor_chromo = chromo_list[n_ind]
                    neighbor_pos = neighbor_chromo.center + rel_img * self.box

                    # Chromophore separation needs converting to m
                    sep = (
                        np.linalg.norm(
                            self.current_chromo.center - neighbor_pos
                        )
                        * 1e-10
                    )

                    hop_rate = hf.get_hop_rate(
                        self.lambda_ij,
                        ti,
                        delta_E_ij,
                        prefactor,
                        self.temp,
                        use_vrh=True,
                        rij=sep,
                        vrh=self.vrh_delocalization,
                        boltz=self.boltz,
                    )
                else:
                    hop_rate = hf.get_hop_rate(
                        self.lambda_ij,
                        ti,
                        delta_E_ij,
                        prefactor,
                        self.temp,
                        boltz=self.boltz,
                    )
                hop_time = hf.get_event_tau(hop_rate)
                # Keep track of the chromophoreid and the corresponding tau
                hop_times.append([n_ind, hop_time, rel_img])
        # Sort by ascending hop time
        hop_times.sort(key=lambda x: x[1])
        if len(hop_times) == 0:
            # We are trapped here, so create a dummy hop with time 1E99
            hop_times = [[self.current_chromo.id, 1e99, [0, 0, 0]]]
        # As long as we're not limiting by the number of hops:
        if self.hop_limit is None:
            # Ensure that the next hop does not put the carrier over its
            # lifetime
            if (self.current_time + hop_times[0][1]) > self.lifetime:
                # Send the termination signal to singleCoreRunKMC.py
                return False
        # Move the carrier and send the contiuation signal to
        # singleCoreRunKMC.py
        # Take the quickest hop
        n_ind, hop_time, rel_img = hop_times[0]

        v_print("\thop_times:", verbose, v_level=1)
        hop_str = "\n".join([f"\t\t{i} {j:.2e} {k}" for (i, j, k) in hop_times])
        v_print(hop_str, verbose, v_level=1)
        v_print(f"\tHopping to {n_ind}", verbose, v_level=1)

        self.perform_hop(chromo_list[n_ind], hop_time, rel_img)
        return True

    def perform_hop(self, destination_chromo, hop_time, rel_image):
        init_id = self.current_chromo.id
        dest_id = destination_chromo.id
        self.image += rel_image
        # Carrier image now sorted, so update its current position
        self.current_chromo = destination_chromo
        # Increment the simulation time
        self.current_time += hop_time
        # Increment the hop counter
        self.n_hops += 1
        # Now update the sparse history matrix
        if self.c_type == "hole" and self.hole_history is not None:
            self.hole_history[init_id, dest_id] += 1
        elif self.c_type == "electron" and self.electron_history is not None:
            self.electron_history[init_id, dest_id] += 1


def run_single_kmc(
    jobs,
    KMC_directory,
    chromo_list,
    snap,
    temp,
    carrier_kwargs={},
    cpu_rank=None,
    seed=None,
    send_end=None,
    verbose=0,
):
    if seed is not None:
        np.random.seed(seed)

    if cpu_rank is not None:
        # If we're running on multiple cpus, don't print to std out
        # print to a log file instead and start fresh: remove it if it exists
        filename = os.path.join(KMC_directory, f"kmc_{cpu_rank:02d}.log")
        if os.path.exists(filename):
            os.remove(filename)
    else:
        filename = None

    v_print(f"Found {len(jobs):d} jobs to run", verbose, filename=filename)

    try:
        use_avg_hoprates = carrier_kwargs["use_avg_hoprates"]
    except KeyError:
        use_avg_hoprates = False

    if use_avg_hoprates:
        # Chosen to split hopping by inter-intra molecular hops, so get
        # molecule data
        mol_id_dict = get_molecule_ids(snap, chromo_list)
        # molidDict is a dictionary where the keys are the chromoids, and
        # the vals are the molids
    else:
        mol_id_dict = None
    t0 = time.perf_counter()
    carrier_list = []
    box = snap.configuration.box[:3]
    for i_job, [carrier_no, lifetime, ctype] in enumerate(jobs):
        v_print(f"starting job {i_job}", verbose, filename=filename)
        t1 = time.perf_counter()
        # Find a random position to start the carrier in
        while True:
            i = np.random.randint(0, len(chromo_list) - 1)
            i_chromo = chromo_list[i]
            if (ctype == "electron") and (i_chromo.species != "acceptor"):
                continue
            elif (ctype == "hole") and (i_chromo.species != "donor"):
                continue
            break
        # Create the carrier instance
        i_carrier = Carrier(
            i_chromo,
            lifetime,
            carrier_no,
            box,
            temp,
            len(chromo_list),
            mol_id_dict=mol_id_dict,
            **carrier_kwargs,
        )
        continue_sim = True
        while continue_sim:
            continue_sim = i_carrier.calculate_hop(chromo_list, verbose=verbose)
        # Now the carrier has finished hopping, let's calculate its vitals
        i_carrier.update_displacement()

        t2 = time.perf_counter()
        elapsed_time = float(t2) - float(t1)
        time_str = hf.time_units(elapsed_time)

        v_print(
            f"\t{i_carrier.c_type} hopped {i_carrier.n_hops} times over "
            + f"{i_carrier.current_time:.2e} seconds "
            + f"into image {i_carrier.image} for a displacement of"
            + f"\n\t{i_carrier.displacement:.2f} (took walltime {time_str})",
            verbose,
            filename=filename,
        )
        carrier_list.append(i_carrier)
    t3 = time.perf_counter()
    elapsed_time = float(t3) - float(t0)
    time_str = hf.time_units(elapsed_time)
    if send_end is not None:
        send_end.send(carrier_list)
    else:
        return carrier_list


def snap_molecule_indices(snap):
    """Find molecule index for each particle.

    Given a snapshot from a trajectory, compute clusters of bonded molecules
    and return an array of the molecule index of each particle.

    Parameters
    ----------
    snap : gsd.hoomd.Snapshot
        Trajectory snapshot.

    Returns
    -------
    numpy array (N_particles,)
    """
    system = freud.AABBQuery.from_system(snap)
    n_query_pts = n_pts = snap.bonds.N
    query_pt_inds = snap.bonds.group[:, 0]
    pt_inds = snap.bonds.group[:, 1]
    distances = system.box.compute_distances(
        system.points[query_pt_inds], system.points[pt_inds]
    )
    nlist = freud.NeighborList.from_arrays(
        n_query_pts, n_pts, query_pt_inds, pt_inds, distances
    )
    cluster = freud.cluster.Cluster()
    cluster.compute(system=system, neighbors=nlist)
    return cluster.cluster_idx


def get_molecule_ids(snap, chromo_list):
    print("Determining molecule IDs...")
    # Determine molecules based on bonds
    molecules = snap_molecule_indices(snap)
    # Convert to chromo_mol_id
    chromo_mol_id = {}
    for i, chromo in enumerate(chromo_list):
        chromo_mol_id[i] = molecules[chromo.atom_ids[0]]
    return chromo_mol_id


def get_jobslist(sim_times, n_holes=0, n_elec=0, nprocs=None, seed=None):
    # Get the random seed now for all the child processes
    if seed is not None:
        np.random.seed(seed)
    if nprocs is None:
        nprocs = mp.cpu_count()
    # Determine the maximum simulation times based on the parameter dictionary
    carriers = []
    # Modification: Rather than being clever here with the carriers, I'm just
    # going to create the master list of jobs that need running and then
    # randomly shuffle it. This will hopefully permit a similar number of holes
    # and electrons and lifetimes to be run simultaneously providing adequate
    # statistics more quickly
    for lifetime in sim_times:
        for carrier_no in range(n_holes):
            carriers.append([carrier_no, lifetime, "hole"])
        for carrier_no in range(n_elec):
            carriers.append([carrier_no, lifetime, "electron"])
    np.random.shuffle(carriers)
    step = math.ceil(len(carriers) / nprocs)
    jobs_list = [carriers[i : i + step] for i in range(0, len(carriers), step)]
    return jobs_list


def run_kmc(
    jobs_list,
    KMC_directory,
    chromo_list,
    snap,
    temp,
    combine_KMC_results=True,
    carrier_kwargs={},
    verbose=0,
    ): # pragma: no cover
    running_jobs = []
    pipes = []

    for cpu_rank, jobs in enumerate(jobs_list):
        child_seed = np.random.randint(0, 2 ** 32)

        recv_end, send_end = mp.Pipe(False)
        p = mp.Process(
            target=run_single_kmc,
            args=(jobs, KMC_directory, chromo_list, snap, temp),
            kwargs={
                "carrier_kwargs": carrier_kwargs,
                "seed": child_seed,
                "send_end": send_end,
                "verbose": verbose,
                "cpu_rank": cpu_rank,
            },
        )
        running_jobs.append(p)
        pipes.append(recv_end)
        p.start()

    # wait for all jobs to finish
    for p in running_jobs:
        p.join()

    carriers_lists = [x.recv() for x in pipes]

    carriers = [item for sublist in carriers_lists for item in sublist]
    # Now combine the carrier data
    v_print("All KMC jobs completed!", verbose)
    if combine_KMC_results:
        v_print("Combining outputs...", verbose)

        combined_data = {}
        for carrier in carriers:
            d = carrier.__dict__
            for key, val in d.items():
                if key in ["initial_chromo", "current_chromo"]:
                    val = val.center
                    key = key.split("_")[0] + "_position"
                if key not in ["hole_history", "electron_history"]:
                    val = [val]
                if key not in combined_data:
                    combined_data[key] = val
                else:
                    try:
                        combined_data[key] += val
                    except TypeError:
                        # catch errors trying to add None and None
                        combined_data[key] = val
        return combined_data
    return carriers
