import math
import multiprocessing as mp
import os
import pickle
import signal
import time
import traceback

import numpy as np
from scipy.sparse import lil_matrix

from morphct import helper_functions as hf


elem_chrg = 1.60217657e-19  # C
k_B = 1.3806488e-23  # m^{2} kg s^{-2} K^{-1}
hbar = 1.05457173e-34  # m^{2} kg s^{-1}
log_file = None


class Carrier:
    def __init__(
        self,
        chromo,
        lifetime,
        carrier_no,
        box,
        temp,
        mol_ID_dict,
        hop_limit = None,
        record_history=True,
        use_average_hop_rates = False,
        average_intra_hop_rate = None,
        average_inter_hop_rate = None,
        use_koopmans = False,
        boltz_penalty = False,
        use_VRH = False,
        hopping_prefactor = 1.0
    ):
        self.id = carrier_no
        self.image = np.array([0, 0, 0])
        self.initial_chromo = chromo
        self.current_chromo = chromo
        self.lambda_ij = self.current_chromo.reorganisation_energy
        self.hop_limit = hop_limit
        self.temp = temp
        self.lifetime = lifetime
        self.current_time = 0.0

        self.hole_history = None
        self.electron_history = None
        n = len(chromo_list)
        if self.current_chromo.species == "donor":
            self.carrier_type = "hole"
            if record_history:
                self.hole_history = lil_matrix((n, n), dtype=int)
        elif self.current_chromo.species == "acceptor":
            self.carrier_type = "electron"
            if record_history:
                self.electron_history = lil_matrix((n, n), dtype=int)

        self.n_hops = 0
        self.box = box
        self.displacement = None
        self.mol_ID_dict = mol_ID_dict

        self.use_average_hop_rates = use_average_hop_rates
        self.average_intra_hop_rate = average_intra_hop_rate
        self.average_inter_hop_rate = average_inter_hop_rate

        # Set the use of Koopmans' approximation to false if the key does not
        # exist in the parameter dict
        self.use_koopmans = use_koopmans
        # Are we using a simple Boltzmann penalty?
        self.boltz_penalty = boltz_penalty
        # Are we applying a distance penalty beyond the transfer integral?
        self.use_VRH = use_VRH
        if self.use_VRH:
            self.VRH_delocalisation = (self.current_chromo.VRH_delocalisation)

        self.hopping_prefactor = hopping_prefactor

    def calculate_hop(self, chromo_list):
        # Terminate if the next hop would be more than the termination limit
        if self.hop_limit is not None:
            if self.no_hops + 1 > self.hop_limit:
                return 1
        # Determine the hop times to all possible neighbors
        hop_times = []
        if self.use_average_hop_rates is True:
            # Use the average hop values given in the parameter dict to pick a
            # hop
            for i,img in self.current_chromo.neighbors:
                neighbor = chromo_list[i]
                assert neighbor.id == i
                if (
                    self.mol_ID_dict[self.current_chromo.id]
                    == self.mol_ID_dict[neighbor.id]
                    ):
                    hop_rate = self.average_intra_hop_rate
                else:
                    hop_rate = self.average_inter_hop_rate
                hop_time = hf.determine_event_tau(hop_rate)
                # Keep track of the chromophoreID and the corresponding tau
                hop_times.append([neighbor.ID, hop_time])
        else:
            # Obtain the reorganisation energy in J (from eV in the parameter
            # file)
            for i_neighbor, ti in enumerate(self.current_chromo.neighbors_ti):
                # Ignore any hops with a NoneType transfer integral (usually
                # due to an orca error)
                if transfer_integral is None:
                    continue
                n_ind = self.current_chromo.neighbors[i_neighbor][0]
                delta_E_ij = self.current_chromo.neighbors_delta_e[i_neighbor]
                # Load the specified hopping prefactor
                prefactor = self.hopping_prefactor
                # Get the relative image so we can update the carrier image
                # after the hop
                rel_image = self.current_chromo.neighbors[i_neighbor][1]
                # All of the energies are in eV currently, so convert them to J
                if self.use_VRH is True:
                    neighbor_chromo = chromo_list[n_ind]
                    neighbor_pos = neighbor_chromo.pos + rel_image * self.box

                    # Chromophore separation needs converting to m
                    chromo_separation = (
                        hf.calculate_separation(
                            self.current_chromo.pos, neighbor_pos
                        )
                        * 1e-10
                    )
                    hop_rate = hf.calculate_carrier_hop_rate(
                        self.lambda_ij * elem_chrg,
                        ti * elem_chrg,
                        delta_E_ij * elem_chrg,
                        prefactor,
                        self.temp,
                        use_VRH=True,
                        rij=chromo_separation,
                        VRH_delocalisation=self.VRH_delocalisation,
                        boltz_pen=self.boltz_penalty,
                    )
                else:
                    hop_rate = hf.calculate_carrier_hop_rate(
                        self.lambda_ij * elem_chrg,
                        ti * elem_chrg,
                        delta_E_ij * elem_chrg,
                        prefactor,
                        self.temp,
                        boltz_pen=self.boltz_penalty,
                    )
                hop_time = hf.determine_event_tau(hop_rate)
                # Keep track of the chromophoreID and the corresponding tau
                hop_times.append([n_ind, hop_time, rel_image])
        # Sort by ascending hop time
        hop_times.sort(key=lambda x: x[1])
        if len(hop_times) == 0:
            # We are trapped here, so create a dummy hop with time 1E99
            hop_times = [[self.current_chromo.ID, 1e99, [0, 0, 0]]]
        # As long as we're not limiting by the number of hops:
        if self.hop_limit is None:
            # Ensure that the next hop does not put the carrier over its
            # lifetime
            if (self.current_time + hop_times[0][1]) > self.lifetime:
                # Send the termination signal to singleCoreRunKMC.py
                return 1
        # Move the carrier and send the contiuation signal to
        # singleCoreRunKMC.py
        # Take the quickest hop
        self.perform_hop(
            chromo_list[hop_times[0][0]], hop_times[0][1], hop_times[0][2]
        )
        return 0

    def perform_hop(self, destination_chromo, hop_time, rel_image):
        initial_ID = self.current_chromo.ID
        destination_ID = destination_chromo.ID
        self.image += rel_image
        # Carrier image now sorted, so update its current position
        self.current_chromo = destination_chromo
        # Increment the simulation time
        self.current_time += hop_time
        # Increment the hop counter
        self.no_hops += 1
        # Now update the sparse history matrix
        if (self.carrier_type.lower() == "hole") and (
            self.hole_history_matrix is not None
        ):
            self.hole_history_matrix[initial_ID, destination_ID] += 1
        elif (self.carrier_type.lower() == "electron") and (
            self.electron_history_matrix is not None
        ):
            self.electron_history_matrix[initial_ID, destination_ID] += 1


class termination_signal:
    kill_sent = False

    def __init__(self):
        signal.signal(signal.SIGINT, self.catch_kill)
        signal.signal(signal.SIGTERM, self.catch_kill)

    def catch_kill(self, signum, frame):
        self.kill_sent = True


class terminate(Exception):
    """This class is raised to terminate a KMC simulation"""

    def __init__(self, string):
        self.string = string

    def __str__(self):
        return self.string


def save_pickle(save_data, save_pickle_name):
    with open(save_pickle_name, "wb+") as pickle_file:
        pickle.dump(save_data, pickle_file)
    hf.write_to_file(
        log_file, [f"Pickle file saved successfully as {save_pickle_name}"]
    )


def calculate_displacement(init_pos, final_pos, final_image, sim_dims):
    displacement = (
            final_pos - init_pos + final_image *
            (sim_dims[:,1] - sim_dims[:,0])
            )
    return np.linalg.norm(displacement)


def initialise_save_data(n_chromos, seed):
    return {
        "seed": seed,
        "ID": [],
        "image": [],
        "lifetime": [],
        "current_time": [],
        "no_hops": [],
        "displacement": [],
        "hole_history_matrix": lil_matrix((n_chromos, n_chromos), dtype=int),
        "electron_history_matrix": lil_matrix(
            (n_chromos, n_chromos), dtype=int
        ),
        "initial_position": [],
        "final_position": [],
        "carrier_type": [],
    }


def split_molecules(input_dictionary, chromo_list):
    # Split the full morphology into individual molecules
    # Create a lookup table `neighbor list' for all connected atoms called
    # {bondedAtoms}
    bonded_atoms = hf.obtain_bonded_list(input_dictionary["bond"])
    molecule_list = [i for i in range(len(input_dictionary["type"]))]
    # Recursively add all atoms in the neighbor list to this molecule
    for mol_ID in range(len(molecule_list)):
        molecule_list = update_molecule(mol_ID, molecule_list, bonded_atoms)
    # Here we have a list of len(atoms) where each index gives the molID
    mol_ID_dict = {}
    for chromo in chromo_list:
        AAID_to_check = chromo.AAIDs[0]
        mol_ID_dict[chromo.ID] = molecule_list[AAID_to_check]
    return mol_ID_dict


def update_molecule(atom_ID, molecule_list, bonded_atoms):
    # Recursively add all neighbors of atom number atomID to this molecule
    try:
        for bonded_atom in bonded_atoms[atom_ID]:
            # If the moleculeID of the bonded atom is larger than that of the
            # current one, update the bonded atom's ID to the current one's to
            # put it in this molecule, then iterate through all of the bonded
            # atom's neighbors
            if molecule_list[bonded_atom] > molecule_list[atom_ID]:
                molecule_list[bonded_atom] = molecule_list[atom_ID]
                molecule_list = update_molecule(
                    bonded_atom, molecule_list, bonded_atoms
                )
            # If the moleculeID of the current atom is larger than that of the
            # bonded one, update the current atom's ID to the bonded one's to
            # put it in this molecule, then iterate through all of the current
            # atom's neighbors
            elif molecule_list[bonded_atom] < molecule_list[atom_ID]:
                molecule_list[atom_ID] = molecule_list[bonded_atom]
                molecule_list = update_molecule(
                    atom_ID, molecule_list, bonded_atoms
                )
            # Else: both the current and the bonded atom are already known to
            # be in this molecule, so we don't have to do anything else.
    except KeyError:
        # This means that there are no bonded CG sites (i.e. it's a single molecule)
        pass
    return molecule_list


def run_single_kmc(
        jobs_to_run,
        KMC_directory,
        AA_morphdict,
        param_dict,
        chromo_list,
        cpu_rank,
        seed,
        send_end,
        use_avg_hop_rates=False,
        record_history=True
        ):
    global log_file

    log_file = os.path.join(KMC_directory, f"KMC_log_{cpu_rank:02d}.log")
    pickle_file = os.path.join(KMC_directory, f"kmc_{cpu_rank:02d}.pickle")
    # Reset the log file
    with open(log_file, "wb+") as log_file_handle:
        pass
    hf.write_to_file(log_file, [f"Found {len(jobs_to_run):d} jobs to run"])
    # Set the affinities for this current process to make sure it's maximising
    # available CPU usage
    current_PID = os.getpid()
    if use_avg_hop_rates is True:
        # Chosen to split hopping by inter-intra molecular hops, so get
        # molecule data
        mol_ID_dict = split_molecules(AA_morphdict, chromo_list)
        # molIDDict is a dictionary where the keys are the chromoIDs, and
        # the vals are the molIDs
    else:
        mol_ID_dict = None
    # Attempt to catch a kill signal to ensure that we save the pickle before
    # termination
    killer = termination_signal()
    # Save the pickle as a list of `saveCarrier' instances that contain the
    # bare minimum
    save_data = initialise_save_data(len(chromo_list), seed)
    if record_history is False:
        save_data["hole_history_matrix"] = None
        save_data["electron_history_matrix"] = None
    t0 = time.perf_counter()
    save_time = time.perf_counter()
    save_slot = "slot1"
    try:
        for job_number, [carrier_no, lifetime, ctype] in enumerate(jobs_to_run):
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
                chromo_list,
                param_dict,
                i_chromo,
                lifetime,
                carrier_no,
                AA_morphdict,
                mol_ID_dict,
            )
            terminate_simulation = False
            while terminate_simulation is False:
                terminate_simulation = bool(
                        i_carrier.calculate_hop(chromo_list)
                )
                if killer.kill_sent is True:
                    raise terminate(
                        "Kill command sent, terminating KMC simulation..."
                    )
            # Now the carrier has finished hopping, let's calculate its vitals
            init_position = this_carrier.initial_chromo.pos
            final_position = this_carrier.current_chromo.pos
            final_image = this_carrier.image
            sim_dims = this_carrier.sim_dims
            this_carrier.displacement = calculate_displacement(
                initial_position, final_position, final_image, sim_dims
            )
            # Now the calculations are completed, create a barebones class
            # containing the save data
            importantData = [
                "ID",
                "image",
                "lifetime",
                "current_time",
                "no_hops",
                "displacement",
                "carrier_type",
            ]
            for name in importantData:
                save_data[name].append(i_carrier.__dict__[name])
            # Update the carrierHistoryMatrix
            if record_history is True:
                if i_carrier.carrier_type == "hole":
                    save_data["hole_history_matrix"] += (
                            i_carrier.hole_history_matrix
                            )
                elif i_carrier.carrier_type == "electron":
                    save_data["electron_history_matrix"] += (
                            i_carrier.electron_history_matrix
                            )
            # Then add in the initial and final positions
            save_data["initial_position"].append(init_position)
            save_data["final_position"].append(final_position)
            t2 = time.perf_counter()
            elapsed_time = float(t2) - float(t1)
            if elapsed_time < 60:
                time_units = "seconds."
            elif elapsed_time < 3600:
                elapsed_time /= 60.0
                time_units = "minutes."
            elif elapsed_time < 86400:
                elapsed_time /= 3600.0
                time_units = "hours."
            else:
                elapsed_time /= 86400.0
                time_units = "days."

            write_lines = [
                    "{0:s} hopped {1:d} times, over {2:.2e} seconds, ".format(
                        i_carrier.carrier_type.capitalize(),
                        i_carrier.no_hops,
                        i_carrier.current_time,
                        ),
                    f"into image {repr(this_carrier.image)}",
                    ", for a displacement of {0:.2f}, in {1:.2f} ".format(
                        i_carrier.displacement,
                        elapsed_time
                        ),
                    f"wall-clock {time_units:s}"
                    ]
            hf.write_to_file(log_file, ["".join(write_lines)])
            # Save the pickle file every hour
            if (t2 - save_time) > 3600:
                print(
                    "Completed {0:d}/{1:d} jobs. Checkpoint at {2:3d}%".format(
                        job_number,
                        len(jobs_to_run),
                        np.round((job_number + 1) / len(jobs_to_run) * 100)
                    )
                )
                write_str= "Completed {0:d}/{1:d} jobs {2:3d}%".format(
                        job_number,
                        len(jobs_to_run),
                        round((job_number + 1) / len(jobs_to_run) * 100)
                        )
                hf.write_to_file(log_file, [write_str])
                save_pickle(save_data, pickle_file)
                if save_slot.lower() == "slot1":
                    save_slot = "slot2"
                elif save_slot.lower() == "slot2":
                    save_slot = "slot1"
                save_time = time.perf_counter()
    except Exception as error_message:
        print(traceback.format_exc())
        print("Saving the pickle file cleanly before termination...")
        hf.write_to_file(log_file, [str(error_message)])
        hf.write_to_file(
            log_file, ["Saving the pickle file cleanly before termination..."]
        )
        save_pickle(save_data, pickle_file)
        print("Pickle saved! Exiting Python...")
        exit()
    t3 = time.perf_counter()
    elapsed_time = float(t3) - float(t0)
    if elapsed_time < 60:
        time_units = "seconds."
    elif elapsed_time < 3600:
        elapsed_time /= 60.0
        time_units = "minutes."
    elif elapsed_time < 86400:
        elapsed_time /= 3600.0
        time_units = "hours."
    else:
        elapsed_time /= 86400.0
        time_units = "days."
    hf.write_to_file(
        log_file,
        ["All jobs completed in {0:.2f}{1:s}".format(elapsed_time, time_units)],
    )
    hf.write_to_file(
            log_file,
            ["Saving the pickle file cleanly before termination..."]
            )
    save_pickle(save_data, pickle_file)
    hf.write_to_file(log_file, ["Exiting normally..."])
    send_end.send(save_data)


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
    jobs_list = [carriers[i : i+step] for i in range(0, len(carriers), step)]
    return jobs_list


def run_kmc(
        jobs_list,
        KMC_directory,
        AA_morphdict,
        param_dict,
        chromo_list
        ):
    running_jobs = []
    pipes = []

    for proc_ID, jobs in enumerate(jobs_list):
        child_seed = np.random.randint(0, 2 ** 32)

        recv_end, send_end = mp.Pipe(False)
        p = mp.Process(target=run_single_kmc, args=(
            jobs,
            KMC_directory,
            AA_morphdict,
            param_dict,
            chromo_list,
            proc_ID,
            child_seed,
            send_end
        ))
        running_jobs.append(p)
        pipes.append(recv_end)
        p.start()

    # wait for all jobs to finish
    for p in running_jobs:
        p.join()

    carrier_data_list = [x.recv() for x in pipes]

    # Now combine the carrier data
    print("All KMC jobs completed!")
    if param_dict["combine_KMC_results"] is True:
        print("Combining outputs...")
        combined_data = {}
        for carrier_data in carrier_data_list:
            for key, val in carrier_data.items():
                    if key not in combined_data:
                        combined_data[key] = val
                    else:
                        combined_data[key] += val
        return combined_data
    return carrier_data_list
