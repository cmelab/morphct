import os
import pickle
import sys
import multiprocessing as mp
import numpy as np
import subprocess as sp
from morphct.definitions import PROJECT_ROOT, SINGLE_ORCA_RUN_FILE
from morphct.code import helper_functions as hf


def create_inputs(chromo_list, AA_morphdict, param_dict):
    # Singles first
    for chromophore in chromo_list:
        # Include the molecule terminating units on the required atoms of the
        # chromophore
        if chromophore.terminate is True:
            terminating_group_positions = terminate_monomers(
                chromophore, param_dict, AA_morphdict
            )
            terminating_group_images = [chromophore.image] * len(
                terminating_group_positions
            )
        else:
            terminating_group_positions = None
            terminating_group_images = None
        chromophore.qcc_input = write_qcc_inp(
            AA_morphdict,
            chromophore.AAIDs,
            [chromophore.image] * len(chromophore.AAIDs),
            terminating_group_positions,
            terminating_group_images,
        )
    # Determine how many pairs there are first:
    n_pairs = np.sum([len(chromo.neighbors) for chromo in chromo_list])
    print(f"There are {n_pairs // 2} total neighbor pairs to consider.")
    # /2 because the forwards and backwards hops are identical
    # Then consider each chromophore against every other chromophore
    qcc_pairs = []
    for chromo1 in chromo_list:
        neighbors_ID = [neighbor[0] for neighbor in chromo1.neighbors]
        neighbors_image = [neighbor[1] for neighbor in chromo1.neighbors]
        for chromo2 in chromo_list:
            # Skip if chromo2 is not one of chromo1's neighbors
            # Also skip if chromo2's ID is < chromo1's ID to prevent
            # duplicates
            if (chromo2.ID not in neighbors_ID) or (chromo2.ID < chromo1.ID):
                continue
            # Update the qcc input name
            pair = (chromo1.ID, chromo2.ID)
            # Find the correct relative image for the neighbor chromophore
            chromo2_rel_image = neighbors_image[neighbors_ID.index(chromo2.ID)]
            chromo2_transformation = list(
                np.array(chromo1.image)
                - np.array(chromo2.image)
                + np.array(chromo2_rel_image)
            )
            # Find the dimer AAIDs and relative images for each atom
            AAIDs = chromo1.AAIDs + chromo2.AAIDs
            images = [[0, 0, 0] for i in range(len(chromo1.AAIDs))] + [
                chromo2_transformation for i in range(len(chromo2.AAIDs))
            ]
            # Now add the terminating groups to both chromophores
            # Note that we would only ever expect both chromophores to require
            # termination or neither
            if chromo1.terminate is True:
                term_group_pos1 = terminate_monomers(
                    chromo1, param_dict, AA_morphdict
                )
                term_group_pos2 = terminate_monomers(
                    chromo2, param_dict, AA_morphdict
                )
                # We don't want to add the terminating hydrogens for adjacent
                # monomers, so remove the ones that are within a particular
                # distance
                term_group_pos1, term_group_pos2 = remove_adjacent_terminators(
                    term_group_pos1, term_group_pos2
                )
                terminating_group_images1 = [
                    [0, 0, 0] for i in range(len(term_group_pos1))
                ]
                terminating_group_images2 = [
                    chromo2_transformation for i in range(len(term_group_pos2))
                ]
                # Write the dimer input file
                qcc_input = write_qcc_inp(
                    AA_morphdict,
                    AAIDs,
                    images,
                    term_group_pos1 + term_group_pos2,
                    terminating_group_images1 + terminating_group_images2,
                )
            else:
                # Write the dimer input file
                qcc_input = write_qcc_inp(
                    AA_morphdict,
                    AAIDs,
                    images,
                    None,
                    None,
                )
            qcc_pairs.append((pair,qcc_input))
    return qcc_pairs


def remove_adjacent_terminators(group1, group2):
    pop_list = [[], []]
    for index1, terminal_hydrogen1 in enumerate(group1):
        for index2, terminal_hydrogen2 in enumerate(group2):
            sep = np.linalg.norm(terminal_hydrogen2 - terminal_hydrogen1)
            if sep < 1.2:
                pop_list[0].append(index1)
                pop_list[1].append(index2)
    for group_no, group in enumerate(pop_list):
        for index in sorted(group, reverse=True):
            try:
                [group1, group2][group_no].pop(index)
            except IndexError:
                raise SystemError(
                    """
                    Tried to pop a termination group that does not exist...
                    are you sure this is an atomistic morphology?
                    """
                )
    return group1, group2


def write_qcc_inp(
    AA_morphdict,
    AAIDs,
    images,
    terminal_pos,
    terminal_images,
):
    qcc_lines = []
    all_atom_types = []
    all_positions = []
    numstr = "0123456789"
    lxyz = [AA_morphdict["lx"], AA_morphdict["ly"], AA_morphdict["lz"]]
    # Format the atom positions ready for qcc
    for index, atom_ID in enumerate(AAIDs):
        # Cut the integer bit off the atomType. To allow atom types like Ca and
        # Br, where the atom is defined by one upper- and one lower-case letter,
        # use iter to return the first element of the list of upper case letters'
        # and the first element of the list of lower case letters in the atom
        # type (if any) and .join them together.
        atom_type = AA_morphdict["type"][atom_ID].strip(numstr).capitalize()
        all_atom_types.append(atom_type)
        # Add in the correct periodic images to the position
        all_positions.append(
            AA_morphdict["unwrapped_position"][atom_ID]
            + np.array([(images[index][i] * lxyz[i]) for i in range(3)])
        )
    # Now add in the terminating Hydrogens if necessary
    if terminal_pos is not None:
        for ind, pos in enumerate(terminal_pos):
            # Cut the integer bit off the atomType
            all_atom_types.append("H")
            # Add in the correct periodic images to the position
            all_positions.append(
                    pos + np.array(
                        [terminal_images[ind][i] * lxyz[i] for i in range(3)]
                        )
                    )
    # Now geometrically centralize all of the atoms that are to be included in
    # this input file to make it easier on qcc
    central_position = np.array(
        [
            np.average(np.array(all_positions)[:, 0]),
            np.average(np.array(all_positions)[:, 1]),
            np.average(np.array(all_positions)[:, 2]),
        ]
    )
    # Create the lines to be written in the input file
    for index, position in enumerate(all_positions):
        qcc_lines.append(
            "{0:s}  {1:.5f}  {2:.5f}  {3:.5f};".format(
                all_atom_types[index],
                position[0] - central_position[0],
                position[1] - central_position[1],
                position[2] - central_position[2],
            )
        )
    qcc_input = " ".join(qcc_lines)
    return qcc_input

def terminate_monomers(chromophore, param_dict, AA_morphdict):
    # No CG morphology, so we will use the UA -> AA code definition of which
    # atoms need to have hydrogens added to them.
    new_hydrogen_positions = []
    for atom_index_chromo, atom_index_morph in enumerate(chromophore.AAIDs):
        atom_type = AA_morphdict["type"][atom_index_morph]
        if atom_type not in param_dict["molecule_terminating_connections"].keys():
            continue
        bonded_AAIDs = []
        # Iterate over all termination connections defined for this atomType (in
        # case we are trying to do something mega complicated)
        for connection_info in param_dict["molecule_terminating_connections"][
            atom_type
        ]:
            for [bond_name, AAID1, AAID2] in chromophore.bonds:
                if AAID1 == atom_index_morph:
                    if AAID2 not in bonded_AAIDs:
                        bonded_AAIDs.append(AAID2)
                elif AAID2 == atom_index_morph:
                    if AAID1 not in bonded_AAIDs:
                        bonded_AAIDs.append(AAID1)
            if len(bonded_AAIDs) != connection_info[0]:
                continue
            new_hydrogen_positions += hf.get_terminating_positions(
                AA_morphdict["unwrapped_position"][atom_index_morph],
                [
                    AA_morphdict["unwrapped_position"][bonded_AAID]
                    for bonded_AAID in bonded_AAIDs
                ],
                1,
            )
    # Return terminatingGroups (positions of those hydrogens to be added to the
    # qcc input)
    return new_hydrogen_positions


def get_qcc_jobs(input_dir, param_dict, proc_IDs):
    # First delete any previous log files as we're about to start again with the
    # ZINDO/S calculations
    try:
        os.unlink(input_dir.replace("/input_orca", "/*.log"))
    except OSError:
        pass
    # Obtain a list of files to run
    single_qcc_file_list = os.listdir(os.path.join(input_dir, "single"))
    pair_qcc_file_list = os.listdir(os.path.join(input_dir, "pair"))
    qcc_files_to_run = []
    for file_name in single_qcc_file_list:
        if file_name[-4:] == ".inp":
            qcc_files_to_run.append(os.path.join(input_dir, "single", file_name))
    for file_name in pair_qcc_file_list:
        if file_name[-4:] == ".inp":
            qcc_files_to_run.append(os.path.join(input_dir, "pair", file_name))
    qcc_files_to_run.sort()
    if param_dict["overwrite_current_data"] is False:
        # Do not run any jobs that have already have an output file (and so have
        # at least started to
        # run if not finished)
        pop_list = []
        for job_no, job in enumerate(qcc_files_to_run):
            try:
                with open(
                    job.replace("input_orca", "output_orca").replace(".inp", ".out"),
                    "r",
                ):
                    pop_list.append(job_no)
            except IOError:
                pass
        pop_list.sort(reverse=True)
        for pop_index in pop_list:
            qcc_files_to_run.pop(pop_index)
    if len(qcc_files_to_run) == 0:
        return []
    # Create a jobslist for each procID
    jobs_list = [
        qcc_files_to_run[
            i : i + (int(np.ceil(len(qcc_files_to_run) / len(proc_IDs))))
        ]
        for i in range(
            0,
            len(qcc_files_to_run),
            int(np.ceil(len(qcc_files_to_run) / float(len(proc_IDs)))),
        )
    ]
    return jobs_list


def main(
    AA_morphdict,
    CG_morphdict,
    CGtoAAID_list,
    param_dict,
    chromo_list,
):
    # Get the random seed now for all the child processes
    if param_dict["random_seed_override"] is not None:
        np.random.seed(param_dict["random_seed_override"])
    create_inputs(chromo_list, AA_morphdict, param_dict)
    input_dir = os.path.join(
        param_dict["output_orca_directory"], "chromophores", "input_orca"
    )
    proc_IDs = param_dict["proc_IDs"]
    jobs_list = get_qcc_jobs(input_dir, param_dict, proc_IDs)
    # Shuffle the jobsList to spread it out over the cores
    np.random.shuffle(jobs_list)
    number_of_inputs = sum([len(qcc_files_to_run) for qcc_files_to_run in jobs_list])
    print("Found", number_of_inputs, "qcc files to run.")
    if number_of_inputs > 0:
        # Create pickle file containing the jobs sorted by ProcID to be picked
        # up by single_core_run_qcc.py
        pickle_name = input_dir.replace("input_orca", "orca_jobs.pickle")
        with open(pickle_name, "wb+") as pickle_file:
            pickle.dump(jobs_list, pickle_file)
        print("Orca jobs list written to", pickle_name)
        if len(jobs_list) <= len(proc_IDs):
            proc_IDs = proc_IDs[: len(jobs_list)]
        running_jobs = []
        # Open the required processes to execute the qcc jobs
        for CPU_rank, jobs in enumerate(jobs_list):
            running_jobs.append(
                sp.Popen(
                    [
                        "python",
                        SINGLE_ORCA_RUN_FILE,
                        param_dict["output_orca_directory"],
                        param_dict["output_morphology_directory"],
                        str(CPU_rank),
                        str(int(param_dict["overwrite_current_data"])),
                        str(int(param_dict["remove_orca_inputs"])),
                    ]
                )
            )
        # Wait for all jobs to complete
        [p.wait() for p in running_jobs]
        # Delete the job pickle
        os.system(" ".join(["rm", pickle_name]))
    return (
        AA_morphdict,
        CG_morphdict,
        CGtoAAID_list,
        param_dict,
        chromo_list,
    )


if __name__ == "__main__":
    pass
