import itertools
import multiprocessing as mp

import numpy as np


# UNIVERSAL CONSTANTS, DO NOT CHANGE!
elem_chrg = 1.60217657e-19  # C
k_B = 1.3806488e-23  # m^{2} kg s^{-2} K^{-1}
hbar = 1.05457173e-34  # m^{2} kg s^{-1}


def box_points(box):
    dims = np.array([(-i / 2, i / 2) for i in box / 10])
    corners = [
        (dims[0, i], dims[1, j], dims[2, k])
        for i, j, k in itertools.product(range(2), repeat=3)
    ]
    corner_ids = [0, 1, 3, 1, 5, 4, 0, 2, 3, 7, 5, 4, 6, 2, 6, 7]
    box_pts = np.array([corners[i] for i in corner_ids])
    return box_pts


def v_print(string, verbosity, v_level=0, filename=None):  # pragma: no cover
    if verbosity > v_level:
        if filename is None:
            print(string)
        else:
            with open(filename, "a") as f:
                f.write(f"{string}\n")


def time_units(elapsed_time, precision=2):
    if elapsed_time < 60:
        time_units = "seconds"
    elif elapsed_time < 3600:
        elapsed_time /= 60.0
        time_units = "minutes"
    elif elapsed_time < 86400:
        elapsed_time /= 3600.0
        time_units = "hours"
    else:
        elapsed_time /= 86400.0
        time_units = "days"
    return f"{elapsed_time:.{precision}f} {time_units}"


def parallel_sort(list1, list2):
    """
    Sort a pair of lists by the first list in ascending order.
    (e.g., given mass and position, it will sort by mass and return lists
    such that mass[i] still corresponds to position[i])
    """
    types = [None, None]

    for i, l in enumerate([list1, list2]):
        if isinstance(l, np.ndarray):
            types[i] = "array"
        elif isinstance(l, list):
            types[i] = "list"

    list1, list2 = zip(*sorted(zip(list1, list2)))
    lists = [list1, list2]

    for i, t in enumerate(types):
        if t == "array":
            lists[i] = np.array(lists[i])
        elif t == "list":
            lists[i] = list(lists[i])

    return lists


def get_hop_rate(
    lambd,
    ti,
    delta_e,
    prefactor,
    temp,
    use_vrh=False,
    rij=0.0,
    vrh=1.0,
    boltz=False,
):
    # Based on the input parameters, can make this the semiclassical Marcus
    # Hopping Rate Equation, or a more generic Miller Abrahams-based hop
    # Firstly, to prevent divide-by-zero errors:
    if ti == 0.0:
        return 0
    # Regardless of hopping type, sort out the prefactor first:
    lambd *= elem_chrg
    ti *= elem_chrg
    delta_e *= elem_chrg

    k = prefactor * (2 * np.pi / hbar) * (ti ** 2)
    k *= np.sqrt(1 / (4 * lambd * np.pi * k_B * temp))

    # VRH?
    if use_vrh is True:
        k *= np.exp(-rij / vrh)
    # Simple Boltzmann energy penalty?
    if boltz is True:
        # Only apply the penalty if delta_e is positive
        if delta_e > 0.0:
            k *= np.exp(-delta_e / (k_B * temp))
        # Otherwise, k *= 1
    else:
        k *= np.exp(-((delta_e + lambd) ** 2) / (4 * lambd * k_B * temp))
    return k


def get_event_tau(
    rate,
    slowest=None,
    fastest=None,
    max_attempts=None,
    log_file=None,
    verbose=0,
):
    if rate == 0:
        # If rate == 0, then make the hopping time extremely long
        return 1e99

    # Use the KMC algorithm to determine the wait time to this hop
    counter = 0
    if max_attempts is None:
        x = np.random.random()
        while x == 0 or x == 1:
            x = np.random.random()
        tau = -np.log(x) / rate
        return tau
    while counter < max_attempts:
        x = np.random.random()
        # Ensure that we don't get exactly 0.0 or 1.0, which would break our
        # logarithm
        if (x == 0.0) or (x == 1.0):
            continue
        tau = -np.log(x) / rate
        if (fastest is not None) and (slowest is not None):
            if fastest < tau < slowest:
                return tau
            counter += 1

    err_msg = f"""
    Attempted {max_attempts:d} times to obtain an event timescale within the
    tolerances: {fastest:.2e} <= tau < {slowest:.2e} with the given rate
    {rate:.2e}, without success.
    Permitting the event anyway with tau={tau:.2e}...
    """
    v_print(err_msg, verbose, filename=log_file)
    return tau


# TODO delete?
# functions below here do not seem to be used


def find_axis(atom1, atom2, normalize=True):
    """
    This function determines the normalized vector from the location of
    atom1 to atom2.
    """
    sep = atom2 - atom1
    if normalize is True:
        norm = np.linalg.norm(sep)
        if norm == 0:
            return sep
        return sep / norm
    return sep


def get_coords_normvec(coords):  # pragma: no cover
    AB = coords[1] - coords[0]
    AC = coords[2] - coords[0]
    return np.cross(AB, AC) / np.linalg.norm(normal)


def get_rotation_matrix(vector1, vector2):  # pragma: no cover
    """
    This function returns the rotation matrix around the origin that maps
    vector1 to vector 2
    """
    cross_product = np.cross(vector1, vector2)
    sin_angle = np.sqrt(np.sum(cross_product ** 2))
    cos_angle = np.dot(vector1, vector2)
    skew_matrix = np.array(
        [
            [0, -cross_product[2], cross_product[1]],
            [cross_product[2], 0, -cross_product[0]],
            [-cross_product[1], cross_product[0], 0],
        ]
    )
    skew_matrix_squared = skew_matrix @ skew_matrix
    rot_matrix = (
        np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        + skew_matrix
        + skew_matrix_squared * ((1 - cos_angle) / (sin_angle ** 2))
    )
    return rot_matrix


def get_nprocs():  # pragma: no cover
    # Determine the number of available processors, either by querying the
    # SLURM_NPROCS environment variable, or by using multiprocessing to count
    # the number of visible CPUs.
    try:
        procs = int(os.environ.get("SLURM_NPROCS"))
    except TypeError:
        # Was not loaded using SLURM, so use all physical processors
        procs = mp.cpu_count()
    return procs


def get_FRET_hop_rate(
    prefactor, lifetime, r_F, rij, delta_e, temp
):  # pragma: no cover
    # Foerster Transport Hopping Rate Equation
    # The prefactor included here is a bit of a bodge to try and get the
    # mean-free paths of the excitons more in line with the 5nm of experiment.
    # Possible citation: 10.3390/ijms131217019 (they do not do the simulation
    # they just point out some limitations of FRET which assumes point-dipoles
    # which does not necessarily work in all cases)
    if delta_e <= 0:
        boltz = 1
    else:
        boltz = np.exp(-elem_chrg * delta_e / (k_B * temp))
    k_FRET = (prefactor / lifetime) * (r_F / rij) ** 6 * boltz
    return k_FRET


def get_miller_abrahams_hop_rate(
    prefactor, separation, radius, delta_e, temp
):  # pragma: no cover
    k = prefactor * np.exp(-2 * separation / radius)
    if delta_e > 0:
        k *= np.exp(-delta_e / (k_B * temp))
    return k
