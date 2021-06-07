import itertools
import multiprocessing as mp

import numpy as np


# UNIVERSAL CONSTANTS, DO NOT CHANGE!
elem_chrg = 1.60217657e-19  # C
k_B = 1.3806488e-23  # m^{2} kg s^{-2} K^{-1}
hbar = 1.05457173e-34  # m^{2} kg s^{-1}


def box_points(box):
    """Arrange the corner coordinates such that following them draws a cube.

    Parameters
    ----------
    box : numpy.ndarray
        Box array containing [Lx, Ly, Lz]

    Returns
    -------
    numpy.ndarray (16,3)
        xyz coordinates of the box corners in order such that if a line is drawn
        between each point in order it draws a cube.
    """
    dims = np.array([(-i / 2, i / 2) for i in box / 10])
    corners = [
        (dims[0, i], dims[1, j], dims[2, k])
        for i, j, k in itertools.product(range(2), repeat=3)
    ]
    corner_ids = [0, 1, 3, 1, 5, 4, 0, 2, 3, 7, 5, 4, 6, 2, 6, 7]
    box_pts = np.array([corners[i] for i in corner_ids])
    return box_pts


def v_print(string, verbosity, v_level=0, filename=None):  # pragma: no cover
    """Print based on verbosity level.

    If verbosity is greater than v_level, v_print will print.

    Parameters
    ----------
    string : str
        The string to print.
    verbosity : int
        The current verbosity level
    v_level : int, default 0
        The level above which string should print
    filename : path, default None
        Path to a file where the vprint will be written instad of printed to
        stdout. If None is given, v_print prints to stdout.
    """
    if verbosity > v_level:
        if filename is None:
            print(string)
        else:
            with open(filename, "a") as f:
                f.write(f"{string}\n")


def time_units(elapsed_time, precision=2):
    """Convert elapsed time in seconds to its largest unit.

    Example
    -------
    >>> time_units(12345)
    '3.43 hours'

    Parameters
    ----------
    elapsed_time : float
        Elapsed time in seconds
    precision : int, default 2
        Number of values after decimal place to display.

    Returns
    -------
    str
        Elapsed time formatted as a string
    """
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
    """Sort a pair of lists by the first list in ascending order.

    For example, given lists of mass and position, it will sort by mass and
    return lists such that mass[i] still corresponds to position[i]

    Parameters
    ----------
    list1 : list
        The list to be sorted by
    list2 : list
        Another list which will be rearranged to maintain order with list1

    Returns
    -------
    list of lists
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
    """Get the hopping rate.

    Parameters
    ----------
    lambd : float
        The reorganization energy in eV.
    ti : float
        The transfer integral between the chromophores in eV.
    delta_e : float
        The energy difference between the frontier orbitals of the chromophores
        in eV.
    prefactor : float
        A prefactor to the rate equation.
    temp : float
        The temperature in Kelvin.
    use_vrh : bool, default False
        Whether to use variable-range hopping.
    rij : float, default 0.0
        The distance between the chromophores in Angstroms. (only used with VRH)
    vrh : float, default 1.0
        A cutoff distance in Angstroms. (only used with VRH)
    boltz : bool, default False
        Whether to apply a simple Boltzmann energy penalty.

    Returns
    -------
    float
        The hopping rate in inverse seconds.
    """
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
    """Get the time an event would take.

    Parameters
    ----------
    rate : float
        The rate in inverse seconds
    slowest : float, default None
        The slowest allowed time. (Only used if max_attempts is not None.)
    fastest : float, default None
        The fastest allowed time. (Only used if max_attempts is not None.)
    max_attempts : int, default None
        Number of attempts allowed to obtain a time between slowest and fastest.
    log_file : path, default None
        Path to a file to which to log output. If None is given, output is
        printed to stdout.
    verbose : int, default 0
        The verbosity level.

    Returns
    -------
    float
        The time in seconds the event would take given its rate.
    """
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


def find_axis(atom1, atom2, normalize=True):
    """Find normalized vector from atom1 to atom2."""
    sep = atom2 - atom1
    if normalize is True:
        norm = np.linalg.norm(sep)
        if norm == 0:
            return sep
        return sep / norm
    return sep
