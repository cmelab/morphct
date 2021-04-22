from collections import defaultdict
import csv
import itertools
import os

import freud
import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.optimize import curve_fit
from scipy.signal import argrelextrema
from scipy.sparse import lil_matrix
from scipy.stats import linregress

from morphct import helper_functions as hf
from morphct.mobility_kmc import get_molecule_ids


plt = None
p3 = None


def split_carriers(combined_data):
    hole_inds = np.where(np.array(combined_data["c_type"]) == "hole")[0]
    elec_inds = np.where(np.array(combined_data["c_type"]) == "electron")[0]
    hole_data = {}
    elec_data = {}
    for key, val in combined_data.items():
        if key.endswith("history"):
            hole_data[key] = val
            elec_data[key] = val
        else:
            hole_data[key] = [val[i] for i in hole_inds]
            elec_data[key] = [val[i] for i in elec_inds]
    return hole_data, elec_data


def get_times_msds(carrier_data):
    total = 0
    total_averaged = 0
    squared_disps = defaultdict(list)
    actual_times = defaultdict(list)
    for i, displacement in enumerate(carrier_data["displacement"]):
        if (
            carrier_data["current_time"][i] > carrier_data["lifetime"][i] * 2
            or carrier_data["current_time"][i] < carrier_data["lifetime"][i] / 2
            or carrier_data["n_hops"][i] == 1
        ):
            total += 1
            continue
        key = carrier_data["lifetime"][i]
        # A -> m
        squared_disps[key] += [(carrier_data["displacement"][i] * 1e-10) ** 2]
        actual_times[key] += [carrier_data["current_time"][i]]

        # Also keep track of whether each carrier is a hole or an electron
        total_averaged += 1
        total += 1
    if total > total_averaged:
        print(f"\tNotice: The data from {total - total_averaged} carriers were")
        print("\tdiscarded due to the carrier lifetime being more than double")
        print("\t(or less than half of) the specified carrier lifetime.")
    times = []
    msds = []
    time_stderr = []
    msd_stderr = []
    for lt, disps in squared_disps.items():
        times.append(lt)
        time_stderr.append(np.std(actual_times[lt]) / len(actual_times[lt]))
        msds.append(np.average(disps))
        msd_stderr.append(np.std(disps) / len(disps))
    return times, msds, time_stderr, msd_stderr


def plot_displacement_dist(carrier_data, c_type, path):  # pragma: no cover
    plt.figure()
    plt.hist(np.array(carrier_data["displacement"]) * 0.1, bins=60, color="b")
    plt.xlabel(f"{c_type.capitalize()} Displacement (nm)")
    plt.ylabel("Frequency (Arb. U.)")
    filename = f"{c_type}_displacement_dist.png"
    filepath = os.path.join(path, filename)
    plt.savefig(filepath, dpi=300)
    print(f"\tFigure saved as {filename}")
    plt.close()


def plot_cluster_size_dist(clusters, path):  # pragma: no cover
    species = ["donor", "acceptor"]
    for i, cl in enumerate(clusters):
        if cl is not None:
            sizes = [np.log10(len(c)) for c in cl.cluster_keys if len(c) > 5]
        else:
            sizes = None
        if sizes:
            plt.figure()
            plt.hist(
                sizes,
                bins=np.logspace(0, np.ceil(np.max(sizes)), 20),
                color="b",
            )
            plt.xscale("log")
            plt.xlim([1, 10 ** np.ceil(np.max(sizes))])
            plt.xlabel("Cluster Size (Arb. U.)")
            plt.ylabel("Frequency (Arb. U.)")
            filename = f"{species[i]}_cluster_dist.png"
            filepath = os.path.join(path, filename)
            plt.savefig(filepath, dpi=300)
            print(f"\tFigure saved as {filename}")
            plt.close()


def get_connections(chromo_list, carrier_history, box):
    """
    Function to create an array of with a starting point, a vector
    and the number of hops that occured.
    Requires:
        chromo_list,
        carrier_history
    Returns:
        7xN array
    """
    # Create an "empty" array to store data.
    connections = np.zeros(7)
    # Iterate through the chromo_list
    for i, chromo in enumerate(chromo_list):
        # Iterate through the neighbors of the chromophore
        for j, image in chromo.neighbors:
            # Only consider one direction.
            if i < j:
                chromoj = chromo_list[j]
                # Get the vector between the two chromophores.
                if np.array_equal(image, [0, 0, 0]):
                    vector = chromoj.center - chromo.center
                # Account for pbc if not in same relative image.
                else:
                    vector = chromoj.center - chromo.center + image * box

                # Get the net number of times the path was travelled.
                forward = carrier_history[j, i]
                reverse = carrier_history[i, j]
                times = abs(forward - reverse)

                # Append the array if the net times travelled is greater than 0
                if times > 0:
                    datum = np.hstack(
                        (chromo.center, vector, np.array([np.log10(times)]))
                    )
                    connections = np.vstack((connections, datum))
    # Return the array excluding the zeros first line.
    return connections[1:]


def plot_connections(
    chromo_list, carrier_history, c_type, path
):  # pragma: no cover
    # A complicated function that shows connections between carriers in 3D
    # that carriers prefer to hop between.
    # Connections that are frequently used are highlighted in black, whereas
    # rarely used connections are more white.
    # Import matplotlib color modules to set up color bar.
    from matplotlib import colors
    import matplotlib.cm as cmx

    # Create a figure class
    fig = plt.figure(figsize=(7, 6))
    # Make a 3D subplot
    ax = fig.add_subplot(111, projection="3d")

    # Create the array for all the chromophore connections
    connections = get_connections(chromo_list, carrier_history, box)

    # Determine the smallest, non-zero number of times two chromophores
    # are connected.
    vmin = np.min(np.array(connections)[:, 6])
    # Determine the max number of times two chormophores are connected.
    vmax = np.max(np.array(connections)[:, 6])

    # Set up the color bar.
    color_map = plt.get_cmap("inferno")
    c_norm = colors.Normalize(vmin=np.floor(vmin), vmax=np.ceil(vmax))
    scalar_map = cmx.ScalarMappable(norm=c_norm, cmap=color_map)
    hop_colors = scalar_map.to_rgba(connections[:, 6])

    # Set up the intensity for the hops so more travelled paths are more intense
    alphas = connections[:, 6] / vmax
    hop_colors[:, 3] = alphas

    # Plot the vectors between two chromophores
    ax.quiver(
        connections[:, 0],
        connections[:, 1],
        connections[:, 2],
        connections[:, 3],
        connections[:, 4],
        connections[:, 5],
        color=hop_colors,
        arrow_length_ratio=0,
        linewidth=0.7,
    )

    # Draw boxlines
    box_pts = hf.box_points(box)
    ax.plot(box_pts[:, 0], box_pts[:, 1], box_pts[:, 2], c="k", linewidth=1)

    # Make the color bar
    scalar_map.set_array(connections[:, 6])
    tick_location = np.arange(0, vmax, 1)
    # Plot the color bar
    cbar = plt.colorbar(scalar_map, ticks=tick_location, shrink=0.8, aspect=20)
    cbar.ax.set_yticklabels([rf"10$^{{{x}}}$" for x in tick_location])

    # Name and save the figure.
    if c_type == "hole":
        species = "Donor"
    else:
        species = "Acceptor"
    fig.title = f"{species} ({c_type.capitalize()}) Network"

    filename = f"3d_{c_type}_network.png"
    filepath = (os.path.join(path, filename),)
    plt.savefig(filepath, bbox_inches="tight", dpi=300)
    print(f"\tFigure saved as {filename}")

    plt.clf()


def calc_mobility(lin_fit_X, lin_fit_Y, time_err, msd_err, temp):
    # YVals have a std error avmsdError associated with them
    # XVals have a std error avTimeError assosciated with them
    numerator = lin_fit_Y[-1] - lin_fit_Y[0]
    denominator = lin_fit_X[-1] - lin_fit_X[0]

    # Diffusion coeff D = d(MSD)/dt * 1/2n (n = 3 = number of dimensions)
    # Ref: Carbone2014a (Carbone and Troisi)
    diff_coeff = numerator / (6 * denominator)

    # The error in the mobility is proportionally the same as the error in the
    # diffusion coefficient as the other variables are constants without error
    diff_err = diff_coeff * np.sqrt(
        (msd_err / numerator) ** 2 + (time_err / denominator) ** 2
    )

    # Use Einstein-Smoluchowski relation
    # This is in m^{2} / Vs
    mobility = hf.elem_chrg * diff_coeff / (hf.k_B * temp)

    # Convert to cm^{2}/ Vs
    mobility *= 100 ** 2
    mob_err = (diff_err / diff_coeff) * mobility
    return mobility, mob_err


def plot_msd(
    times, msds, time_stderr, msd_stderr, c_type, temp, path
):  # pragma: no cover
    fit_X = np.linspace(np.min(times), np.max(times), 100)
    gradient, intercept, r_val, p_val, std_err = linregress(times, msds)
    print(f"Standard Error {std_err}")
    print(f"Fitting r_val = {r_val}")
    fit_Y = (fit_X * gradient) + intercept
    mobility, mob_error = calc_mobility(
        fit_X, fit_Y, np.average(time_stderr), np.average(msd_stderr), temp
    )
    plt.plot(times, msds)
    plt.errorbar(times, msds, xerr=time_stderr, yerr=msd_stderr)
    plt.plot(fit_X, fit_Y, "r")
    plt.xlabel("Time (s)")
    plt.ylabel(r"MSD (m$^{2}$)")
    plt.title(rf"$\mu_{{0, {c_type[0]}}}$ = {mobility:.3e} cm$^{2}$/Vs", y=1.1)
    filename = f"lin_MSD_{c_type}.png"
    filepath = os.path.join(path, filename)
    plt.savefig(filepath, dpi=300)
    plt.clf()
    print(f"\tFigure saved as {filename}")

    plt.semilogx(times, msds)
    plt.errorbar(times, msds, xerr=time_stderr, yerr=msd_stderr)
    plt.semilogx(fit_X, fit_Y, "r")
    plt.xlabel("Time (s)")
    plt.ylabel(r"MSD (m$^{2}$)")
    plt.title(rf"$\mu_{{0, {c_type[0]}}}$ = {mobility:.3e} cm$^{2}$/Vs", y=1.1)
    filename = f"semi_log_MSD_{c_type}.png"
    filepath = os.path.join(path, filename)
    plt.savefig(filepath, dpi=300)
    plt.clf()
    print(f"\tFigure saved as {filename}")

    plt.plot(times, msds)
    plt.errorbar(times, msds, xerr=time_stderr, yerr=msd_stderr)
    plt.plot(fit_X, fit_Y, "r")
    plt.xlabel("Time (s)")
    plt.ylabel(r"MSD (m$^{2}$)")
    plt.xscale("log")
    plt.yscale("log")
    plt.title(rf"$\mu_{{0,{c_type[0]}}}$ = {mobility:.3e} cm$^{{2}}$/Vs", y=1.1)
    filename = f"log_MSD_{c_type}.png"
    filepath = os.path.join(path, filename)
    plt.savefig(filepath, dpi=300)
    plt.clf()
    print(f"\tFigure saved as {filename}")
    return mobility, mob_error, r_val ** 2


def get_anisotropy(xyzs):
    # First calculate the `center of position' for the particles
    center = np.mean(xyzs, axis=0)
    # First calculate the gyration tensor:
    sxx = sxy = sxz = syy = syz = szz = 0

    shift_xyzs = xyzs - center
    sxx = np.sum(shift_xyzs[:, 0] ** 2)
    sxy = np.sum(shift_xyzs[:, 0] * shift_xyzs[:, 1])
    sxz = np.sum(shift_xyzs[:, 0] * shift_xyzs[:, 2])
    syy = np.sum(shift_xyzs[:, 1] ** 2)
    syz = np.sum(shift_xyzs[:, 1] * shift_xyzs[:, 2])
    szz = np.sum(shift_xyzs[:, 2] ** 2)

    S = np.array([[sxx, sxy, sxz], [sxy, syy, syz], [sxz, syz, szz]])

    eigenvals, eigenvecs = np.linalg.eig(S)
    # We only need the eigenvalues though, no more matrix multiplication
    # Then calculate the relative shape anisotropy (kappa**2)
    return 3 / 2 * np.sum(eigenvals ** 2) / np.sum(eigenvals) ** 2 - 1 / 2


def plot_hop_vectors(
    carrier_data, chromo_list, snap, c_type, path
):  # pragma: no cover
    from matplotlib import colors
    import matplotlib.cm as cmx

    box = snap.configuration.box[:3]
    if c_type == "hole":
        carrier_history = carrier_data["hole_history"]
    else:
        carrier_history = carrier_data["electron_history"]

    nonzero = list(zip(*carrier_history.nonzero()))
    hop_vectors = []
    intensities = []
    for i, j in nonzero:
        # Find the normal of chromo1 plane (first three atoms)
        ichromo_i = chromo_list[i]
        jchromo_j = chromo_list[j]
        normal = hf.get_chromo_normvec(ichromo, snap)

        # Calculate rotation matrix required to map this onto (0, 0, 1)
        rotation_matrix = hf.get_rotation_matrix(normal, np.array([0, 0, 1]))

        # Find the vector from chromoi to chromoj (using correct relative image)
        relative_image = [img for ind, img in ichromo.neighbors if ind == j][0]
        jchromo_center = jchromo.center + relative_image * box

        unrotated_hop_vector = jchromo_center - ichromo_center

        # Apply rotation matrix to chromo1-chromo2 vector
        hop_vector = np.matmul(rotation_matrix, unrotated_hop_vector)

        # Store the hop vector and the intensity
        intensity = np.log10(carrier_history[i, j])
        hop_vectors.append(hop_vector)
        intensities.append(intensity)

    # Convert hop_vectors to an np.array
    hop_vectors = np.array(hop_vectors)
    # Create a figure class
    fig = plt.figure(figsize=(7, 6))
    # Make a 3D subplot
    ax = fig.add_subplot(111, projection="3d")
    # Determine the smallest, non-zero n_times two chromophores are connected.
    I_min = np.min(intensities)
    # Determine the max number of times two chormophores are connected.
    I_max = np.max(intensities)

    # Set up the color bar.
    color_map = plt.get_cmap("Greys")
    c_norm = colors.Normalize(vmin=np.floor(I_min), vmax=np.ceil(I_max))
    scalar_map = cmx.ScalarMappable(norm=c_norm, cmap=color_map)
    hop_colors = scalar_map.to_rgba(intensities)

    n = len(hop_colors)
    # Plot the vectors between two chromophores
    # X-vals for start of each arrow
    # Y-vals for start of each arrow
    # Z-vals for start of each arrow
    # X-vecs for direction of each arrow
    # Y-vecs for direction of each arrow
    # Z-vecs for direction of each arrow
    ax.quiver(
        np.zeros(n),
        np.zeros(n),
        np.zeros(n),
        hop_vectors[:, 0],
        hop_vectors[:, 1],
        hop_vectors[:, 2],
        color=hop_colors,
        arrow_length_ratio=0,
        linewidth=0.7,
    )
    ax.set_xlim([-np.max(hop_vectors[:, 0]), np.max(hop_vectors[:, 0])])
    ax.set_ylim([-np.max(hop_vectors[:, 1]), np.max(hop_vectors[:, 1])])
    ax.set_zlim([-np.max(hop_vectors[:, 2]), np.max(hop_vectors[:, 2])])

    # Make the color bar
    scalar_map.set_array(intensities)
    # tick_location = np.arange(0, int(np.ceil(vmax)) + 1, 1)
    tick_location = np.arange(0, I_max, 1)
    # Plot the color bar
    cbar = plt.colorbar(scalar_map, ticks=tick_location, shrink=0.8, aspect=20)
    cbar.ax.set_yticklabels([rf"10$^{{{x}}}$" for x in tick_location])

    # Name and save the figure.
    if c_type == "hole":
        species = "donor"
    else:
        species = "acceptor"
    figure_title = f"{species.capitalize()} ({c_type.capitalize()}) Network"

    filename = f"hop_vec_{c_type}.png"
    filepath = (os.path.join(path, filename),)
    plt.savefig(filepath, bbox_inches="tight", dpi=300)
    print(f"\tFigure saved as {filename}")
    plt.clf()


def plot_anisotropy(carrier_data, c_type, three_d, path):  # pragma: no cover
    box = carrier_data["box"][0]
    xyzs = []
    # Get the indices of the carriers that travelled the furthest
    # only do the first 1000, in case there's a lot
    for i, pos in enumerate(carrier_data["current_position"][:1000]):
        image = carrier_data["image"][i]
        position = image * box + pos
        xyzs.append(position / 10.0)

    colors = ["b"] * len(xyzs)
    xyzs = np.array(xyzs)

    anisotropy = get_anisotropy(xyzs)
    print("----------------------------------------")
    print(f"{c_type} charge transport anisotropy: {anisotropy:.3f}")
    print("----------------------------------------")

    if three_d:
        # Reduce number of plot markers
        fig = plt.gcf()
        ax = p3.Axes3D(fig)
        plt.scatter(xyzs[:, 0], xyzs[:, 1], zs=xyzs[:, 2], c=colors, s=20)
        plt.scatter(0, 0, zs=0, c="r", s=50)

        # Draw boxlines
        box_pts = hf.box_points(box)
        ax.plot(box_pts[:, 0], box_pts[:, 1], box_pts[:, 2], c="k", linewidth=1)

        ax.set_xlabel("X (nm)", fontsize=20, labelpad=40)
        ax.set_ylabel("Y (nm)", fontsize=20, labelpad=40)
        ax.set_zlabel("Z (nm)", fontsize=20, labelpad=40)
        maximum = np.max(xyzs)
        ax.set_xlim([-maximum, maximum])
        ax.set_ylim([-maximum, maximum])
        ax.set_zlim([-maximum, maximum])
        ticks = [
            getattr(ax, f"{xyz}axis").get_major_ticks()
            for xyz in ["x", "y", "z"]
        ]
        ticks = [i for l in ticks for i in l]
        for tick in ticks:
            tick.label.set_fontsize(16)
        ax.dist = 11
        plt.title(f"Anisotropy ({c_type.capitalize()})", y=1.1)
        filename = f"anisotropy_{c_type}.png"
        filepath = os.path.join(path, filename)
        plt.savefig(filepath, bbox_inches="tight", dpi=300)
        plt.clf()
        print(f"\tFigure saved as {filename}")
    return anisotropy


def plot_temp_progression(
    temp, mobility, mob_err, anisotropy, c_type, path
):  # pragma: no cover
    plt.gcf()
    plt.clf()
    xvals = temp
    yvals = mobility
    yerrs = mob_err
    plt.xlabel("Temperature (Arb. U)")
    plt.ylabel(r"Mobility (cm$^{2}$ / Vs)")
    plt.errorbar(xvals, yvals, xerr=0, yerr=yerrs)
    plt.yscale("log")
    filename = f"mobility_{c_type}.png"
    filepath = os.path.join(path, filename)
    plt.savefig(filepath, dpi=300)
    plt.clf()
    print(f"\tFigure saved as {filename}")

    plt.plot(temp, anisotropy, c="r")
    plt.xlabel("Temperature (Arb. U)")
    plt.ylabel(r"$\kappa$ (Arb. U)")
    filename = f"anisotropy_{c_type}.png"
    filepath = os.path.join(path, filename)
    plt.savefig(filepath, dpi=300)
    plt.clf()
    print(f"\tFigure saved as {filename}")


def get_lambda_ij(chromo_length):
    # The equation for the internal reorganisation energy was obtained from
    # the data given in
    # Johansson, E and Larsson, S; 2004, Synthetic Metals 144: 183-191.
    # External reorganisation energy obtained from
    # Liu, T and Cheung, D. L. and Troisi, A; 2011, Phys. Chem. Chem. Phys.
    # 13: 21461-21470
    lambda_external = 0.11  # eV
    if chromo_length < 12:
        lambda_internal = 0.20826 - (chromo_length * 0.01196)
    else:
        lambda_internal = 0.06474
    lambdae_V = lambda_external + lambda_internal
    return lambdae_V


def gaussian(x, a, x0, sigma):
    return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))


def gauss_fit(data):
    mean = np.mean(data)
    std = np.std(data)
    hist, bin_edges = np.histogram(data, bins=100)
    try:
        fit_args, fit_conv = curve_fit(
            gaussian, bin_edges[:-1], hist, p0=[1, mean, std]
        )
    except RuntimeError:
        fit_args = None
    return bin_edges, fit_args, mean, std


def plot_neighbor_hist(
    chromo_list, chromo_mol_id, box, sepcut, path
):  # pragma: no cover
    seps_donor = []
    seps_acceptor = []
    for i, ichromo in enumerate(chromo_list):
        for j, img in ichromo.neighbors:
            jchromo = chromo_list[j]
            # Skip any chromophores that are part of the same molecule
            if chromo_mol_id[i] == chromo_mol_id[j]:
                continue
            sep = np.linalg.norm(jchromo.center + img * box - ichromo.center)
            if ichromo.species == "donor":
                seps_donor.append(sep)
            else:
                seps_acceptor.append(sep)
    species = ["donor", "acceptor"]
    seps = [seps_donor, seps_acceptor]
    for sp_i, sp in enumerate(species):
        sep = seps[sp_i]
        if len(sep) == 0:
            continue

        plt.figure()
        n, bin_edges, _ = plt.hist(sep, bins=40, color="b")
        bin_centers = (bin_edges[1:] + bin_edges[:-1]) / 2.0
        bin_centers = np.insert(bin_centers, 0, 0)
        n = np.insert(n, 0, 0)
        smooth_n = gaussian_filter(n, 1.0)
        plt.plot(bin_centers, smooth_n, color="r")
        if sepcut[sp_i] is None:
            sepcut[sp_i] = get_dist_cutoff(
                bin_centers, smooth_n, min_i=0, at_least=100
            )
        if sepcut[sp_i] is not None:
            print(f"{sp} separation cluster cut-off set to {sepcut[sp_i]}")
            plt.axvline(sepcut[sp_i], c="k")
        plt.xlabel(rf"{sp.capitalize()} r$_{{i,j}}$ ($\AA$)")
        plt.ylabel("Frequency (Arb. U.)")
        filename = f"neighbor_hist_{sp}.png"
        filepath = os.path.join(path, filename)
        plt.savefig(filepath, dpi=300)
        plt.close()
        print(f"Neighbor histogram figure saved as {filename}")
        plt.close()


def plot_orientation_hist(
    chromo_list, chromo_mol_id, orientations, ocut, path
):  # pragma: no cover
    orientations_donor = []
    orientations_acceptor = []
    for i, ichromo in enumerate(chromo_list):
        ivec = orientations[i]
        for j, img in ichromo.neighbors:
            jchromo = chromo_list[j]
            # Skip any chromophores that are part of the same molecule
            if chromo_mol_id[i] == chromo_mol_id[j]:
                continue
            jvec = orientations[j]
            dot_product = np.dot(ivec, jvec)
            # in radians
            angle = np.arccos(np.abs(dot_product))
            if ichromo.species == "donor":
                orientations_donor.append(angle)
            elif ichromo.species == "acceptor":
                orientations_acceptor.append(angle)
    species = ["donor", "acceptor"]
    orients = [orientations_donor, orientations_acceptor]
    for sp_i, sp in enumerate(species):
        orient = orients[sp_i]
        if len(orient) == 0:
            continue

        plt.figure()
        n, bin_edges, _ = plt.hist(orient, bins=40, color="b")
        bin_centers = (bin_edges[1:] + bin_edges[:-1]) / 2.0
        bin_centers = np.insert(bin_centers, 0, 0)
        n = np.insert(n, 0, 0)
        smooth_n = gaussian_filter(n, 1.0)
        plt.plot(bin_centers, smooth_n, color="r")
        if ocut[sp_i] is None:
            ocut[sp_i] = get_dist_cutoff(
                bin_centers, smooth_n, max_i=0, at_least=100
            )
        if ocut[sp_i] is not None:
            print(f"{sp} orientation cluster cut-off set to {ocut[sp_i]}")
            plt.axvline(ocuts[sp_i], c="k")
        plt.xlabel(f"{sp.capitalize()} Orientations (rad)")
        plt.xlim([0, np.deg2rad(90)])
        plt.xticks(np.arange(0, np.deg2rad(91), np.deg2rad(15)))
        plt.ylabel("Frequency (Arb. U.)")
        filename = f"orientation_hist_{sp}.png"
        filepath = os.path.join(path, filename)
        plt.savefig(filepath, dpi=300)
        plt.close()
        print(f"Orientation histogram figure saved as {filename}")


def create_cutoff_dict(sepcut, ocut, ticut, freqcut):
    cutoff_dict = {
        "separation": sepcut,
        "orientation": ocut,
        "ti": ticut,
        "freq": freqcut,
    }
    return cutoff_dict


def get_clusters(chromo_list, snap, rmax=None):
    clusters = []
    box = snap.configuration.box
    species = ["donor", "acceptor"]
    for sp_i, sp in enumerate(species):
        chromo_ids = [c.id for c in chromo_list if c.species == sp]
        positions = np.array([chromo_list[i].center for i in chromo_ids])
        print(f"Examining the {sp} material...")
        if len(positions) == 0:
            print("No material found. Continuing...")
            clusters.append(None)
            continue
        print("Calculating clusters...")
        if rmax is None:
            rmax = max(box) / 4
            print(f"No cutoff provided: cluster cutoff set to {rmax:.3f}")

        cl = freud.cluster.Cluster()

        cl.compute((box, positions), neighbors={"r_max": rmax})

        large = sum([1 for c in cl.cluster_keys if len(c) > 6])
        biggest = max([len(c) for c in cl.cluster_keys])
        psi = large / cl.num_clusters
        print("----------------------------------------")
        print(f"{sp.capitalize()}: Detected {cl.num_clusters} total")
        print(f"and {large} large clusters (size > 6).")
        print(f"Largest cluster size: {biggest} chromophores.")
        print(f'Ratio in "large" clusters: {psi:.2f}')
        print("----------------------------------------")
        clusters.append(cl)
    return clusters


def get_orientations(chromo_list, snap):
    orientations = []
    for chromo in chromo_list:
        positions = snap.particles.position[chromo.atom_ids]

        # There is a really cool way to do this with single value composition
        # but on the time crunch I didn't have time to learn how to implement
        # it properly. Check https://goo.gl/jxuhvJ for more details.
        plane = get_plane(positions)
        orientations.append(np.array(plane) / np.linalg.norm(plane))
    return orientations


def get_plane(positions):
    ## See https://goo.gl/jxuhvJ for details on this methodology.
    vec1 = hf.find_axis(positions[0], positions[1])
    vec2 = hf.find_axis(positions[0], positions[2])
    return np.cross(vec1, vec2)


def update_cluster(atom_ID, cluster_list, neighbor_dict): # pragma: no cover
    try:
        for neighbor in neighbor_dict[atom_ID]:
            if cluster_list[neighbor] > cluster_list[atom_ID]:
                cluster_list[neighbor] = cluster_list[atom_ID]
                cluster_list = update_cluster(
                    neighbor, cluster_list, neighbor_dict
                )
            elif cluster_list[neighbor] < cluster_list[atom_ID]:
                cluster_list[atom_ID] = cluster_list[neighbor]
                cluster_list = update_cluster(
                    neighbor, cluster_list, neighbor_dict
                )
    except KeyError:
        pass
    return cluster_list


def cluster_tcl_script(clusters, large_cluster, path): # pragma: no cover
    """
    Create a tcl script for each identified cluster.
    """
    # Obtain the IDs of the cluster sizes, sorted by largest first
    print("Sorting the clusters by size...")
    species = ["donor", "acceptor"]
    for i, cl in enumerate(clusters):
        sp = species[i]

        if cl is None:
            print(f"No clusters found for {sp}, continuing...")

        sorted_cl = sorted(cl.cluster_keys, key=lambda i: len(i), reverse=True)
        colors = list(range(int(1e6)))
        count = 0

        print("Creating tcl header...")
        tcl_text = ["mol delrep 0 0;"]
        tcl_text += ["pbc wrap -center origin;"]
        tcl_text += ["pbc box -color black -center origin -width 4;"]
        tcl_text += ["display resetview;"]
        tcl_text += ["color change rgb 9 1.0 0.29 0.5;"]
        tcl_text += ["color change rgb 16 0.25 0.25 0.25;"]

        if sp == "donor":
            for i, chromo_ids in enumerate(sorted_cl):
                print(f"Creating tcl for cluster {i+1}/{len(sorted_cl)}")
                if len(chromo_ids) > large_cluster:
                    # Only make clusters that are ``large''
                    inclust = ""
                    for chromo in chromo_ids:
                        inclust += f"{inclust}{chromo:d} "
                tcl_text += ["mol material AOEdgy;"]
                # Use AOEdgy if donor
                # The +1 makes the largest cluster red rather than
                # blue (looks better with AO, DoF, shadows)
                tcl_text += [f"mol color ColorID {colors[count + 1 % 32]:d};"]
                # VMD has 32 unique colors
                tcl_text += ["mol representation VDW 4.0 8.0;"]
                tcl_text += [f"mol selection resid {inclust:s};"]
                tcl_text += ["mol addrep 0;"]
                count += 1
        else:
            for i, chromo_ids in enumerate(sorted_cl):
                print(f"Creating tcl for cluster {i+1}/{len(sorted_cl)}")
                if len(chromo_ids) > large_cluster:
                    inclust = ""
                    for chromo in chromo_ids:
                        inclust += f"{inclust}{chromo:d} "
                tcl_text += ["mol material Glass2;"]
                # Use Glass2 if acceptor
                # The +1 makes the largest cluster red rather than blue (looks
                # better with AO, DoF, shadows)
                tcl_text += [f"mol color ColorID {colors[count + 1 % 32]:d};"]
                tcl_text += ["mol representation VDW 4.0 8.0;"]
                tcl_text += [f"mol selection resid {inclust};"]
                tcl_text += ["mol addrep 0;"]
                count += 1
    tcl_file_path = os.path.join(path, "{sp}_cluster_colors.tcl")
    with open(tcl_file_path, "w+") as tcl_file:
        tcl_file.writelines("".join(tcl_text))
    print(f"\nClusters coloring written to {tcl_file_path}")


def get_lists_for_3d_clusters(
        clusters,
        chromo_list,
        colors,
        large_cluster
): # pragma: no cover
    data = []
    species = ["donor", "acceptor"]
    for i, cl in enumerate(clusters):
        sp = species[i]

        if cl is None:
            print(f"No clusters found for {sp}, continuing...")

        if sp == "donor":
            for clust_id, chromo_ids in enumerate(cl.cluster_keys):
                if len(chromo_ids) > large_cluster:
                    for c_id in chromo_ids:
                        chromo = chromo_list[c_id]
                        data.append([*chromo.center, "w", colors[clust_id % 7]])
        else:
            for clust_id, chromo_ids in enumerate(cl.cluster_keys):
                if len(chromo_ids) > large_cluster:
                    for c_id in chromo_ids:
                        chromo = chromo_list[c_id]
                        data.append(
                            [*chromo.center, colors[clust_id % 7], "none"]
                        )

    # Needs a sort so that final image is layered properly
    data = list(sorted(data, key=lambda x: x[0]))
    # Split up list into sublists
    xyzs = np.array([row[:3] for row in data])

    face_colors = [row[3] for row in data]
    edge_colors = [row[4] for row in data]
    return xyzs, face_colors, edge_colors


def plot_clusters_3D(
    chromo_list, clusters, box, generate_tcl, path
):  # pragma: no cover
    fig = plt.figure()
    ax = p3.Axes3D(fig)
    colors = ["r", "g", "b", "c", "m", "y", "k"]
    large_cluster = 6
    if generate_tcl:
        cluster_tcl_script(clusters, large_cluster, path)

    xyzs, face_colors, edge_colors = get_lists_for_3d_clusters(
        clusters, chromo_list, colors, large_cluster
    )
    ax.scatter(
        xyzs[:, 0],
        xyzs[:, 1],
        xyzs[:, 2],
        facecolors=face_colors,
        edgecolors=edge_colors,
        alpha=0.6,
        s=40,
    )

    # Draw boxlines
    box_pts = hf.box_points(box)
    ax.plot(box_pts[:, 0], box_pts[:, 1], box_pts[:, 2], c="k", linewidth=1)

    ax.set_xlim([-box[0] / 2, box[0] / 2])
    ax.set_ylim([-box[1] / 2, box[1] / 2])
    ax.set_zlim([-box[2] / 2, box[2] / 2])

    filename = "clusters.png"
    filepath = os.path.join(path, "figures", filename)
    plt.savefig(filepath, bbox_inches="tight", dpi=300)
    print(f"3D cluster figure saved as {filename}")
    plt.close()


def plot_energy_levels(chromo_list, data_dict, path):  # pragma: no cover
    homo_levels = []
    lumo_levels = []
    donor_delta_eij = []
    acceptor_delta_eij = []
    donor_lambda_ij = None
    acceptor_lambda_ij = None
    for chromo in chromo_list:
        if chromo.species == "donor":
            homo_levels.append(chromo.homo)
            for i, delta_eij in enumerate(chromo.neighbors_delta_e):
                ti = chromo.neighbors_ti[i]
                if delta_eij is not None and ti is not None:
                    donor_delta_eij.append(delta_eij)
                donor_lambda_ij = chromo.reorganization_energy
        else:
            lumo_levels.append(chromo.lumo)
            for i, delta_eij in enumerate(chromo.neighbors_delta_e):
                ti = chromo.neighbors_ti[i]
                if delta_eij is not None and ti is not None:
                    acceptor_delta_eij.append(delta_eij)
                acceptor_lambda_ij = chromo.reorganization_energy
    if donor_delta_eij:
        (donor_bin_edges, donor_fit_args, donor_mean, donor_std) = gauss_fit(
            donor_delta_eij
        )
        donor_err = donor_std / np.sqrt(len(donor_delta_eij))

        data_dict["donor_delta_eij_mean"] = donor_mean
        data_dict["donor_delta_eij_std"] = donor_std
        data_dict["donor_delta_eij_err"] = donor_err

        homo_av = np.average(homo_levels)
        homo_std = np.std(homo_levels)
        homo_err = homo_std / np.sqrt(len(homo_levels))
        data_dict["donor_homo_mean"] = homo_av
        data_dict["donor_homo_std"] = homo_std
        data_dict["donor_homo_err"] = homo_err
        print(f"Donor HOMO Level = {homo_av:.3f} +/- {homo_err:.3f}")
        print(f"Donor Delta E_ij mean = {donor_mean:.3f} +/- {donor_err:.3f}")
        plot_delta_eij(
            donor_delta_eij,
            donor_bin_edges,
            donor_fit_args,
            "donor",
            donor_lambda_ij,
            path,
        )
    if acceptor_delta_eij:
        (
            acceptor_bin_edges,
            acceptor_fit_args,
            acceptor_mean,
            acceptor_std,
        ) = gauss_fit(acceptor_delta_eij)
        acceptor_err = acceptor_std / np.sqrt(len(acceptor_delta_eij))

        data_dict["acceptor_delta_eij_mean"] = acceptor_mean
        data_dict["acceptor_delta_eij_std"] = acceptor_std
        data_dict["acceptor_delta_eij_err"] = acceptor_err
        lumo_av = np.average(lumo_levels)
        lumo_std = np.std(lumo_levels)
        lumo_err = lumo_std / np.sqrt(len(lumo_levels))
        data_dict["acceptor_lumo_mean"] = LUMO_av
        data_dict["acceptor_lumo_std"] = LUMO_std
        data_dict["acceptor_lumo_err"] = LUMO_err
        print(
            f"Acceptor LUMO Level = {lumo_av} +/- {lumo_err}\n",
            f"Acceptor Delta E_ij mean = {acceptor_mean} +/-{acceptor_err}",
        )
        plot_delta_eij(
            acceptor_delta_eij,
            acceptor_bin_edges,
            acceptor_fit_args,
            "acceptor",
            acceptor_lambda_ij,
            path,
        )


def plot_delta_eij(
    delta_eij, gauss_bins, fit_args, species, lambda_ij, path
):  # pragma: no cover
    plt.figure()
    n, bins, patches = plt.hist(
        delta_eij,
        np.arange(np.min(delta_eij), np.max(delta_eij), 0.05),
        color=["b"],
    )
    if fit_args is not None:
        gauss_Y = gaussian(gauss_bins[:-1], *fit_args)
        scale_factor = max(n) / max(gauss_Y)
        plt.plot(gauss_bins[:-1], gauss_Y * scale_factor, "ro:")
    else:
        print("No Gaussian found (probably zero-width delta function)")
    if lambda_ij is not None:
        plt.axvline(-float(lambda_ij), c="k")
    plt.ylabel("Frequency (Arb. U.)")
    plt.xlabel(rf"{species.capitalize()} $\Delta E_{{i,j}}$ (eV)")

    filename = f"{species}_delta_E_ij.png"
    filepath = os.path.join(path, filename)
    plt.savefig(filepath, dpi=300)
    plt.close()
    print(f"\tFigure saved as {filename}")


def plot_mixed_hopping_rates(
    chromo_list,
    clusters,
    chromo_mol_id,
    data_dict,
    cutoff_dict,
    temp,
    path,
    use_vrh=False,
    koopmans=None,
    boltz=False,
):  # pragma: no cover
    # Create all the empty lists we need
    hop_types = ["intra", "inter"]
    hop_targets = ["c", "m"]
    # hop_properties = ["r", "T"]
    species = ["d", "a"]
    prop_lists = defaultdict(list)
    for i, ichromo in enumerate(chromo_list):
        imol = chromo_mol_id[i]
        for ineighbor, ti in enumerate(ichromo.neighbors_ti):
            if ti is None or ti == 0:
                continue
            j = ichromo.neighbors[ineighbor][0]
            jchromo = chromo_list[j]
            jmol = chromo_mol_id[j]
            delta_e = ichromo.neighbors_delta_e[ineighbor]
            if ichromo.species == jchromo.species:
                lambda_ij = ichromo.reorganization_energy
            else:
                lambda_ij = (
                    ichromo.reorganization_energy
                    + jchromo.reorganization_energy
                ) / 2
            # Now take into account the various behaviours we can have from
            # the parameter file
            prefactor = 1.0
            # Apply the koopmans prefactor
            if koopmans is not None:
                prefactor *= koopmans
            # Apply the distance penalty due to VRH
            if use_vrh:
                vrh = 1.0 / ichromo.vrh_delocalization
                rel_image = ichromo.neighbors[ineighbor][1]
                jchromo_center = jchromo.center + rel_image * box
                rij = np.linalg.norm(ichromo.center - jchromo_center) * 1e-10
                rate = hf.get_hop_rate(
                    lambda_ij,
                    ti,
                    delta_e,
                    prefactor,
                    temp,
                    use_vrh=use_vrh,
                    rij=rij,
                    vrh=vrh,
                    boltz=boltz,
                )
            else:
                rate = hf.get_hop_rate(
                    lambda_ij, ti, delta_e, prefactor, temp, boltz=boltz
                )
            if j < i:
                continue
            # Do intra- / inter- clusters
            if ichromo.species == "acceptor":
                if clusters[1] is not None:
                    for cluster in clusters[1].cluster_keys:
                        if i in cluster and j in cluster:
                            prop_lists["intra_cra"].append(rate)
                            prop_lists["intra_cTa"].append(ti)
                            break
                        else:
                            prop_lists["inter_cra"].append(rate)
                            prop_lists["inter_cTa"].append(ti)
            else:
                if clusters[0] is not None:
                    for cluster in clusters[0].cluster_keys:
                        if i in cluster and j in cluster:
                            prop_lists["intra_crd"].append(rate)
                            prop_lists["intra_cTd"].append(ti)
                            break
                        else:
                            prop_lists["inter_crd"].append(rate)
                            prop_lists["inter_cTd"].append(ti)
            # Now do intra- / inter- molecules
            if imol == jmol:
                if ichromo.species == "acceptor":
                    prop_lists["intra_mra"].append(rate)
                    prop_lists["intra_mTa"].append(ti)
                else:
                    prop_lists["intra_mrd"].append(rate)
                    prop_lists["intra_mTd"].append(ti)
            else:
                if ichromo.species == "acceptor":
                    prop_lists["inter_mra"].append(rate)
                    prop_lists["inter_mTa"].append(ti)
                else:
                    prop_lists["inter_mrd"].append(rate)
                    prop_lists["inter_mTd"].append(ti)
    # Donor cluster Plots:
    if prop_lists["intra_crd"]:
        val = prop_lists["intra_crd"]
        mean = np.mean(val)
        avg_std = np.std(val) / len(val)
        print(f"Mean intra-cluster donor rate: {mean:.3e}+/-{avg_std:.3e}")

    if prop_lists["inter_crd"]:
        val = prop_lists["inter_crd"]
        mean = np.mean(val)
        avg_std = np.std(val) / len(val)
        print(f"Mean inter-cluster donor rate: {mean:.3e}+/-{avg_std:.3e}")

    if prop_lists["intra_crd"] or prop_lists["inter_crd"]:
        plot_stacked_hist_rates(
            prop_lists["intra_crd"],
            prop_lists["inter_crd"],
            ["Intra-cluster", "Inter-cluster"],
            "donor",
            path,
        )
        plot_stacked_hist_tis(
            prop_lists["intra_cTd"],
            prop_lists["inter_cTd"],
            ["Intra-cluster", "Inter-cluster"],
            "donor",
            cutoff_dict["ti"][0],
            path,
        )
    # Acceptor cluster Plots:
    if prop_lists["intra_cra"]:
        val = prop_lists["intra_cra"]
        mean = np.mean(val)
        avg_std = np.std(val) / len(val)
        print(f"Mean intra-cluster acceptor rate: {mean:.3e}+/-{avg_std:.3e}")

    if prop_lists["inter_cra"]:
        val = prop_lists["inter_cra"]
        mean = np.mean(val)
        avg_std = np.std(val) / len(val)
        print(f"Mean inter-cluster acceptor rate: {mean:.3e}+/-{avg_std:.3e}")

    if prop_lists["intra_cra"] or prop_lists["inter_cra"]:
        plot_stacked_hist_rates(
            prop_lists["intra_cra"],
            prop_lists["inter_cra"],
            ["Intra-cluster", "Inter-cluster"],
            "acceptor",
            path,
        )
        plot_stacked_hist_tis(
            ["intra_cTa"],
            prop_lists["inter_cTa"],
            ["Intra-cluster", "Inter-cluster"],
            "acceptor",
            cutoff_dict["ti"][1],
            path,
        )
    # Donor Mol Plots:
    if prop_lists["intra_mrd"]:
        val = prop_lists["intra_mrd"]
        mean = np.mean(val)
        avg_std = np.std(val) / len(val)
        print(f"Mean intra-molecular donor rate: {mean:.3e}+/-{avg_std:.3e}")

    if prop_lists["inter_mrd"]:
        val = prop_lists["inter_mrd"]
        mean = np.mean(val)
        avg_std = np.std(val) / len(val)
        print(f"Mean inter-molecular donor rate: {mean:.3e}+/-{avg_std:.3e}")

    if prop_lists["intra_mrd"] or prop_lists["inter_mrd"]:
        plot_stacked_hist_rates(
            prop_lists["intra_mrd"],
            prop_lists["inter_mrd"],
            ["Intra-mol", "Inter-mol"],
            "donor",
            path,
        )
    # Acceptor Mol Plots:
    if prop_lists["intra_mra"]:
        val = prop_lists["intra_mra"]
        mean = np.mean(val)
        avg_std = np.std(val) / len(val)
        print(f"Mean intra-molecular acceptor rate: {mean:.3e}+/-{avg_std:.3e}")

    if prop_lists["inter_mra"]:
        val = prop_lists["inter_mra"]
        mean = np.mean(val)
        avg_std = np.std(val) / len(val)
        print(f"Mean inter-molecular acceptor rate: {mean:.3e}+/-{avg_std:.3e}")

    if prop_lists["intra_mra"] or prop_lists["inter_mra"]:
        plot_stacked_hist_rates(
            prop_lists["intra_mra"],
            prop_lists["inter_mra"],
            ["Intra-mol", "Inter-mol"],
            "acceptor",
            path,
        )
    # Update the dataDict
    for hop_type, target, sp in itertools.product(
        hop_types, hop_targets, species
    ):
        hop_name = f"{hop_type}_{target}r{sp}"
        n_hops = len(prop_lists[hop_name])
        if n_hops == 0:
            continue

        other_hop = hop_types[hop_types.index(hop_type) * -1 + 1]
        other_hop_name = f"{other_hop}_{target}r{sp}"

        total_hops = n_hops + len(prop_lists[other_hop_name])
        proportion = n_hops / total_hops

        mean_rate = np.mean(prop_lists[hop_name])
        stdev_rate = np.std(prop_lists[hop_name])
        data_dict[f"{hop_name}_hops"] = n_hops
        data_dict[f"{hop_name}_proportion"] = proportion
        data_dict[f"{hop_name}_rate_mean"] = mean_rate
        data_dict[f"{hop_name}_rate_std"] = stdev_rate


def plot_stacked_hist_rates(
    data1, data2, labels, species, path
):  # pragma: no cover
    plt.figure()
    n, bins, patches = plt.hist(
        [data1, data2],
        bins=np.logspace(1, 18, 40),
        stacked=True,
        color=["r", "b"],
        label=labels,
    )
    plt.ylabel("Frequency (Arb. U.)")
    plt.xlabel(rf"{species.capitalize()} k$_{{i,j}}$ (s$^{-1}$)")
    plt.xlim([1, 1e18])
    plt.xticks([1e0, 1e3, 1e6, 1e9, 1e12, 1e15, 1e18])
    plt.ylim([0, np.max(n) * 1.02])
    plt.legend(loc=2, prop={"size": 18})
    plt.gca().set_xscale("log")
    filename = f"{species}_hopping_rate_clusters.png"
    filepath = os.path.join(path, filename)
    plt.savefig(filepath, dpi=300)
    plt.close()
    print(f"\tFigure saved as {filename}")


def plot_stacked_hist_tis(
    data1, data2, labels, species, cutoff, path
):  # pragma: no cover
    plt.figure()
    n, bins, patches = plt.hist(
        [data1, data2],
        bins=np.linspace(0, 1.2, 20),
        stacked=True,
        color=["r", "b"],
        label=labels,
    )
    plt.ylabel("Frequency (Arb. U.)")
    plt.xlabel(rf"{species.capitalize()} J$_{{i,j}}$ (eV)")
    plt.ylim([0, np.max(n) * 1.02])
    if cutoff is not None:
        plt.axvline(float(cutoff), c="k")
    plt.legend(loc=0, prop={"size": 18})

    filename = f"{species}_transfer_integral_clusters.png"
    filepath = os.path.join(path, filename)
    plt.savefig(filepath, dpi=300)
    plt.close()
    print(f"\tFigure saved as {filename}")


def write_csv(data_dict, path): # pragma: no cover
    filepath = os.path.join(path, "results.csv")
    with open(filepath, "w") as f:
        w = csv.writer(f)
        for key, val in data_dict.items():
            w.writerow([key, val])
    print(f"\tCSV file written to {filepath}")


# TODO
def get_dist_cutoff(
    bin_centers, dist, min_i=None, max_i=None, at_least=100, log=False
):
    try:
        if min_i is not None:
            # Looking for minima
            minima = argrelextrema(dist, np.less)[0]
            if min_i < 0:
                # Sometimes a tiny minimum at super RHS breaks this
                cutoff = 0.0
                while True:
                    selected_min = minima[min_i]
                    cutoff = bin_centers[selected_min]
                    if dist[selected_min] > at_least:
                        break
                    min_i -= 1
            else:
                # Sometimes a tiny maximum at super LHS breaks this
                cutoff = 0.0
                while True:
                    selected_min = minima[min_i]
                    cutoff = bin_centers[selected_min]
                    if dist[selected_min] > at_least:
                        break
                    min_i += 1
        elif max_i is not None:
            # Looking for maxima
            maxima = argrelextrema(dist, np.greater)[0]
            if max_i < 0:
                # Sometimes a tiny maximum at super RHS breaks this
                cutoff = 0.0
                while True:
                    selected_max = maxima[max_i]
                    cutoff = bin_centers[selected_max]
                    if dist[selected_max] > at_least:
                        break
                    max_i -= 1
            else:
                # Sometimes a tiny maximum at super LHS breaks this
                cutoff = 0.0
                while True:
                    selected_max = maxima[max_i]
                    cutoff = bin_centers[selected_max]
                    if dist[selected_max] > at_least:
                        break
                    max_i += 1
        if log is True:
            cutoff = 10 ** cutoff
        return cutoff
    except IndexError:
        print("Notice: No minima found in distribution. Cutoff set to None.")
        return None


def plot_ti_hist(chromo_list, chromo_mol_id, ticut, path):  # pragma: no cover
    # ti_dist [[DONOR], [ACCEPTOR]]
    ti_intra = [[], []]
    ti_inter = [[], []]
    species = ["donor", "acceptor"]
    labels = ["Intra-mol", "Inter-mol"]
    for sp_i, sp in enumerate(species):
        for i, ichromo in enumerate(chromo_list):
            if ichromo.species != sp:
                continue
            for i_neighbor, (j, img) in enumerate(ichromo.neighbors):
                jchromo = chromo_list[j]
                if jchromo.species != sp:
                    continue
                if i >= j:
                    continue
                ti = ichromo.neighbors_ti[i_neighbor]
                if ti is not None:
                    if chromo_mol_id[i] == chromo_mol_id[j]:
                        ti_intra[sp_i].append(ti)
                    else:
                        ti_inter[sp_i].append(ti)
        if not (ti_intra[sp_i] and ti_inter[sp_i]):
            continue
        plt.figure()
        maxti = np.max(ti_intra[sp_i] + ti_inter[sp_i])
        n, bin_edges, _ = plt.hist(
            [ti_intra[sp_i], ti_inter[sp_i]],
            bins=np.linspace(0, maxti, 20),
            color=["r", "b"],
            stacked=True,
            label=labels,
        )
        bin_centers = (bin_edges[1:] + bin_edges[:-1]) / 2
        smooth_n = gaussian_filter(n[0] + n[1], 1.0)
        plt.plot(bin_centers, smooth_n, color="r")
        if ticut[sp_i] is None:
            ticut[sp_i] = get_dist_cutoff(
                bin_centers, smooth_n, min_i=-1, at_least=100
            )
        if ticut[sp_i] is not None:
            print(f"{sp} TI cluster cut-off set to {ticuts[sp_i]}")
            plt.axvline(ticut[sp_i], c="k")
        plt.xlim([0, np.max(ti_intra[sp_i] + ti_inter[sp_i])])
        plt.ylim([0, np.max(n) * 1.02])
        plt.ylabel("Frequency (Arb. U.)")
        plt.xlabel(rf"{sp.capitalize()} J$_{{i,j}}$ (eV)")
        plt.legend(loc=1, prop={"size": 18})
        filename = f"{sp}_transfer_integral_mols.png"
        filepath = os.path.join(path, filename)
        plt.savefig(filepath, dpi=300)
        plt.close()
        print(f"\tFigure saved as {filename}")


def plot_frequency_dist(
    c_type, carrier_history, freqcut, path
):  # pragma: no cover
    c_ind = ["hole", "electron"].index(c_type)
    nonzero = list(zip(*carrier_history.nonzero()))
    frequencies = []
    for i, j in nonzero:
        if i < j:
            # Only consider hops in one direction
            continue
        frequency = carrier_history[i, j]
        frequencies.append(np.log10(frequency))
    plt.figure()
    n, bin_edges, _ = plt.hist(frequencies, bins=60, color="b")
    bin_centers = (bin_edges[1:] + bin_edges[:-1]) / 2.0
    smooth_n = gaussian_filter(n, 1.0)
    plt.plot(bin_centers, smooth_n, color="r")
    if freqcut[c_ind] is None:
        print("DYNAMIC CUT")
        freqcut[c_ind] = get_dist_cutoff(
            bin_centers, smooth_n, min_i=-1, at_least=100, log=True
        )
        print(f"Cluster cut-off based on hop frequency set to {freqcut[c_ind]}")
    if freqcut[c_ind] is not None:
        plt.axvline(np.log10(freqcut[c_ind]), c="k")

    plt.xlabel(f"Total {c_type} hops (Arb. U.)")
    ax = plt.gca()
    tick_labels = np.arange(0, np.ceil(np.max(frequencies)) + 1, 1)
    plt.xlim([0, np.ceil(np.max(frequencies))])
    plt.xticks(tick_labels, [rf"10$^{{{int(x)}}}$" for x in tick_labels])
    plt.ylabel("Frequency (Arb. U.)")
    filename = f"total_hop_freq_{c_type}.png"
    filepath = os.path.join(path, filename)
    plt.savefig(filepath, dpi=300)
    plt.close()
    print(f"\tFigure saved as {filename}")


def plot_net_frequency_dist(c_type, carrier_history, path):  # pragma: no cover
    nonzero = list(zip(*carrier_history.nonzero()))
    frequencies = []
    for i, j in nonzero:
        if i < j:
            # Only consider hops in one direction
            continue
        frequency = np.abs(carrier_history[i, j] - carrier_history[j, i])
        if frequency > 0:
            frequencies.append(np.log10(frequency))
    plt.figure()
    plt.hist(frequencies, bins=60, color="b")
    plt.xlabel(f"Net {c_type} hops (Arb. U.)")
    ax = plt.gca()
    tick_labels = np.arange(0, np.ceil(np.max(frequencies)) + 1, 1)
    plt.xlim([0, np.ceil(np.max(frequencies))])
    plt.xticks(tick_labels, [rf"10$^{{{x}}}$" for x in tick_labels])
    plt.ylabel("Frequency (Arb. U.)")
    filename = f"net_hop_freq_{c_type}.png"
    filepath = os.path.join(path, filename)
    plt.savefig(filepath, dpi=300)
    plt.close()
    print(f"\tFigure saved as {filename}")


def plot_discrepancy_frequency_dist(
    c_type, carrier_history, path
):  # pragma: no cover
    nonzero = list(zip(*carrier_history.nonzero()))
    frequencies = []
    net_equals_total = 0
    net_near_total = 0
    for i, j in nonzero:
        if i < j:
            # Only consider hops in one direction
            continue
        total_hops = carrier_history[i, j] + carrier_history[j, i]
        net_hops = np.abs(carrier_history[i, j] - carrier_history[j, i])
        frequency = total_hops - net_hops
        if frequency == 0:
            net_equals_total += 1
            net_near_total += 1
        else:
            if frequency < 10:
                net_near_total += 1
            frequencies.append(np.log10(frequency))
    plt.figure()
    plt.hist(frequencies, bins=60, color="b")
    plt.xlabel("Discrepancy (Arb. U.)")
    ax = plt.gca()
    tick_labels = np.arange(0, np.ceil(np.max(frequencies)) + 1, 1)
    plt.xlim([0, np.ceil(np.max(frequencies))])
    plt.xticks(tick_labels, [rf"10$^{{{x}}}$" for x in tick_labels])
    plt.ylabel("Frequency (Arb. U.)")
    filename = f"hop_discrepancy_{c_type}.png"
    filepath = os.path.join(path, filename)
    plt.savefig(filepath, dpi=300)
    plt.close()
    print(f"There are {net_equals_total} paths with one-way transport.")
    print(f"There are {net_near_total} paths with total - net < 10.")
    print(f"\tFigure saved as {filename}")


def plot_mobility_msd(
    c_type, times, msds, time_stderr, msd_stderr, temp, path
):  # pragma: no cover
    # Create the first figure that will be replotted each time
    plt.figure()
    times, msds = hf.parallel_sort(times, msds)
    mobility, mob_error, r_squared = plot_msd(
        times, msds, time_stderr, msd_stderr, c_type, temp, path
    )
    print("----------------------------------------")
    print(
        f"{c_type.capitalize()} mobility = {mobility:.2E} ",
        f"+/- {mob_error:.2E} cm^2 V^-1 s^-1",
    )
    print("----------------------------------------")
    plt.close()
    return mobility, mob_error, r_squared


def carrier_plots(
    c_type, carrier_data, chromo_list, snap, freqcut, three_d, temp, path
): # pragma: no cover
    print(f"Considering the transport of {c_type}...")
    if c_type == "hole":
        carrier_history = carrier_data["hole_history"]
    else:
        carrier_history = carrier_data["electron_history"]

    print("Obtaining mean squared displacements...")
    times, msds, time_stderr, msd_stderr = get_times_msds(carrier_data)

    print(f"Plotting distribution of {c_type} displacements")
    plot_displacement_dist(carrier_data, c_type, path)

    print("Calculating mobility...")
    mobility, mob_error, r_squared = plot_mobility_msd(
        c_type, times, msds, time_stderr, msd_stderr, temp, path
    )

    print(f"Calculating {c_type} trajectory anisotropy...")
    anisotropy = plot_anisotropy(carrier_data, c_type, three_d, path)

    if three_d:
        print("Plotting hop vector distribution")
        plot_hop_vectors(carrier_data, chromo_list, snap, c_type, path)

        if carrier_history is not None:
            print("Determining carrier hopping connections...")
            plot_connections(chromo_list, carrier_history, c_type, path)

    print(f"Plotting {c_type} hop frequency distribution...")
    plot_frequency_dist(c_type, carrier_history, freqcut, path)

    print(f"Plotting {c_type} net hop frequency distribution...")
    plot_net_frequency_dist(c_type, carrier_history, path)

    print("Plotting (total - net hops) discrepancy distribution...")
    plot_discrepancy_frequency_dist(c_type, carrier_history, path)

    return anisotropy, mobility, mob_error, r_squared


def main(
    combined_data,
    temp,
    chromo_list,
    snap,
    path,
    three_d=False,
    freqcut=[None, None],
    sepcut=[None, None],
    ocut=[None, None],
    ticut=[None, None],
    generate_tcl=False,
    sequence_donor=None,
    sequence_acceptor=None,
    backend=None,
    use_vrh=False,
    koopmans=None,
    boltz=False,
):  # pragma: no cover
    # Load the matplotlib backend and the plotting subroutines
    global plt
    global p3
    if backend is not None:
        import matplotlib

        matplotlib.use(backend.strip())
    import matplotlib.pyplot as plt

    try:
        import mpl_toolkits.mplot3d as p3
    except ImportError:
        print("Could not import 3D plotting engine!")
        three_d = False

    # Create the figures path if it doesn't already exist
    fig_dir = os.path.join(path, "figures")
    os.makedirs(fig_dir, exist_ok=True)
    # Load in all the required data
    data_dict = {}
    print("---------- KMC_ANALYZE ----------")
    print(f"All figures saved in {fig_dir}")
    print("---------------------------------")

    box = snap.configuration.box[:3]
    hole_data, elec_data = split_carriers(combined_data)
    # Calculate the mobilities
    if hole_data["id"]:
        c_type = "hole"
        carrier_data = hole_data

        anisotropy, mobility, mob_error, r_squared = carrier_plots(
            c_type,
            carrier_data,
            chromo_list,
            snap,
            freqcut,
            three_d,
            temp,
            fig_dir,
        )

        data_dict[f"{c_type}_anisotropy"] = anisotropy
        data_dict[f"{c_type}_mobility"] = mobility
        data_dict[f"{c_type}_mobility_err"] = mob_error
        data_dict[f"{c_type}_mobility_r_squared"] = r_squared

    if elec_data["id"]:
        c_type = "electron"
        carrier_data = elec_data

        anisotropy, mobility, mob_error, r_squared = carrier_plots(
            c_type,
            carrier_data,
            chromo_list,
            snap,
            freqcut,
            three_d,
            temp,
            fig_dir,
        )

        data_dict[f"{c_type}_anisotropy"] = anisotropy
        data_dict[f"{c_type}_mobility"] = mobility
        data_dict[f"{c_type}_mobility_err"] = mob_error
        data_dict[f"{c_type}_mobility_r_squared"] = r_squared

    # Now plot the distributions!
    chromo_mol_id = get_molecule_ids(snap, chromo_list)

    plot_energy_levels(chromo_list, data_dict, fig_dir)

    orientations = get_orientations(chromo_list, snap)

    plot_neighbor_hist(chromo_list, chromo_mol_id, box, sepcut, fig_dir)
    plot_orientation_hist(
        chromo_list, chromo_mol_id, orientations, ocut, fig_dir
    )
    plot_ti_hist(chromo_list, chromo_mol_id, ticut, fig_dir)
    cutoff_dict = create_cutoff_dict(sepcut, ocut, ticut, freqcut)
    print("Cut-offs: ('value', [donor, acceptor])")
    print(*[f"\t{i}" for i in cutoff_dict.items()], sep="\n")

    clusters = get_clusters(chromo_list, snap, rmax=None)

    if three_d:
        print("Plotting 3D cluster location plot...")
        plot_clusters_3D(chromo_list, clusters, box, generate_tcl, path)

    plot_mixed_hopping_rates(
        chromo_list,
        clusters,
        chromo_mol_id,
        data_dict,
        cutoff_dict,
        temp,
        fig_dir,
        use_vrh=use_vrh,
        koopmans=koopmans,
        boltz=boltz,
    )

    print("Plotting cluster size distribution...")
    plot_cluster_size_dist(clusters, fig_dir)

    print("Writing CSV Output File...")
    write_csv(data_dict, path)

    print("Plotting Mobility and Anisotropy progressions...")
    if sequence_donor is not None:
        if data_dict["hole_anisotropy"]:
            plot_temp_progression(
                sequence_donor,
                data_dict["hole_mobility"],
                data_dict["hole_mobility_err"],
                data_dict["hole_anisotropy"],
                "hole",
                fig_dir,
            )
    if sequence_acceptor is not None:
        data_dict["electron_anisotropy"]
        data_dict["electron_mobility"]
        if data_dict["electron_anisotropy"]:
            plot_temp_progression(
                sequence_acceptor,
                data_dict["electron_mobility"],
                data_dict["electron_mobility_err"],
                data_dict["electron_anisotropy"],
                "electron",
                fig_dir,
            )
    else:
        print("Skipping plotting mobility evolution.")
