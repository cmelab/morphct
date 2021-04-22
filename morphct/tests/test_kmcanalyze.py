import numpy as np
import pytest

from base_test import BaseTest


class TestKMC(BaseTest):
    def test_split_carriers(self, p3ht_combined_carriers):
        from morphct.kmc_analyze import split_carriers

        holes,elecs = split_carriers(p3ht_combined_carriers)

        assert len(holes["id"]) == 20
        assert len(elecs["id"]) == 0

    def test_msds_mobility(self, p3ht_combined_carriers):
        from scipy.stats import linregress
        from morphct.kmc_analyze import calc_mobility, get_times_msds

        lts, msds, lt_std, msd_std = get_times_msds(p3ht_combined_carriers)

        assert lts == [1e-13, 1e-12]
        assert np.allclose(msds, [4.16367929357e-18, 5.4571824610e-18])
        assert np.allclose(lt_std, [3.35540746198e-16, 1.2514617195e-16])
        assert np.allclose(msd_std, [6.9193894210e-19, 4.6082781975e-19])

        temp = 300
        fit_X = np.linspace(np.min(lts), np.max(lts), 100)
        gradient, intercept, _, _, _ = linregress(lts, msds)
        lt_err = np.average(lt_std)
        msd_err = np.average(msd_std)
        fit_Y = (fit_X * gradient) + intercept
        mob, mob_err = calc_mobility(fit_X, fit_Y, lt_err, msd_err, temp)

        assert np.allclose(mob, 0.092657299520)
        assert np.allclose(mob_err, 0.041287982679)

    def test_get_connections(
            self, p3ht_chromo_list_energies, p3ht_combined_carriers
            ):
        from morphct.kmc_analyze import get_connections

        history = p3ht_combined_carriers["hole_history"]
        chromo_list = p3ht_chromo_list_energies
        box = np.array([85.18963, 85.18963, 85.18963])

        connections = get_connections(chromo_list, history, box)

        assert connections.shape == (25,7)
        assert np.allclose(connections[0][0], 8.53228681737)

    def test_get_anisotropy(self, p3ht_combined_carriers):
        from morphct.kmc_analyze import get_anisotropy

        box = p3ht_combined_carriers["box"][0]
        xyzs = []
        for i, pos in enumerate(p3ht_combined_carriers["current_position"]):
            image = p3ht_combined_carriers["image"][i]
            position = image * box + pos
            xyzs.append(position / 10.0)
        xyzs = np.array(xyzs)
        anisotropy = get_anisotropy(xyzs)

        assert np.allclose(anisotropy, 0.416638588839)

    def test_get_lambda_ij(self):
        from morphct.kmc_analyze import get_lambda_ij

        assert get_lambda_ij(10) == 0.19866
        assert get_lambda_ij(15) == 0.17474

    def test_gauss_fit(self, p3ht_chromo_list_energies):
        from morphct.kmc_analyze import gauss_fit

        chromo_list = p3ht_chromo_list_energies

        delta_es = []
        for chromo in chromo_list:
            if chromo.species == "donor":
                    for i, delta_e in enumerate(chromo.neighbors_delta_e):
                        ti = chromo.neighbors_ti[i]
                        if delta_e is not None and ti is not None:
                            delta_es.append(delta_e)

        bin_edges, fit_args, mean, std = gauss_fit(delta_es)

        assert len(bin_edges) == 101
        assert np.allclose(bin_edges[0], -0.450059)
        assert np.allclose(fit_args, [ 1.47390e+1, -4.500590e-3,  7.866818e-2])
        assert mean == 0
        assert np.allclose(std, 0.142674)

    def test_get_clusters(self, p3ht_chromo_list_energies, p3ht_snap):
        from morphct.kmc_analyze import get_clusters

        chromo_list = p3ht_chromo_list_energies
        cluster = get_clusters(chromo_list, p3ht_snap)[0]

        assert len(cluster.cluster_idx) == 30
        assert len(cluster.cluster_keys) == 2

    def test_get_orientations(self, p3ht_chromo_list_energies, p3ht_snap):
        from morphct.kmc_analyze import get_orientations

        chromo_list = p3ht_chromo_list_energies
        orients = get_orientations(chromo_list, p3ht_snap)

        assert len(orients) == 30
        assert np.allclose(orients[0], [0.42126867, 0.55094784,-0.720409])
