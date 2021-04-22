import numpy as np
import pytest

from base_test import BaseTest


class TestKMC(BaseTest):
    def test_init_carrier(self, p3ht_chromo_list_energies):
        from morphct.mobility_kmc import Carrier

        chromo_list = p3ht_chromo_list_energies
        n = len(chromo_list)
        chromo = chromo_list[0]
        box = np.array([85.18963, 85.18963, 85.18963])
        carrier = Carrier(chromo, 1e-13, 0, box, 300, n)

        assert carrier.displacement == 0

        np.random.seed(42)
        carrier.calculate_hop(chromo_list)

        assert carrier.n_hops == 1
        assert carrier.current_chromo.id == 1
        assert carrier.current_time == 7.685088279604407e-15

        carrier.update_displacement()

        assert carrier.displacement == 9.193875570863462

    def test_carrier_errors(self, p3ht_chromo_list_energies):
        from morphct.mobility_kmc import Carrier

        chromo_list = p3ht_chromo_list_energies
        n = len(chromo_list)
        chromo = chromo_list[0]
        box = np.array([85.18963, 85.18963, 85.18963])

        with pytest.raises(ValueError):
            Carrier(chromo, 1e-12, 0, box, 300, n, use_avg_hoprates=True)

        with pytest.warns(UserWarning):
            Carrier(chromo, 1e-12, 0, box, 300, n, avg_inter_rate=0)

        with pytest.raises(ValueError):
            Carrier(
                    chromo,
                    1e-12,
                    0,
                    box,
                    300,
                    n,
                    avg_inter_rate=0,
                    avg_intra_rate=0,
                    use_avg_hoprates=True
                    )

    def test_getjobslist(self):
        from morphct.mobility_kmc import get_jobslist

        n = 4
        lts = [1.0e-13, 1.0e-12]
        jobs = get_jobslist(lts, n_holes=10, n_elec=10, nprocs=n, seed=42)

        assert len(jobs) == n
        assert jobs[0][0] == [9, 1e-13, 'electron']

    def test_runsinglekmc(self, tmpdir, p3ht_chromo_list_energies, p3ht_snap):
        from morphct.mobility_kmc import run_single_kmc

        jobs = [[0, 1e-12, 'hole']]
        chromo_list = p3ht_chromo_list_energies
        box = np.array([85.18963, 85.18963, 85.18963])

        carrier = run_single_kmc(
                jobs, tmpdir, chromo_list, p3ht_snap, 300, seed=42
                )[0]

        assert carrier.n_hops == 212
        assert carrier.current_chromo.id == 14

        c_kwargs = {
            "use_avg_hoprates": True,
            "avg_intra_rate":2.243e+14,
            "avg_inter_rate":1.681e+11,
        }
        carrier = run_single_kmc(
                jobs,
                tmpdir,
                chromo_list,
                p3ht_snap,
                300,
                seed=42,
                carrier_kwargs=c_kwargs
                )[0]

        assert carrier.n_hops == 1553
        assert carrier.current_chromo.id == 19

        jobs = [[0, 1e-11, 'hole']]
        c_kwargs = {
                "use_koopmans" : True,
                "boltz": True,
                "use_vrh": True
                }

        carrier = run_single_kmc(
                jobs,
                tmpdir,
                chromo_list,
                p3ht_snap,
                300,
                seed=42,
                carrier_kwargs=c_kwargs
                )[0]

        assert carrier.n_hops == 975
        assert carrier.current_chromo.id == 11
