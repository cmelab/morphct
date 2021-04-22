import numpy as np
from base_test import BaseTest


class TestEQCC(BaseTest):
    def test_gets_homolumo(self, p3ht_chromo_list):
        from morphct.execute_qcc import get_homolumo

        chromo = p3ht_chromo_list[0]
        data = get_homolumo(chromo.qcc_input)
        assert np.allclose(
            data, np.array([[-9.01337182, -8.5404688, 0.17193304, 0.86523495]])
        )

    def test_singles_homolumo(self, p3ht_chromo_list):
        from morphct.execute_qcc import singles_homolumo

        chromo = p3ht_chromo_list[0]
        data = singles_homolumo([chromo])
        assert np.allclose(
            data, np.array([[-9.01337182, -8.5404688, 0.17193304, 0.86523495]])
        )

    def test_dimer_homolumo(self, p3ht_qcc_pairs):
        from morphct.execute_qcc import dimer_homolumo

        qcc_pair = p3ht_qcc_pairs[0]
        pair, energies = dimer_homolumo([qcc_pair])[0]

        assert pair == (0, 1)
        assert np.allclose(
            energies,
            np.array([-8.70917208, -8.21756382, -0.26434042, 0.35184321]),
        )

    def test_set_energyvalues(
        self, p3ht_chromo_list_neighbors, p3ht_sfilename, p3ht_dfilename
    ):
        from morphct.execute_qcc import set_energyvalues

        set_energyvalues(
            p3ht_chromo_list_neighbors, p3ht_sfilename, p3ht_dfilename
        )

        chromo = p3ht_chromo_list_neighbors[0]

        assert chromo.homo_1 == -9.013371824203956
        assert chromo.homo == -8.54046879758057
        assert chromo.lumo == 0.1719330393321203
        assert chromo.lumo_1 == 0.8652349542720108
        assert chromo.neighbors_delta_e[0] == -0.016112646653095197
        assert chromo.neighbors_ti[0] == 0.2456720694088973
