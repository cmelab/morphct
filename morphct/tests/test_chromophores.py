import numpy as np
import pytest

from base_test import BaseTest


class TestChromophores(BaseTest):
    def test_bad_init_chromophore(self, p3ht_snap):
        from morphct.chromophores import Chromophore, conversion_dict

        with pytest.raises(TypeError):
            Chromophore()

        atom_ids = np.array([1, 0, 4, 3, 2, 5, 6, 7, 8, 9, 10])
        with pytest.raises(TypeError):
            Chromophore(0, p3ht_snap, atom_ids, "bad_species", conversion_dict)

    def test_init_chromophore(self, p3ht_snap):
        from morphct.chromophores import Chromophore, conversion_dict

        atom_ids = np.array([1, 0, 4, 3, 2, 5, 6, 7, 8, 9, 10])
        chromo = Chromophore(0, p3ht_snap, atom_ids, "donor", conversion_dict)

        assert np.array_equal(chromo.atom_ids, atom_ids)
        assert np.allclose(
            chromo.center, np.array([8.53228682, -27.96646586, -35.07523658])
        )
        assert chromo.neighbors == chromo.dissociation_neighbors == []
        assert chromo.neighbors_ti == chromo.neighbors_delta_e == []
        assert chromo.homo == chromo.homo_1 == None
        assert chromo.lumo == chromo.lumo_1 == None
        assert chromo.id == 0
        assert np.array_equal(chromo.image, np.array([0, 0, 0]))
        assert chromo.n_atoms == 11
        assert len(chromo.qcc_input) == 1649
        assert isinstance(chromo.qcc_input, str)
        assert (
            chromo.qcc_input.split(";")[0]
            == "C 0.5295219754255722 0.9320423792455799 1.655608689036562"
        )
        assert chromo.reorganization_energy == 0.3064
        assert chromo.species == "donor"
        assert np.allclose(
            chromo.unwrapped_center,
            np.array([8.53228682, -27.96646586, -35.07523658]),
        )
        assert chromo.vrh_delocalization == 2e-10

    def test_chromos_from_smiles(self, p3ht_snap):
        from morphct.chromophores import get_chromo_ids_smiles, conversion_dict

        p3ht_smarts = "[#6]1[#6][#16][#6][#6]1CCCCCC"
        aaids = get_chromo_ids_smiles(p3ht_snap, p3ht_smarts, conversion_dict)
        assert len(aaids) == 29
        assert np.all([len(i) == 11 for i in aaids])

    def test_set_neighbors_voronoi(self, p3ht_snap, p3ht_chromo_list):
        from morphct.chromophores import set_neighbors_voronoi, conversion_dict

        box = p3ht_snap.configuration.box[:3]
        qcc_pairs = set_neighbors_voronoi(
            p3ht_chromo_list, p3ht_snap, conversion_dict, d_cut=min(box) / 2
        )
        pair, qcc_input = qcc_pairs[0]
        chromo = p3ht_chromo_list[0]

        assert (
            len(chromo.neighbors)
            == len(chromo.neighbors_ti)
            == len(chromo.neighbors_delta_e)
            == 14
        )
        assert len(qcc_pairs) == 181
        assert pair == (0, 1)
        assert (
            qcc_input.split(";")[0]
            == "C -3.7486698030653614 -0.26928230395247255 -1.6724275287954669"
        )

    def test_repr(self, p3ht_chromo_list):
        chromo = p3ht_chromo_list[0]
        assert (
            repr(chromo)
            == "Chromophore 0 (donor): 11 atoms at 8.532 -27.966 -35.075"
        )
