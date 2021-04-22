import os
import pytest


test_dir = os.path.dirname(__file__)


class BaseTest:
    @pytest.fixture
    def p3ht_snap(self):
        import gsd.hoomd

        filepath = os.path.join(test_dir, "assets/p3ht_2_15mers.gsd")
        with gsd.hoomd.open(name=filepath, mode="rb") as f:
            snap = f[0]

        return snap

    @pytest.fixture
    def p3ht_chromo_list(self):
        import pickle

        filepath = os.path.join(test_dir, "assets/chromo_list.pkl")
        with open(filepath, mode="rb") as f:
            chromo_list = pickle.load(f)
        return chromo_list

    @pytest.fixture
    def p3ht_chromo_list_neighbors(self):
        import pickle

        filepath = os.path.join(test_dir, "assets/chromo_list_neighbors.pkl")
        with open(filepath, mode="rb") as f:
            chromo_list = pickle.load(f)
        return chromo_list

    @pytest.fixture
    def p3ht_qcc_pairs(self):
        import pickle

        filepath = os.path.join(test_dir, "assets/qcc_pairs.pkl")
        with open(filepath, mode="rb") as f:
            qcc_pairs = pickle.load(f)
        return qcc_pairs

    @pytest.fixture
    def p3ht_sfilename(self):
        filepath = os.path.join(test_dir, "assets/singles_energies.txt")
        return filepath

    @pytest.fixture
    def p3ht_dfilename(self):
        filepath = os.path.join(test_dir, "assets/dimer_energies.txt")
        return filepath

    @pytest.fixture
    def p3ht_chromo_list_energies(self):
        import pickle

        filepath = os.path.join(test_dir, "assets/chromo_list_energies.pkl")
        with open(filepath, mode="rb") as f:
            chromo_list = pickle.load(f)
        return chromo_list

    @pytest.fixture
    def p3ht_combined_carriers(self):
        import pickle

        filepath = os.path.join(test_dir, "assets/combined_carriers.pkl")
        with open(filepath, mode="rb") as f:
            carriers = pickle.load(f)
        return carriers
