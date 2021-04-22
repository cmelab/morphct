import numpy as np
import pytest

from base_test import BaseTest


class TestHelperFunctions(BaseTest):
    def test_box_points(self):
        from morphct.helper_functions import box_points

        box = np.array([10, 10, 10])
        assert np.array_equal(
            box_points(box),
            np.array(
                [
                    [-0.5, -0.5, -0.5],
                    [-0.5, -0.5, 0.5],
                    [-0.5, 0.5, 0.5],
                    [-0.5, -0.5, 0.5],
                    [0.5, -0.5, 0.5],
                    [0.5, -0.5, -0.5],
                    [-0.5, -0.5, -0.5],
                    [-0.5, 0.5, -0.5],
                    [-0.5, 0.5, 0.5],
                    [0.5, 0.5, 0.5],
                    [0.5, -0.5, 0.5],
                    [0.5, -0.5, -0.5],
                    [0.5, 0.5, -0.5],
                    [-0.5, 0.5, -0.5],
                    [0.5, 0.5, -0.5],
                    [0.5, 0.5, 0.5],
                ]
            ),
        )

    def test_time_units(self):
        from morphct.helper_functions import time_units

        assert time_units(10) == "10.00 seconds"
        assert time_units(1000, precision=0) == "17 minutes"
        assert time_units(10000) == "2.78 hours"
        assert time_units(100000) == "1.16 days"

    def test_find_axis(self):
        from morphct.helper_functions import find_axis

        assert np.array_equal(
            find_axis(np.array([0, 0, 0]), np.array([1, 0, 0])),
            np.array([1, 0, 0]),
        )
        assert np.allclose(
            find_axis(np.array([0, 0, 0]), np.array([1, 1, 1])),
            np.array([0.57735027, 0.57735027, 0.57735027]),
        )
        assert np.array_equal(
            find_axis(
                np.array([0, 0, 0]), np.array([1, 1, 1]), normalize=False
            ),
            np.array([1, 1, 1]),
        )
        assert np.array_equal(
            find_axis(np.array([0, 0, 0]), np.array([0, 0, 0])),
            np.array([0, 0, 0]),
        )

    def test_parallel_sort(self):
        from morphct.helper_functions import parallel_sort

        list1 = [3, 2, 5, 1]
        list2 = ["A", "B", "C", "D"]

        list1, list2 = parallel_sort(list1, list2)

        assert list1 == [1, 2, 3, 5]
        assert list2 == ["D", "B", "A", "C"]

        list1 = np.array([3, 2, 5, 1])
        list2 = ("A", "B", "C", "D")

        list1, list2 = parallel_sort(list1, list2)

        assert np.array_equal(list1, np.array([1, 2, 3, 5]))
        assert list2 == ("D", "B", "A", "C")

    def get_hop_rate(self):
        from morphct.helper_functions import get_hop_rate

        lamda = 0.3064
        ti = 0
        delta = 0.06218457533310762
        factor = 1
        temp = 300

        # ti is zero so hop rate should be zero
        assert get_hop_rate(lamda, ti, delta, factor, temp) == 0

        ti = 0.2456720694088973
        delta = 0.016112646653095197

        assert get_hop_rate(lamda, ti, delta, factor, temp) == pytest.approx(
            68518361827044, 1
        )

        ti = 0.0013270585750558073
        delta = -0.021286665397703075
        rij = 3.659072209672184e-09
        vrh = 2e-10

        assert get_hop_rate(
            lamda,
            ti,
            delta,
            factor,
            temp,
            use_vrh=True,
            rij=rij,
            vrh=2e-10,
            boltz=True,
        ) == pytest.approx(603.98144350, 1e-8)

    def get_event_tau(self):
        from morphct.helper_functions import get_event_tau

        rate = 0

        assert get_event_tau(rate) == 1e99

        np.random.seed(42)
        rate = 10

        assert get_event_tau(rate) == pytest.approx(0.098205635, 1e-8)

        slowest = 3
        fastest = 0.5
        max_attempts = 100
        assert get_event_tau(
            rate, slowest=slowest, fastest=fastest, max_attempts=max_attempts
        ) == pytest.approx(0.519899395, 1e-8)
