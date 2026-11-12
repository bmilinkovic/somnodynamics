"""Small synthetic checks for the public, data-independent workflow."""

import unittest

import numpy as np

from multiscale import run_grid_parallel, run_time_changing_dd, minimum_dd_by_window
from multiscale.model import (iss2cak, prepare_var, ss2trfun, var_to_cak,
                              var_to_ss, var_to_trfun)


class CompanionCalculationTests(unittest.TestCase):
    def test_block_elimination_matches_dense_companion(self):
        rng = np.random.default_rng(91)
        a = rng.normal(scale=0.025, size=(4, 4, 3))
        ass, c, k = var_to_ss(a)
        np.testing.assert_allclose(var_to_cak(a), iss2cak(ass, c, k),
                                   rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(var_to_trfun(a, 16),
                                   ss2trfun(ass, c, k, 16),
                                   rtol=1e-12, atol=1e-12)
        prepared = prepare_var(a, np.eye(4), 128, 16)
        np.testing.assert_allclose(prepared[3], iss2cak(ass, c, k),
                                   rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(prepared[4], ss2trfun(ass, c, k, 16),
                                   rtol=1e-12, atol=1e-12)

    def test_selected_grid_block_keeps_full_grid_seed(self):
        a = np.zeros((3, 3, 1))
        a[:, :, 0] = np.array([[.5, .1, 0], [0, .4, .1], [0, 0, .3]])
        v = np.eye(3)
        kwargs = dict(fs=64, fres=16, dimensions=(1,),
                      bands={"low": (2, 8), "high": (10, 20)},
                      runs=2, seed=17, workers=1, pre_iters=4, final_iters=4)
        full = run_grid_parallel(a, v, **kwargs)
        selected = run_grid_parallel(a, v, pairs=[(1, "high")], **kwargs)
        np.testing.assert_allclose(selected[(1, "high")]["final"].values,
                                   full[(1, "high")]["final"].values,
                                   rtol=0, atol=0)


class TimeChangingDDTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(8)
        self.signals = np.zeros((3, 320))
        a = np.array([[.6, .12, 0], [0, .5, .1], [0, 0, .4]])
        for t in range(1, 320):
            self.signals[:, t] = a @ self.signals[:, t-1] + rng.normal(size=3)
        self.windows = {"early": (0, 160), "late": (160, 320)}
        self.common = dict(fs=64, fres=16, dimensions=[1], pmax=2,
                           criterion="aic", runs=2, pre_iters=4,
                           final_iters=4, state_space_iters=3)

    def test_both_domains_use_each_window(self):
        for domain, extra in (("temporal", {}),
                              ("spectral", {"bands": {"low": (2, 12)}})):
            with self.subTest(domain=domain):
                result = run_time_changing_dd(
                    self.signals, self.windows, domain=domain,
                    **self.common, **extra)
                rows = minimum_dd_by_window(result, domain=domain)
                self.assertEqual([row["window"] for row in rows], ["early", "late"])
                self.assertTrue(all(np.isfinite(row["dd_min"]) for row in rows))
                self.assertTrue(all(row["model_order"] > 0 for row in rows))

    def test_rejects_invalid_windows(self):
        with self.assertRaisesRegex(ValueError, "too short"):
            run_time_changing_dd(self.signals, {"short": (0, 3)},
                                 domain="temporal", **self.common)


if __name__ == "__main__":
    unittest.main()
