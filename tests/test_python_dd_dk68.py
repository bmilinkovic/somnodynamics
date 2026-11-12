"""End-to-end check of the Python DK68 adapter on synthetic source epochs."""

import argparse
import importlib.util
import tempfile
import unittest
from pathlib import Path

import numpy as np
from scipy.io import loadmat, savemat


SCRIPT = (Path(__file__).resolve().parents[1] / "code" /
          "standalone_sleep_dk5_dk68_lrt30_current_epoching_2s256" /
          "ssdi" / "run_python_dd_dk68.py")
SPEC = importlib.util.spec_from_file_location("run_python_dd_dk68", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class PythonDDDK68Tests(unittest.TestCase):
    def test_source_to_compatible_result(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source_dir = root / "source/EPCTL01/mat"
            source_dir.mkdir(parents=True)
            rng = np.random.default_rng(42)
            data = np.zeros((3, 180, 2))
            dynamics = np.array([[.7, .12, 0], [0, .6, .08], [0, 0, .5]])
            for trial in range(2):
                for t in range(1, data.shape[1]):
                    data[:, t, trial] = dynamics @ data[:, t-1, trial] + rng.normal(size=3)
            source = source_dir / "EPCTL01_W__source_time_series_68.mat"
            savemat(source, {"data": data, "labels": np.array(["R1", "R2", "R3"],
                                                               dtype=object)})
            jobs = MODULE.discover(root / "source")
            self.assertEqual([(j[0], j[1]) for j in jobs], [("EPCTL01", "W_")])
            normalised, labels = MODULE.read_source(source, 3)
            np.testing.assert_allclose(normalised.mean(axis=(1, 2)), 0, atol=1e-12)
            np.testing.assert_allclose(normalised.std(axis=(1, 2), ddof=1), 1,
                                       atol=1e-12)
            self.assertEqual(labels, ["R1", "R2", "R3"])

            args = argparse.Namespace(
                output_root=root / "results", required_regions=3, dimensions=[1],
                pmax=2, criterion="LRT", regmode="LWR", fs=64, fres=16,
                runs=2, pre_iters=4, final_iters=4, workers=1, seed=0,
                overwrite=False)
            bands = {"delta": (2, 12), "theta": (14, 24)}
            path = MODULE.run_one("EPCTL01", "W_", source, args, bands)
            result = loadmat(path, simplify_cells=True)["results"]
            self.assertTrue(np.isfinite(np.asarray(result["delta_dopto_m1"])).all())
            self.assertIn("delta_gopto_m1", result)
            self.assertIn("theta_dopto_m1", result)
            self.assertEqual(int(result["meta"]["morder"]), 1)
            self.assertTrue((args.output_root / "logs/status/EPCTL01_W_.json").exists())
            self.assertEqual(MODULE.run_one("EPCTL01", "W_", source, args, bands),
                             path)
            checkpoint = MODULE.checkpoint_directory("EPCTL01", "W_", source,
                                                     args, bands)
            model_file = checkpoint / "model.npz"
            retained = checkpoint / "delta_m1.npz"
            interrupted = checkpoint / "theta_m1.npz"
            timestamps = (model_file.stat().st_mtime_ns, retained.stat().st_mtime_ns)
            original_values = np.asarray(result["delta_dopto_m1"])
            path.unlink()
            interrupted.unlink()
            resumed = MODULE.run_one("EPCTL01", "W_", source, args, bands)
            np.testing.assert_array_equal(
                np.asarray(loadmat(resumed, simplify_cells=True)["results"]["delta_dopto_m1"]),
                original_values)
            self.assertEqual(timestamps, (model_file.stat().st_mtime_ns,
                                          retained.stat().st_mtime_ns))
            self.assertTrue(interrupted.exists())

            compact_root = root / "compact"
            compact_file = compact_root / "W_/EPCTL01_var_model.npz"
            MODULE.atomic_npz(
                compact_file, A=np.asarray(result["model"]["A"]).reshape(3, 3, -1),
                V=result["model"]["V"],
                order=np.asarray(result["meta"]["morder"]),
                labels=np.asarray(labels, dtype=str),
                subject=np.asarray("EPCTL01"), stage=np.asarray("W_"),
                fs=np.asarray(64.0), pmax=np.asarray(2),
                criterion=np.asarray("LRT"), regmode=np.asarray("LWR"))
            self.assertEqual(len(MODULE.discover_models(compact_root)), 1)
            args.output_root = root / "prefit_results"
            prefit = MODULE.run_one("EPCTL01", "W_", compact_file, args, bands)
            np.testing.assert_allclose(
                np.asarray(loadmat(prefit, simplify_cells=True)["results"]["delta_dopto_m1"]),
                original_values, rtol=0, atol=0)


if __name__ == "__main__":
    unittest.main()
