import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]


class OfficialScoringTests(unittest.TestCase):
    def score(self, data, genes, ctrl_size=50, seed=0):
        with tempfile.TemporaryDirectory() as directory:
            work = Path(directory)
            data.write_h5ad(work / "input.h5ad")
            environment = dict(os.environ, AUTONOMICS_INPUT0=str(work / "input.h5ad"),
                AUTONOMICS_OUTPUT0=str(work / "output.h5ad"), AUTONOMICS_OUTPUT1=str(work / "report.json"),
                SC_P_JSON_GENE_SETS=json.dumps({"score": genes}),
                SC_P_INT_CTRL_SIZE=str(ctrl_size), SC_P_INT_RANDOM_STATE=str(seed))
            result = subprocess.run([sys.executable, str(ROOT / "scripts/score_genes.py")],
                env=environment, capture_output=True, text=True, timeout=120)
            if result.returncode:
                self.assertFalse((work / "output.h5ad").exists())
                self.assertFalse((work / "report.json").exists())
                return result, None, None
            return result, ad.read_h5ad(work / "output.h5ad"), json.loads((work / "report.json").read_text())

    def fixture(self):
        generator = np.random.default_rng(72)
        return ad.AnnData(generator.poisson(3, (12, 100)).astype(float),
            obs=pd.DataFrame(index=[f"cell{index}" for index in range(12)]),
            var=pd.DataFrame(index=[f"G{index}" for index in range(100)]))

    def compare(self, data, genes, ctrl_size=50, seed=0):
        reference = data.copy()
        sc.tl.score_genes(reference, genes, score_name="score", ctrl_size=ctrl_size,
            random_state=seed, ctrl_as_ref=True, n_bins=25, gene_pool=None,
            copy=False, use_raw=data.raw is not None, layer=None)
        result, scored, report = self.score(data, genes, ctrl_size, seed)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(list(scored.obs_names), list(data.obs_names))
        self.assertEqual(list(scored.var_names), list(data.var_names))
        np.testing.assert_allclose(scored.obs["score"], reference.obs["score"], rtol=0, atol=1e-12)
        original = data.X.toarray() if sparse.issparse(data.X) else data.X
        retained = scored.X.toarray() if sparse.issparse(scored.X) else scored.X
        np.testing.assert_array_equal(retained, original)
        self.assertEqual(report["ctrl_size"], ctrl_size)
        self.assertFalse(report["scientific_signoff_granted"])
        return scored, report

    def test_dense_unclamped_control(self):
        self.compare(self.fixture(), ["G0", "G1"], ctrl_size=500)

    def test_sparse(self):
        data = self.fixture()
        data.X = sparse.csr_matrix(data.X)
        self.compare(data, ["G0", "G1"])

    def test_raw_axis(self):
        original = self.fixture()
        original.raw = original.copy()
        data = original[:, :40].copy()
        scored, report = self.compare(data, ["G80", "G81"])
        self.assertEqual(report["matrix_axis"], "raw")
        np.testing.assert_array_equal(scored.raw.X, data.raw.X)

    def test_missing_and_duplicate_members(self):
        self.compare(self.fixture(), ["G0", "G0", "G1", "MISSING"])

    def test_nonfinite_all_cell_summary(self):
        data = self.fixture()
        data.X[0, :] = np.nan
        scored, report = self.compare(data, ["G0", "G1"])
        self.assertTrue(np.isnan(scored.obs["score"].iloc[0]))
        self.assertEqual(report["gene_sets"]["score"]["status"], "NONFINITE_OFFICIAL_SCORES")
        self.assertIsNone(report["gene_sets"]["score"]["mean"])

    def test_no_controls_is_failure(self):
        result, scored, report = self.score(self.fixture(), [f"G{index}" for index in range(100)])
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("No control genes found", result.stderr)

    def test_axis_collision_is_failure(self):
        data = self.fixture()
        data.obs["score"] = 0
        result, scored, report = self.score(data, ["G0", "G1"])
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("collides with obs", result.stderr)


if __name__ == "__main__":
    unittest.main()
