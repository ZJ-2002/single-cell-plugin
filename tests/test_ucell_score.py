"""Golden tests for the h5ad_ucell_score node against official UCell.

Expected values in ucell_golden.json were produced by the official UCell R
package (2.16.0) via tests/generate_ucell_golden.R -- never by this plugin.
Regeneration: python tests/make_ucell_fixtures.py && Rscript
tests/generate_ucell_golden.R (see tests/README.md).
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

TESTS_DIR = Path(__file__).resolve().parent
PLUGIN_ROOT = TESTS_DIR.parent
sys.path.insert(0, str(PLUGIN_ROOT))
sys.path.insert(0, str(TESTS_DIR))

import workflow  # noqa: E402
from ucell_fixtures import all_fixtures  # noqa: E402

GOLDEN_PATH = TESTS_DIR / "ucell_golden.json"
EPSILON = 1e-14


def load_golden() -> dict[str, Any]:
    payload = json.loads(GOLDEN_PATH.read_text(encoding="utf-8"))
    assert payload["_provenance"]["package"] == "UCell"
    assert payload["_provenance"]["entrypoint"] == "ScoreSignatures_UCell"
    return payload


def build_adata(fixture) -> ad.AnnData:
    return ad.AnnData(
        X=sp.csr_matrix(fixture.as_array()),
        obs=pd.DataFrame(index=pd.Index(fixture.cells, name="cell_id")),
        var=pd.DataFrame(index=pd.Index(fixture.genes, name="gene_id")),
    )


def run_ucell(fixture, tmp_path: Path, extra_params: dict[str, Any] | None = None):
    """Drive workflow.ucell_score exactly as the container does: env contract
    files plus the params.json payload the run_workflow.py shim builds."""
    input_path = tmp_path / "input.h5ad"
    output_path = tmp_path / "output.h5ad"
    report_path = tmp_path / "report.json"
    build_adata(fixture).write_h5ad(input_path, compression="lzf")
    params: dict[str, Any] = {"gene_sets": fixture.gene_sets}
    params.update(fixture.params)
    if extra_params:
        params.update(extra_params)
    params_path = tmp_path / "params.json"
    params_path.write_text(json.dumps(params), encoding="utf-8")

    environ = {
        "AUTONOMICS_INPUT0": str(input_path),
        "AUTONOMICS_OUTPUT0": str(output_path),
        "AUTONOMICS_OUTPUT1": str(report_path),
        "AUTONOMICS_SINGLE_CELL_PARAMS": str(params_path),
    }
    original = {name: os.environ.get(name) for name in environ}
    os.environ.update(environ)
    try:
        adata, h5ad_out, report_out = workflow.ucell_score(params)
    finally:
        for name, value in original.items():
            if value is None:
                os.environ.pop(name, None)
            else:
                os.environ[name] = value
    report = json.loads(report_out.read_text(encoding="utf-8"))
    return adata, report


class TestAgainstOfficialUCell:
    @pytest.mark.parametrize("fixture", all_fixtures(), ids=lambda f: f.name)
    def test_matches_official_ucell(self, fixture, tmp_path: Path) -> None:
        golden = load_golden()["fixtures"][fixture.name]
        adata, report = run_ucell(fixture, tmp_path)
        assert list(adata.obs_names) == golden["cells"]
        for name, expected in golden["scores"].items():
            assert name in adata.obs, f"missing obs column {name}"
            actual = adata.obs[name].to_numpy(dtype=float)
            assert actual.shape == (len(golden["cells"]),)
            deviations = np.abs(actual - np.asarray(expected))
            worst = float(deviations.max())
            assert worst <= EPSILON, f"{fixture.name}/{name}: max |delta| = {worst}"

    def test_audit_counterexample_is_exactly_one(self, tmp_path: Path) -> None:
        """Workorder F01 counterexample: 1 cell, 10 genes, top-2 signature.
        The pre-fix rational heuristic returned 0.7333333; official UCell
        with maxRank=10 returns 1.0."""
        from ucell_fixtures import fixture_by_name

        adata, _ = run_ucell(fixture_by_name("audit_counterexample"), tmp_path)
        assert adata.obs["top2"].to_numpy(dtype=float)[0] == 1.0

    def test_full_dispatch_roundtrip(self, tmp_path: Path) -> None:
        """End to end through workflow.run(): env dispatch, params.json,
        H5AD output written back, report next to it."""
        from ucell_fixtures import fixture_by_name

        fixture = fixture_by_name("missing_skip")
        input_path = tmp_path / "input.h5ad"
        build_adata(fixture).write_h5ad(input_path, compression="lzf")
        params_path = tmp_path / "params.json"
        params_path.write_text(
            json.dumps({"gene_sets": fixture.gene_sets, **fixture.params}),
            encoding="utf-8",
        )
        environ = {
            workflow.WORKFLOW_ENV: "ucell_score",
            workflow.PARAMS_PATH: str(params_path),
            "AUTONOMICS_INPUT0": str(input_path),
            "AUTONOMICS_OUTPUT0": str(tmp_path / "output.h5ad"),
            "AUTONOMICS_OUTPUT1": str(tmp_path / "report.json"),
        }
        original = {name: os.environ.get(name) for name in environ}
        os.environ.update(environ)
        try:
            assert workflow.run() is None
        finally:
            for name, value in original.items():
                if value is None:
                    os.environ.pop(name, None)
                else:
                    os.environ[name] = value
        roundtrip = ad.read_h5ad(tmp_path / "output.h5ad")
        assert roundtrip.uns["ucell_params"]["missing_genes"] == "skip"
        golden = load_golden()["fixtures"]["missing_skip"]
        actual = roundtrip.obs["with_absent"].to_numpy(dtype=float)
        assert float(np.abs(actual - np.asarray(golden["scores"]["with_absent"])).max()) <= EPSILON


class TestParamsAndReport:
    def test_report_records_algorithm_and_params(self, tmp_path: Path) -> None:
        from ucell_fixtures import fixture_by_name

        fixture = fixture_by_name("truncated_ranks")
        adata, report = run_ucell(fixture, tmp_path)
        assert report["algorithm"] == "ucell"
        assert report["algorithm_version"] == "UCell 2.16.0"
        assert report["max_rank"] == 6
        assert report["effective_max_rank"] == 6
        assert report["ties_method"] == "average"
        assert report["missing_genes"] == "impute"
        assert report["operation"] == "ucell_score"
        assert "mean" in report["gene_sets"]["high"]
        uns = adata.uns["ucell_params"]
        assert uns["algorithm"] == "ucell"
        assert uns["algorithm_version"] == "UCell 2.16.0"
        assert uns["max_rank"] == 6
        assert uns["effective_max_rank"] == 6
        assert uns["ties_method"] == "average"
        assert uns["missing_genes"] == "impute"

    def test_effective_max_rank_caps_at_gene_count(self, tmp_path: Path) -> None:
        from ucell_fixtures import fixture_by_name

        fixture = fixture_by_name("multi_cell_default_params")
        _, report = run_ucell(fixture, tmp_path)
        assert report["max_rank"] == 1500
        assert report["effective_max_rank"] == 10  # 10 genes in the H5AD

    def test_present_missing_summary_semantics(self, tmp_path: Path) -> None:
        from ucell_fixtures import fixture_by_name

        fixture = fixture_by_name("missing_impute")
        _, report = run_ucell(fixture, tmp_path)
        summary = report["gene_sets"]["with_absent"]
        assert summary["genes_requested"] == 4
        assert summary["genes_present"] == 2
        assert summary["missing_genes"] == ["g11", "g12"]
        assert report["gene_sets"]["all_present"]["missing_genes"] == []


class TestContractErrors:
    def _run_with(self, tmp_path: Path, gene_sets, extra=None):
        from ucell_fixtures import fixture_by_name

        fixture = fixture_by_name("multi_cell_default_params")
        input_path = tmp_path / "input.h5ad"
        build_adata(fixture).write_h5ad(input_path, compression="lzf")
        os.environ["AUTONOMICS_INPUT0"] = str(input_path)
        os.environ["AUTONOMICS_OUTPUT0"] = str(tmp_path / "out.h5ad")
        os.environ["AUTONOMICS_OUTPUT1"] = str(tmp_path / "report.json")
        params = {"gene_sets": gene_sets}
        params.update(extra or {})
        try:
            workflow.ucell_score(params)
        finally:
            for name in ("AUTONOMICS_INPUT0", "AUTONOMICS_OUTPUT0", "AUTONOMICS_OUTPUT1"):
                os.environ.pop(name, None)

    def test_signature_longer_than_max_rank_rejected(self, tmp_path: Path) -> None:
        from ucell_fixtures import fixture_by_name

        fixture = fixture_by_name("multi_cell_default_params")
        with pytest.raises(workflow.ContractError, match="max_rank"):
            self._run_with(tmp_path, {"big": list(fixture.genes)}, {"max_rank": 5})

    def test_no_present_genes_rejected(self, tmp_path: Path) -> None:
        with pytest.raises(workflow.ContractError, match="no genes present"):
            self._run_with(tmp_path, {"ghost": ["zz1", "zz2"]})

    def test_invalid_ties_method_rejected(self, tmp_path: Path) -> None:
        with pytest.raises(workflow.ContractError, match="ties_method"):
            self._run_with(tmp_path, {"sig_a": ["g01", "g02"]}, {"ties_method": "ordinal"})

    def test_invalid_missing_mode_rejected(self, tmp_path: Path) -> None:
        with pytest.raises(workflow.ContractError, match="missing_genes"):
            self._run_with(tmp_path, {"sig_a": ["g01", "g02"]}, {"missing_genes": "zero"})

    def test_nonpositive_max_rank_rejected(self, tmp_path: Path) -> None:
        with pytest.raises(workflow.ContractError, match="max_rank"):
            self._run_with(tmp_path, {"sig_a": ["g01", "g02"]}, {"max_rank": 0})
