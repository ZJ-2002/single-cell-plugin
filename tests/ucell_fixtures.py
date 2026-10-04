"""Deterministic UCell fixtures shared by the golden generator and pytest.

Each fixture is pure data: an expression matrix (cells x genes), gene and
cell names, gene-set definitions and the UCell parameters to run with.
`make_ucell_fixtures.py` serialises them to H5AD (python input) and CSV
(R input); `generate_ucell_golden.R` scores the CSVs with the official
UCell R package. Expected values live in `ucell_golden.json` and are
produced only by the R package, never by the plugin implementation.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np


@dataclass(frozen=True)
class UCellFixture:
    name: str
    matrix: list[list[float]]  # cells x genes, row = cell
    genes: list[str]
    gene_sets: dict[str, list[str]]
    params: dict[str, Any] = field(default_factory=dict)  # max_rank / ties_method / missing_genes

    @property
    def cells(self) -> list[str]:
        return [f"cell_{index:04d}" for index in range(len(self.matrix))]

    def as_array(self) -> np.ndarray:
        return np.asarray(self.matrix, dtype=np.float64)


def _gene_names(count: int, prefix: str = "g") -> list[str]:
    return [f"{prefix}{index:02d}" for index in range(1, count + 1)]


def _fixture_definitions() -> list[UCellFixture]:
    fixtures: list[UCellFixture] = []

    # Audit counterexample (workorder WO-5 / audit F01): 1 cell, 10 genes,
    # signature = the two highest-expressed genes. Official UCell with
    # maxRank=10 scores this 1.0; the pre-fix rational heuristic gave 0.7333.
    fixtures.append(
        UCellFixture(
            name="audit_counterexample",
            matrix=[[8.0, 1.0, 0.0, 3.0, 0.0, 2.0, 9.0, 0.0, 1.0, 0.0]],
            genes=_gene_names(10),
            gene_sets={"top2": ["g01", "g07"]},
            params={"max_rank": 10},
        )
    )

    # Distinct expression values, two signatures, maxRank left at the UCell
    # default (1500) so the gene-count cap (10 genes -> maxRank 10) applies.
    fixtures.append(
        UCellFixture(
            name="multi_cell_default_params",
            matrix=[
                [7.0, 0.0, 1.0, 4.0, 2.0, 0.0, 5.0, 3.0, 6.0, 8.0],
                [0.0, 6.0, 2.0, 0.0, 5.0, 1.0, 8.0, 3.0, 4.0, 7.0],
                [3.0, 3.0, 0.0, 1.0, 0.0, 2.0, 9.0, 8.0, 0.0, 1.0],
                [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0],
                [10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            ],
            genes=_gene_names(10),
            gene_sets={"sig_a": ["g03", "g08", "g10"], "sig_b": ["g01", "g05"]},
            params={},
        )
    )

    # Heavy ties inside and outside the signatures; scored under every
    # supported ties method (UCell default is "average").
    ties_matrix = [
        [5.0, 5.0, 5.0, 1.0, 0.0, 0.0, 2.0, 2.0, 4.0],
        [0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 3.0, 3.0, 3.0],
        [9.0, 8.0, 7.0, 7.0, 6.0, 6.0, 0.0, 0.0, 0.0],
        [2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0],
    ]
    for ties_method in ("average", "min", "max", "dense"):
        fixtures.append(
            UCellFixture(
                name=f"ties_{ties_method}",
                matrix=ties_matrix,
                genes=_gene_names(9),
                gene_sets={"tied": ["g01", "g02", "g07"], "tied_zero": ["g05", "g06"]},
                params={"max_rank": 9, "ties_method": ties_method},
            )
        )

    # Rank truncation: 20 genes but maxRank 6, so any signature gene ranked
    # 6th or lower contributes the clamped maxRank value.
    truncated = [[float((gene * 7 + cell * 3) % 21) for gene in range(20)] for cell in range(5)]
    fixtures.append(
        UCellFixture(
            name="truncated_ranks",
            matrix=truncated,
            genes=_gene_names(20),
            gene_sets={
                "high": ["g06", "g13", "g20"],  # top-ranked in most cells
                "low": ["g01", "g02", "g03"],  # mostly beyond maxRank
            },
            params={"max_rank": 6},
        )
    )

    # Partially absent signature genes: g11/g12 are not in the matrix.
    missing_matrix = [
        [4.0, 2.0, 0.0, 1.0, 3.0, 0.0, 5.0, 2.0, 0.0, 1.0],
        [0.0, 0.0, 1.0, 0.0, 0.0, 2.0, 0.0, 1.0, 3.0, 0.0],
        [9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0],
    ]
    for missing_genes in ("impute", "skip"):
        fixtures.append(
            UCellFixture(
                name=f"missing_{missing_genes}",
                matrix=missing_matrix,
                genes=_gene_names(10),
                gene_sets={
                    "with_absent": ["g01", "g02", "g11", "g12"],
                    "all_present": ["g03", "g09"],
                },
                params={"max_rank": 10, "missing_genes": missing_genes},
            )
        )

    # Several signatures, overlapping genes, default parameters.
    multi_matrix = [
        [float((gene * 5 + cell * 11) % 13) for gene in range(12)] for cell in range(4)
    ]
    fixtures.append(
        UCellFixture(
            name="multi_signature",
            matrix=multi_matrix,
            genes=_gene_names(12),
            gene_sets={
                "alpha": ["g01", "g04", "g09"],
                "beta": ["g02", "g07", "g13"],  # g13 absent
                "gamma": ["g01", "g07", "g09", "g12"],  # overlaps alpha/beta
            },
            params={"max_rank": 8},
        )
    )

    # Crosses the plugin's 1000-cell chunking boundary twice (2100 cells).
    rng = np.random.default_rng(20261005)
    big = rng.poisson(lam=2.0, size=(2100, 24)).astype(np.float64)
    fixtures.append(
        UCellFixture(
            name="chunked_large",
            matrix=big.tolist(),
            genes=_gene_names(24),
            gene_sets={"poisson_sig": ["g02", "g09", "g17", "g24"]},
            params={"max_rank": 24},
        )
    )

    return fixtures


def all_fixtures() -> list[UCellFixture]:
    return _fixture_definitions()


def fixture_by_name(name: str) -> UCellFixture:
    for fixture in _fixture_definitions():
        if fixture.name == name:
            return fixture
    raise KeyError(f"unknown UCell fixture `{name}`")
