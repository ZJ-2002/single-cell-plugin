#!/usr/bin/env python3
"""Materialise the deterministic UCell fixtures to disk.

Writes, per fixture: an H5AD (what the plugin reads), a CSV of the same
matrix (what the R golden generator reads; %.17g keeps float64 exact),
and a spec JSON describing gene sets and UCell parameters. Run from the
plugin root:

    python tests/make_ucell_fixtures.py [output_dir]

Default output dir: tests/golden_inputs (gitignored; regenerate on
demand -- ucell_fixtures.py is fully deterministic).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ucell_fixtures import all_fixtures


def main() -> int:
    out_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent / "golden_inputs"
    out_dir.mkdir(parents=True, exist_ok=True)
    spec = []
    for fixture in all_fixtures():
        matrix = fixture.as_array()
        cells = fixture.cells
        h5ad_path = out_dir / f"{fixture.name}.h5ad"
        import anndata as ad

        ad.AnnData(
            X=sp.csr_matrix(matrix),
            obs=pd.DataFrame(index=pd.Index(cells, name="cell_id")),
            var=pd.DataFrame(index=pd.Index(fixture.genes, name="gene_id")),
        ).write_h5ad(h5ad_path, compression="lzf")
        csv_path = out_dir / f"{fixture.name}.csv"
        frame = pd.DataFrame(matrix, index=pd.Index(cells, name="cell_id"),
                             columns=fixture.genes)
        # %.17g round-trips float64 exactly, so R scores the identical matrix.
        frame.to_csv(csv_path, float_format="%.17g")
        spec.append({
            "name": fixture.name,
            "csv": csv_path.name,
            "cells": cells,
            "genes": fixture.genes,
            "gene_sets": fixture.gene_sets,
            "params": fixture.params,
        })
    spec_path = out_dir / "fixtures.json"
    spec_path.write_text(json.dumps(spec, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {len(spec)} fixtures to {out_dir} (spec: {spec_path})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
