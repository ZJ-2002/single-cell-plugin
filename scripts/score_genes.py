"""Original-kind, bounded official Scanpy score_genes wrapper."""

import json
import os
from pathlib import Path
import sys


def unique_object(pairs):
    result = {}
    for name, value in pairs:
        if name in result:
            raise ValueError("Duplicate JSON key: " + name)
        result[name] = value
    return result


def parse_sets(text):
    sets = json.loads(text, object_pairs_hook=unique_object,
                      parse_constant=lambda value: (_ for _ in ()).throw(ValueError("Nonstandard JSON constant: " + value)))
    if not isinstance(sets, dict) or not sets:
        raise ValueError("gene_sets must be a nonempty JSON object")
    for name, genes in sets.items():
        if not isinstance(name, str) or not name.strip():
            raise ValueError("Gene-set names must be nonempty strings")
        if not isinstance(genes, list) or not genes or any(not isinstance(gene, str) or not gene.strip() for gene in genes):
            raise ValueError("Each gene set must contain nonempty gene strings")
    return sets


def integer(text, name, minimum, maximum=None):
    try:
        value = int(text)
    except (ValueError, TypeError) as error:
        raise ValueError(name + " must be an integer") from error
    if str(value) != text or value < minimum or (maximum is not None and value > maximum):
        raise ValueError(name + " outside the frozen integer domain")
    return value


def summarize(scores):
    import numpy as np

    finite = np.isfinite(scores)
    count = int(finite.sum())
    result = {"finite_cells": count, "nonfinite_cells": int(scores.size - count),
              "status": "FINITE" if count == scores.size else "NONFINITE_OFFICIAL_SCORES"}
    if count != scores.size:
        result.update(dict.fromkeys(("mean", "min", "p25", "median", "p75", "max")))
    else:
        values = np.percentile(scores, [0, 25, 50, 75, 100])
        result.update(dict(zip(("min", "p25", "median", "p75", "max"), map(float, values))))
        result["mean"] = float(np.mean(scores))
    return result


def main():
    import anndata as ad
    import numpy as np
    import pandas as pd
    import scanpy as sc

    if sc.__version__ != "1.11.3" or ad.__version__ != "0.11.4":
        raise ValueError("Frozen official package version differs")
    sets = parse_sets(os.environ["SC_P_JSON_GENE_SETS"])
    ctrl_size = integer(os.environ["SC_P_INT_CTRL_SIZE"], "ctrl_size", 1)
    seed = integer(os.environ["SC_P_INT_RANDOM_STATE"], "random_state", 0, 4294967295)
    data = ad.read_h5ad(os.environ["AUTONOMICS_INPUT0"])
    if data.n_obs == 0 or data.n_vars == 0:
        raise ValueError("Empty AnnData axis")
    if not data.obs_names.is_unique or not data.var_names.is_unique:
        raise ValueError("Duplicate AnnData axis")
    if not data.obs.columns.is_unique:
        raise ValueError("Duplicate obs column")
    use_raw = data.raw is not None
    axis = data.raw.var_names if use_raw else data.var_names
    if not axis.is_unique:
        raise ValueError("Duplicate raw gene axis")
    requests = {}
    for name, genes in sets.items():
        if name in data.obs.columns:
            raise ValueError("Score name collides with obs: " + name)
        present = pd.Index(genes).intersection(axis)
        if len(present) < 2:
            raise ValueError("Each gene set needs at least two distinct present genes: " + name)
        requests[name] = {"genes_requested": len(genes), "genes_present": len(present),
                          "duplicate_members": len(genes) - len(set(genes)),
                          "missing_genes": sorted(set(genes) - set(axis))}
    for name, genes in sets.items():
        sc.tl.score_genes(data, genes, score_name=name, ctrl_size=ctrl_size, random_state=seed,
                          ctrl_as_ref=True, n_bins=25, gene_pool=None, copy=False, use_raw=use_raw, layer=None)
        requests[name].update(summarize(np.asarray(data.obs[name], dtype=np.float64)))
    report = {"schema_version": "1.0", "operation": "gene_set_score", "cells": data.n_obs,
              "ctrl_size": ctrl_size, "random_state": seed, "gene_sets": requests,
              "matrix_axis": "raw" if use_raw else "X", "selected_axis_genes": len(axis),
              "method": "scanpy.tl.score_genes", "official_version": "1.11.3",
              "fixed_parameters": {"ctrl_as_ref": True, "n_bins": 25, "gene_pool": None, "layer": None},
              "summary_policy": "ALL_CELLS_OR_NULL_NO_MISSING_CELL_DELETION",
              "scientific_signoff_granted": False}
    content = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ad.settings.allow_write_nullable_strings = True
    data.write_h5ad(os.environ["AUTONOMICS_OUTPUT0"], compression="lzf")
    Path(os.environ["AUTONOMICS_OUTPUT1"]).write_text(content)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print("gene_set_score: " + str(error), file=sys.stderr)
        sys.exit(2)
