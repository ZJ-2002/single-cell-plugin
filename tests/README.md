# UCell golden tests

`test_ucell_score.py` pins the `h5ad_ucell_score` node to the official
[UCell](https://github.com/carmonalab/UCell) R package algorithm.

## Layout

| File | Role |
|---|---|
| `ucell_fixtures.py` | Deterministic fixture definitions (data only; no scoring code) |
| `make_ucell_fixtures.py` | Serialises fixtures to `golden_inputs/` (H5AD + CSV + spec) |
| `generate_ucell_golden.R` | Scores the CSVs with official `ScoreSignatures_UCell`, writes `ucell_golden.json` |
| `ucell_golden.json` | Committed expected values; `_provenance` records the UCell + R versions used |
| `test_ucell_score.py` | Pytest: plugin output vs golden, `epsilon = 1e-14` |

## Independence guarantee

Expected values are produced **only** by the official R package. The plugin
implementation must never generate, adjust, or round these values — deriving
expectations from the code under test would be circular validation. The audit
counterexample (1 cell / 10 genes / top-2 signature, maxRank 10 -> exactly
1.0) is additionally hardcoded from the audit workorder.

## Regenerating (e.g. after editing fixtures)

```sh
python tests/make_ucell_fixtures.py          # deterministic; any python3 with anndata/pandas/scipy
R_LIBS=~/R/library Rscript tests/generate_ucell_golden.R
~/.cache/fix-venvs/ucell/bin/python -m pytest tests/ -q
```

`tests/golden_inputs/` is generated state (gitignored); the fixture source of
truth is `ucell_fixtures.py`, so regeneration is reproducible. CSVs are
written at `%.17g` so R reads the exact float64 matrix the plugin scores.

Current golden provenance: UCell 2.16.0, R 4.6.1 (`ucell_golden.json` ->
`_provenance`). A UCell package upgrade requires regenerating the golden file
and re-reviewing any value changes before committing.
