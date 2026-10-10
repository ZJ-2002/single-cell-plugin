# single-cell plugin

Migrated from the legacy `single_cell_container` (10x preprocessor) and
`single_cell_h5ad` (ten H5AD workflow variants) wrappers in nodes-io. One
directory = one plugin family = one git-able unit, one image
(`single-cell-preprocessor`), eleven `[[nodes]]` entries.

## Layout

- `manifest.toml` — the eleven-node contract: params, ports, panel, image
  provenance
- `scripts/run_workflow.py` — shared runner shim for nine H5AD nodes;
  rebuilds the legacy `params.json` from `SC_P_*` env vars and execs the
  workflow runner baked into the image
- `Dockerfile`, `.dockerignore`, `requirements.txt`, `preprocess.py`,
  `workflow.py` — the image build tree, moved from
  `containers/single-cell-preprocessor/`
- `preprocess.py`, `workflow.py` — the two baked runner programs
  (`/opt/autonomics/preprocess.py`, `/opt/autonomics/workflow.py`)
- `test_preprocess.py` — unit tests for `preprocess.py` (`pytest`)
- `test_single_cell_preprocessor.sh`, `test_single_cell_workflow.sh` —
  real-Podman image baselines (build + inspect/ingest + QC/embed/subset)
- `IMAGE.md` — the runner-level contract inherited from the container
  README (env vars, file inputs/outputs, inspect vs ingest semantics)

## Nodes

| kind | operation | timeout | resources | outputs |
| --- | --- | --- | --- | --- |
| `single_cell_preprocessor` | inspect/ingest (MatrixMarket) | 3600 | 2 CPU / 8Gi / shm 1Gi | `preprocess_report.json`, `preprocessed.h5ad` |
| `h5ad_qc_filter` | `qc_filter` | 3600 | 2 CPU / 8Gi / shm 1Gi | `output.h5ad`, `report.json` |
| `h5ad_pca_neighbors_umap_leiden` | `pca_neighbors_umap_leiden` | 7200 | 4 CPU / 16Gi / shm 2Gi | `output.h5ad`, `report.json` |
| `h5ad_celltypist_annotate` | `celltypist_annotate` | 3600 | 2 CPU / 8Gi / shm 1Gi | `output.h5ad`, `report.json` |
| `h5ad_subset_by_obs` | `subset_by_obs` | 3600 | 2 CPU / 8Gi / shm 1Gi | `output.h5ad`, `report.json` |
| `sc_dense_ingest` | `dense_ingest` | 3600 | 2 CPU / 8Gi / shm 1Gi | `output.h5ad`, `report.json` |
| `h5ad_rank_genes_groups` | `rank_genes_groups` | 3600 | 2 CPU / 8Gi / shm 1Gi | `rank_genes_groups.parquet`, `report.json`, `output.h5ad` |
| `h5ad_cluster_mean_expression` | `cluster_mean_expression` | 3600 | 2 CPU / 8Gi / shm 1Gi | `cluster_mean_expression.parquet` |
| `gene_set_score` | `gene_set_score` | 3600 | 2 CPU / 8Gi / shm 1Gi | `output.h5ad`, `report.json` |
| `h5ad_marker_annotate` | `marker_annotate` | 3600 | 2 CPU / 8Gi / shm 1Gi | `output.h5ad`, `report.json` |
| `h5ad_ucell_score` | `ucell_score` | 3600 | 2 CPU / 8Gi / shm 1Gi | `output.h5ad`, `report.json` |

Every node runs isolated (no network) with a read-only root filesystem and
`PullPolicy::Missing`, exactly like the legacy wrappers.

## Execution shape

### Scoped scoring repair (2026-10-10)

`gene_set_score` now injects `scripts/score_genes.py` directly instead of
executing the older scorer baked into the image. This is necessary: editing
the repository's `workflow.py` alone cannot change a digest-pinned runtime.
The image, family panel, input/output ports and other ten node contracts
are unchanged. The scorer calls official Scanpy 1.11.3 with explicit
parameters, selects the raw gene axis when raw exists, does not clamp
`ctrl_size`, and preserves nonfinite official scores rather than silently
dropping cells from summary statistics. It rejects duplicate axes, output
column collisions and invalid signatures/seeds.

The supported wrapper requires two distinct present genes per signature;
this is not a claim to support every official Scanpy input or parameter.
Contract tests are `python -m unittest discover -s tests -p
'test_score_genes_contract.py' -v`. Run `tests/test_score_genes_container.py`
inside the manifest's pinned image for official-package comparisons.
These tests do not certify the complete family, the resident daemon, or
scientific interpretation. In particular, the legacy `h5ad_ucell_score`
implementation is not certified as official UCell. The family-level
CellTypist panel is still a deployment prerequisite for all nodes.

- `single_cell_preprocessor` is a **baked-runner** node: the compiled argv
  is `python /opt/autonomics/preprocess.py` (no injected script), and the
  five params travel through the same `AUTONOMICS_SINGLE_CELL_*` env vars
  the legacy wrapper set. The inspect-vs-filters cross-field rule is
  enforced inside `preprocess.py`, so no DSL-side validation is needed.
- Nine H5AD nodes are **script-injected** nodes with one shared shim;
  `gene_set_score` injects its dedicated official-package wrapper instead.
  The compiled command is `python` plus the staged shim (the legacy shape
  was `python` plus the staged `workflow.py` — byte-identical argv, a
  different program at argv[1]). The shim writes `params.json` to
  `/work/.autonomics/files/params.json` (the legacy path) from typed
  `SC_P_<KIND>_<name>` env vars (`INT`, `NUM`, `BOOL`, `LST`, `JSON`,
  `STR`; empty means unset and the workflow applies its own default), then
  runs `/opt/autonomics/workflow.py` from the image via `runpy`, so the
  Python implementation is the digest-pinned one, not a plugin copy.

## Install

```sh
export AUTONOMICS_PLUGIN_ROOT=/mnt/projects/node-plugins
cargo test -p container-plugin --test single_cell_migration
```

Declare the plugin in `~/.autonomics/plugins.toml` (pin `rev` to a commit
SHA once published):

```toml
[[plugin]]
name = "single-cell"
git = "git@github.com:auto-nomics/single-cell-plugin.git"
rev = "<commit sha>"
```

## Historical migration parity

The following describes the pre-repair baseline at `1d681eac`. The scoped
scoring repair deliberately changes `gene_set_score`'s injected script and
seed bounds, so the old eleven-node byte-equivalence statement is not an
acceptance claim for this revision. The other ten contracts remain exact.

The golden test (`container-plugin/tests/single_cell_migration.rs`) compares
each compiled `ContainerCommandSpec` against the legacy Rust wrapper output:
image, panels, outputs (path + format + order), resources, timeouts,
artifact prefixes, and argv are byte-exact; scripts are compared
semantically. Deliberate deltas, all preserving runtime behaviour:

- **Per-instance `artifact_prefix`, `timeout_secs`, `cpus`, `memory`,
  `pids_limit` params are gone.** The manifest DSL carries them as static
  node fields (set to the legacy defaults); the plugin schema no longer
  accepts per-instance overrides.
- **`gene_sets` / `marker_sets` are JSON strings, not objects.** The v0
  param DSL has scalar types only (the coloc flattening precedent). The
  shim `json.loads` them; malformed JSON fails at runtime instead of at
  spec build.
- **Enum params became strings** (`operation`, `metadata_separator`,
  `orientation`, `method`, `normalize`): the workflow runners enforce the
  same value sets the legacy Rust `validate` did, now at runtime.
- **`include_percent_expressed` defaults to `true`.** The legacy Rust spec
  declared `false`, but the params.json it shipped never materialized
  Rust-side defaults, and `workflow.py` applies `True` — the plugin keeps
  the effective value.
- **`h5ad_celltypist_annotate` takes its model from the catalog panel**
  (`wjixiang/catalog-celltypist-models-pan-immune` mounted at
  `/panels/celltypist_model`) instead of an optional second input port:
  the DSL cannot declare optional inputs or per-node data bundles. The
  shim injects the mounted model as the legacy `AUTONOMICS_INPUT1`.
  Dropped with it: `use_catalog_model` (always true now),
  `model_bundle` (fixed by the family panel), and `model_path` (arbitrary
  host paths cannot be staged by the manifest pipeline).
- **`sc_dense_ingest` input is now a required port** and the `path`
  param is gone: optional input ports and host-side file resolution have
  no manifest equivalent, so the matrix must arrive through the port.
- **Output ports gained labels and formats** (derived from the declared
  outputs); the legacy ports were bare Files.
- **Cross-field validation moved into the runners** (nonempty patterns,
  positive bounds, enum membership, nonempty gene sets): the same rules
  the legacy `validate` enforced, now enforced where the value is
  consumed.

The family panel (CellTypist model) is attached to all eleven nodes
because the DSL binds panels per family, not per node — the ldsc munge
precedent; ten of the nodes never read the mount.

## Build and baselines

```sh
test_single_cell_preprocessor.sh    # image build + inspect/ingest baseline
test_single_cell_workflow.sh        # image build + qc/embed/subset baseline
pytest test_preprocess.py           # preprocess.py unit tests
```

`SINGLE_CELL_PRODUCTION_BASELINE=1 test_single_cell_preprocessor.sh` adds
the deterministic 20,000-gene by 50,000-cell gzip MatrixMarket header-only
baseline. Set `BUILD_IMAGE=0` to reuse an existing local tag.
