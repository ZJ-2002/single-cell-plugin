import hashlib
import json
from pathlib import Path
import tomllib
import unittest


ROOT = Path(__file__).resolve().parents[1]


class ScoringContractTests(unittest.TestCase):
    def setUp(self):
        self.manifest = tomllib.loads((ROOT / "manifest.toml").read_text())
        self.node = next(node for node in self.manifest["nodes"] if node["kind"] == "gene_set_score")

    def test_other_ten_nodes_unchanged(self):
        others = {node["kind"]: node for node in self.manifest["nodes"] if node["kind"] != "gene_set_score"}
        self.assertEqual(len(others), 10)
        digest = hashlib.sha256(json.dumps(others, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
        self.assertEqual(digest, "c164978fc64814fb376a2047c5ff0f689a9ef4e83162054a3a788abb5ce0ac3c")

    def test_injected_script_and_ports(self):
        self.assertEqual(self.node["command"]["script_file"], "scripts/score_genes.py")
        self.assertEqual(len(self.node["ports"]["inputs"]), 1)
        self.assertEqual([port["path"] for port in self.node["ports"]["outputs"]], ["output.h5ad", "report.json"])
        self.assertEqual(set(self.node["command"]["env"]), {"SC_P_JSON_GENE_SETS", "SC_P_INT_CTRL_SIZE", "SC_P_INT_RANDOM_STATE"})

    def test_proved_script_bytes(self):
        digest = hashlib.sha256((ROOT / "scripts/score_genes.py").read_bytes()).hexdigest()
        self.assertEqual(digest, "9cad19728810b8e9c9933fed33f556870419a985a12f1485ce61d4f4da8c58e8")

    def test_seed_domain(self):
        self.assertEqual(self.node["params"]["random_state"]["min"], 0)
        self.assertEqual(self.node["params"]["random_state"]["max"], 4294967295)

    def test_retained_image_and_family_panel(self):
        self.assertEqual(self.manifest["image"]["reference"], "ghcr.io/auto-nomics/autonomics/single-cell-preprocessor@sha256:c34c26428d13804c2528bc734805c1ccfe909353604ac08da0ce9c68890e7e18")
        self.assertEqual(self.manifest["panels"][0]["bundle"], "wjixiang/catalog-celltypist-models-pan-immune")


if __name__ == "__main__":
    unittest.main()
