import csv
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from raccoon.utils import plotly_baltic


class FakeNode:
    def __init__(self, name, x, y, branch_type="leaf", children=None):
        self.name = name
        self.x = x
        self.y = y
        self.branchType = branch_type
        self.children = children or []
        self.traits = {}

    def is_leaflike(self):
        return self.branchType == "leaf"

    def is_node(self):
        return self.branchType == "internal"


class FakeTree:
    def __init__(self):
        self.root = FakeNode("root", x=0.0, y=1.5, branch_type="internal")
        self.tip_a = FakeNode("tipA", x=1.0, y=1.0, branch_type="leaf")
        self.tip_b = FakeNode("tipB", x=1.1, y=2.0, branch_type="leaf")
        self.root.children = [self.tip_a, self.tip_b]
        self.Objects = [self.root, self.tip_a, self.tip_b]

    def drawTree(self):
        return None

    def getInternal(self):
        return [self.root]


def _write_branch_snps_csv(path: Path) -> None:
    rows = [
        {"parent": "root", "child": "tipA", "site": "10", "snp": "G->A", "dimer": "GA"},
        {"parent": "root", "child": "tipB", "site": "20", "snp": "C->T", "dimer": "TC"},
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["parent", "child", "site", "snp", "dimer"])
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


@pytest.fixture(autouse=True)
def patch_tree(monkeypatch):
    monkeypatch.setattr(plotly_baltic, "load_tree", lambda *_args, **_kwargs: FakeTree())
    monkeypatch.setattr(plotly_baltic, "ensure_node_label", lambda node: node.name)
    yield


def test_build_tree_plot_includes_branch_mutations_and_labels(tmp_path):
    csv_path = tmp_path / "branch_snps.csv"
    _write_branch_snps_csv(csv_path)

    html = plotly_baltic.build_tree_plot("tree", branch_snps_path=str(csv_path))

    assert "Branch: root_tipA" in html
    assert ("10: G->A (GA)" in html) or ("10: G-\\u003eA (GA)" in html)
    assert "tipA" in html
    assert "Reset view" in html
    assert "Expand tree" in html


def test_build_tree_plot_without_branch_snps(tmp_path):
    html = plotly_baltic.build_tree_plot("tree", branch_snps_path=str(tmp_path / "missing.csv"))
    assert "Reset view" in html
    assert "Branch:" not in html
