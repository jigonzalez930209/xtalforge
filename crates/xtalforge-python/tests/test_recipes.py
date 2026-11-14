"""Executable recipes on the Python surface (Task 17.1.5).

Every symbol named in the Python column of ``docs/guide/recipes.md`` must be
exported and callable, and the cheap non-MOF recipes run against the NaCl
corpus structure. The four MOF recipes are covered by the symbol check here
and by ``test_mof_platform.py``; they are not re-run because the recipe table
only promises the entry point.
"""

from __future__ import annotations

import math
from pathlib import Path

import pytest

xtalforge = pytest.importorskip("xtalforge")

NACL = (
    Path(__file__).resolve().parents[3]
    / "tests"
    / "corpus"
    / "simple"
    / "NaCl.cif"
).read_text(encoding="utf-8")

# The Python column of the recipe table, in table order.
RECIPE_SYMBOLS = (
    "read_structure",
    "write",
    "detect_symmetry",
    "to_p1",
    "supercell",
    "slab",
    "bonds",
    "distance",
    "pxrd",
    "ff_energy",
    "ff_optimize",
    "ff_charges",
    "match_structures",
    "group_structures",
    "deform",
    "rattle",
    "interpolate",
    "optimal_supercell",
    "voronoi",
    "coordination",
    "oxidation_states",
    "radial_distribution",
    "unique_miller_indices",
    "wulff",
    "adsorption_sites",
    "mof_check",
    "mof_mofid",
    "mof_pores",
    "mof_widom",
    "build_polymer",
)


def test_every_documented_symbol_exists() -> None:
    missing = [name for name in RECIPE_SYMBOLS if not hasattr(xtalforge, name)]
    assert not missing, f"missing recipe symbols: {missing}"
    not_callable = [
        name for name in RECIPE_SYMBOLS if not callable(getattr(xtalforge, name))
    ]
    assert not not_callable, f"recipe symbols are not callable: {not_callable}"


def test_cheap_recipes_run_on_nacl(tmp_path: Path) -> None:
    # The runtime signature needs the format hint; an empty hint sniffs the
    # CIF content, which is what the recipe promises.
    summary = xtalforge.read_structure(NACL, "")
    assert "Na" in summary["formula"] and "Cl" in summary["formula"]

    written = xtalforge.write(str(tmp_path / "a.vasp"), NACL)
    assert written["bytes"] > 0
    assert (tmp_path / "a.vasp").is_file()

    detected = xtalforge.detect_symmetry(NACL)
    assert int(detected["number"]) > 0

    assert xtalforge.to_p1(NACL)["n_atoms"] == 8
    assert xtalforge.supercell(NACL, 2, 2, 2)["n_atoms"] == 64
    assert xtalforge.slab(NACL, 1, 1, 1, 3, 15)["n_atoms"] > 0
    assert len(xtalforge.bonds(NACL)) > 0

    distance = xtalforge.distance(NACL, 0, 1)
    assert math.isfinite(distance) and distance > 0.0

    assert len(xtalforge.pxrd(NACL)["reflections"]) > 0
    assert math.isfinite(xtalforge.ff_energy(NACL))
    assert "Na" in xtalforge.ff_optimize(NACL)["formula"]

    charges = xtalforge.ff_charges(NACL, "eeq")
    assert len(charges["per_atom"]) == 8
    assert math.isfinite(charges["total"])

    assert xtalforge.match_structures(NACL, NACL)["matched"] is True
    assert xtalforge.group_structures([NACL, NACL])["n_groups"] == 1
    assert xtalforge.deform(NACL, voigt=[0.01, 0.0, 0.0, 0.0, 0.0, 0.0])[
        "n_atoms"
    ] == 8
    assert xtalforge.rattle(NACL, seed=7)["n_atoms"] == 8
    assert len(xtalforge.interpolate(NACL, NACL, 5)) == 5
    assert xtalforge.optimal_supercell(NACL, 8)["n_atoms"] > 8

    assert len(xtalforge.voronoi(NACL)["cells"]) > 0
    assert int(xtalforge.coordination(NACL, 0)["number"]) > 0
    assert len(xtalforge.oxidation_states(NACL)["sites"]) > 0
    assert len(xtalforge.radial_distribution(NACL)["total"]) > 0
    assert len(xtalforge.unique_miller_indices(NACL)) > 0
    assert len(xtalforge.wulff(NACL, [[1.0, 0.0, 0.0, 1.0]])["facets"]) > 0
    assert len(xtalforge.adsorption_sites(NACL)) > 0
    assert "C" in xtalforge.build_polymer("c1ccc(cc1)C=C", 10)["formula"]
