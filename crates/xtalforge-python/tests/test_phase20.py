"""Phase 20 Python surfaces: generation, prototypes and search."""

from pathlib import Path

import pytest

xtalforge = pytest.importorskip("xtalforge")

ROOT = Path(__file__).resolve().parents[3]
NACL = (ROOT / "tests" / "corpus" / "simple" / "NaCl.cif").read_text()


def test_generate_crystal() -> None:
    crystal = xtalforge.generate_crystal(
        {"Na": 4, "Cl": 4}, space_group=225, seed=7
    )
    assert crystal["space_group"] == 225
    assert crystal["structure"]["name"] == "random-225"
    assert len(crystal["structure"]["sites"]) == 2


def test_generate_crystal_range() -> None:
    results = xtalforge.generate_crystal_range(
        {"Na": 4, "Cl": 4}, [216, 225], seed=5
    )
    assert len(results) == 2
    assert all(result["status"] == "ok" for result in results)


def test_prototype_library_and_identification() -> None:
    library = xtalforge.prototype_library()
    assert any(entry["label"] == "B1" for entry in library)
    built = xtalforge.build_prototype("B1")
    assert built["hall"]
    identification = xtalforge.identify_prototype(NACL)
    assert identification["label"] == "B1"


def test_substitution_candidates() -> None:
    report = xtalforge.substitution_candidates(NACL, radius_tolerance=0.4)
    assert report["neutral"]
    assert any(
        candidate["from"] == "Na" and candidate["to"] == "K"
        for candidate in report["candidates"]
    )


def test_random_search() -> None:
    report = xtalforge.random_search(
        {"Na": 4, "Cl": 4}, [225], structures=2, seed=3, relax=False
    )
    assert report["generated"] == 2
    assert report["hits"]
    assert report["calculator"] == "uff"


def test_search_workflow() -> None:
    document = xtalforge.search_workflow(
        "random", "seed.cif", {"n_structures": 1}, None, "out.json", "wf"
    )
    assert document["workflow_version"] == 1
    assert document["steps"][1]["verb"] == "search"


def test_generation_errors_are_typed() -> None:
    with pytest.raises(xtalforge.XtalforgeError):
        xtalforge.generate_crystal({"Na": 4}, space_group=999)
