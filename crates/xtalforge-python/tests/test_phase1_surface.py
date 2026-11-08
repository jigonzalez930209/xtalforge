"""Integration tests for the Phase 1 structure and analysis surface."""

import json
import math
from pathlib import Path

import pytest
import xtalforge

ROOT = Path(__file__).parents[3]
SI = (ROOT / "tests" / "corpus" / "simple" / "Si.cif").read_text(encoding="utf-8")
MOF5 = (ROOT / "tests" / "corpus" / "prepared" / "MOF-5.cif").read_text(
    encoding="utf-8"
)

METHANE = """4
methane
C 0.000 0.000 0.000
H 0.629 0.629 0.629
H -0.629 -0.629 0.629
H -0.629 0.629 -0.629
"""

HEXAGON = "\n".join(
    ["6", "benzene"]
    + [
        f"C {1.4 * math.cos(math.tau * index / 6):.6f} "
        f"{1.4 * math.sin(math.tau * index / 6):.6f} 0.0"
        for index in range(6)
    ]
) + "\n"

POSCAR = """Si
1.0
5.430941 0.0 0.0
0.0 5.430941 0.0
0.0 0.0 5.430941
Si
1
Direct
0.0 0.0 0.0
"""


def test_read_cif_returns_the_cell_and_p1_cell() -> None:
    data = xtalforge.read_cif(SI)
    assert data["name"] == "9011998"
    assert data["formula"] == "Si"
    assert data["hall"] == "F 4d 2 3 -1d"
    assert data["n_sites"] == 1
    assert data["n_atoms"] == 8
    assert data["volume"] == pytest.approx(160.18625729801875, abs=1e-9)
    assert data["density"] > 0.0
    assert data["total_charge"] == 0.0
    assert data["is_p1"] is False
    assert data["lattice"]["a"] == pytest.approx(5.430941, abs=1e-12)
    assert data["p1"]["n_atoms"] == 8
    assert data["p1"]["elements"] == ["Si"] * 8
    assert len(data["p1"]["frac"]) == 8
    assert len(data["p1"]["cart"]) == 8
    assert len(data["p1"]["site_idx"]) == 8


def test_read_structure_dispatches_and_returns_the_same_summary() -> None:
    by_hint = xtalforge.read_structure(SI, "cif")
    by_extension = xtalforge.read_structure(SI, ".cif")
    direct = xtalforge.read_cif(SI)
    for key in (
        "name",
        "formula",
        "hall",
        "n_sites",
        "n_atoms",
        "volume",
        "density",
        "total_charge",
        "is_p1",
    ):
        assert by_hint[key] == direct[key]
        assert by_extension[key] == direct[key]


def test_structure_dictionary_embeds_the_versioned_json_document() -> None:
    data = xtalforge.read_cif(SI)
    document = json.loads(data["json"])
    assert document["schema_version"] == 1
    assert document["structure"]["name"] == "9011998"
    assert len(document["structure"]["sites"]) == 1


def test_supercell_expands_the_asymmetric_unit() -> None:
    data = xtalforge.supercell(SI, 2, 2, 2)
    assert data["repeat"] == [2, 2, 2]
    assert data["n_sites"] == 64
    assert data["n_atoms"] == 64
    assert data["is_p1"] is True
    assert data["volume"] == pytest.approx(8.0 * 160.18625729801875, rel=1e-12)


def test_slab_has_the_expected_surface_cell() -> None:
    data = xtalforge.slab(SI, 1, 0, 0, 1, 10.0)
    assert data["hkl"] == [1, 0, 0]
    assert data["layers"] == 1
    assert data["thickness"] == pytest.approx(5.430941, abs=1e-9)
    assert data["vacuum"] == 10.0
    assert data["n_sites"] == 8
    assert data["n_atoms"] == 8
    assert data["is_p1"] is True
    assert data["volume"] == pytest.approx(455.13745875282876, rel=1e-9)


def test_bonds_match_the_periodic_neighbors() -> None:
    bonds = xtalforge.bonds(SI)
    neighbors = xtalforge.neighbors(SI, 2.5)
    assert len(bonds) == 16
    assert len(neighbors) == 16
    assert all(bond["distance"] == pytest.approx(2.3516664362272315, abs=1e-9) for bond in bonds)
    assert all(bond["explicit"] is False for bond in bonds)
    assert all(bond["order"] == 1.0 for bond in bonds)


def test_rings_on_a_molecular_hexagon() -> None:
    found = xtalforge.rings(HEXAGON, 8)
    assert found == [[0, 1, 2, 3, 4, 5]]


def test_groups_on_mof5() -> None:
    data = xtalforge.read_cif(MOF5)
    assert data["name"] == "MOF-5_prepared"
    assert data["n_sites"] == 106
    assert data["n_atoms"] == 106
    assert data["volume"] > 0.0
    groups = xtalforge.groups(MOF5)
    assert any(group["kind"] == "carboxylate" for group in groups)
    assert all(
        0 <= atom < data["n_atoms"] for group in groups for atom in group["atoms"]
    )


def test_measurements_on_a_molecule() -> None:
    distance = xtalforge.distance(METHANE, 0, 1)
    assert distance == pytest.approx(math.sqrt(3.0) * 0.629, abs=1e-9)
    angle = xtalforge.angle(METHANE, 1, 0, 2)
    assert angle == pytest.approx(math.degrees(math.acos(-1.0 / 3.0)), abs=1e-9)
    torsion = xtalforge.torsion(METHANE, 0, 1, 2, 3)
    assert -180.0 <= torsion <= 180.0


def test_format_autodetection_accepts_xyz() -> None:
    data = xtalforge.read_structure(HEXAGON, "xyz")
    assert data["n_atoms"] == 6
    assert len(xtalforge.bonds(HEXAGON)) == 6


def test_read_structure_accepts_poscar() -> None:
    data = xtalforge.read_structure(POSCAR, "poscar")
    assert data["n_sites"] == 1
    assert data["n_atoms"] == 1
    assert data["volume"] == pytest.approx(160.18625729801875, abs=1e-9)
    assert data["is_p1"] is True


def test_invalid_format_raises_input_error() -> None:
    with pytest.raises(xtalforge.InputError) as caught:
        xtalforge.read_structure(SI, "unknown")
    assert caught.value.args[0] == "invalid_input"


def test_invalid_cif_raises_input_error() -> None:
    with pytest.raises(xtalforge.InputError) as caught:
        xtalforge.read_cif("data_broken\n")
    assert caught.value.args[0] == "invalid_input"


def test_bad_supercell_raises_input_error() -> None:
    with pytest.raises(xtalforge.InputError) as caught:
        xtalforge.supercell(SI, 0, 1, 1)
    assert caught.value.args[0] == "invalid_input"


def test_bad_slab_raises_input_error() -> None:
    with pytest.raises(xtalforge.InputError) as caught:
        xtalforge.slab(SI, 0, 0, 0, 1, 10.0)
    assert caught.value.args[0] == "invalid_input"
    with pytest.raises(xtalforge.InputError) as caught:
        xtalforge.slab(SI, 1, 0, 0, 1, -1.0)
    assert caught.value.args[0] == "invalid_input"


def test_bad_neighbor_and_ring_arguments_raise_input_error() -> None:
    with pytest.raises(xtalforge.InputError) as caught:
        xtalforge.neighbors(SI, 0.0)
    assert caught.value.args[0] == "invalid_input"
    with pytest.raises(xtalforge.InputError) as caught:
        xtalforge.rings(HEXAGON, 2)
    assert caught.value.args[0] == "invalid_input"


def test_bad_measurement_indices_raise_input_error() -> None:
    with pytest.raises(xtalforge.InputError) as caught:
        xtalforge.distance(METHANE, 0, 99)
    assert caught.value.args[0] == "invalid_input"
    with pytest.raises(xtalforge.InputError) as caught:
        xtalforge.angle(METHANE, -1, 0, 1)
    assert caught.value.args[0] == "invalid_input"


def test_python_classes() -> None:
    text = (Path(__file__).parents[3] / "tests" / "corpus" / "simple" / "Si.cif").read_text(
        encoding="utf-8"
    )
    lattice = xtalforge.Lattice(5.430941, 5.430941, 5.430941, 90, 90, 90)
    assert abs(lattice.volume - 5.430941**3) < 1e-9
    assert lattice.matrix.shape == (3, 3)
    assert lattice.to_cart([0.5, 0.5, 0.5]).shape == (3,)

    site = xtalforge.Site("C1", "C", 0.0, 0.0, 0.0).with_occupancy(0.5)
    assert site.occupancy == 0.5
    assert site.frac.shape == (3,)

    structure = xtalforge.Structure.from_cif(text)
    assert structure.formula == "Si"
    assert structure.n_sites == 1
    assert structure.p1().n_atoms == 8
    assert structure.supercell(2, 2, 2).p1().n_atoms == 64
    assert structure.to_format("xyz").startswith("8")
    assert "schema_version" in structure.to_json()


def test_detect_symmetry_surface() -> None:
    text = (Path(__file__).parents[3] / "tests" / "corpus" / "simple" / "Si.cif").read_text(
        encoding="utf-8"
    )
    dataset = xtalforge.detect_symmetry(text, symprec=0.001)
    assert dataset["number"] == 227
    assert dataset["hall_number"] == 525
    assert len(dataset["operations"]) == 192
    assert dataset["wyckoffs"] == ["b"] * 8
    assert dataset["orbits"] == [0] * 8
    assert len(dataset["std_lattice"]) == 3
