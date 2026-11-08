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


def test_symmetry_workflow_surface() -> None:
    text = (Path(__file__).parents[3] / "tests" / "corpus" / "simple" / "Si.cif").read_text(
        encoding="utf-8"
    )
    expanded = xtalforge.to_p1(text)
    assert expanded["is_p1"] is True
    assert expanded["n_sites"] == 8

    asymmetric = xtalforge.asymmetrize(text, symprec=0.001)
    assert asymmetric["n_sites"] == 1
    assert asymmetric["hall"] == "F 4d 2 3 -1d"

    report = xtalforge.set_space_group(asymmetric["json"], "227", "reinterpret")
    assert report["set_report"]["n_output"] == 1
    assert report["hall"] == "F 4d 2 3 -1d"


def test_subgroups_surface() -> None:
    groups = xtalforge.subgroups(525)
    assert groups
    assert any(group["subgroup_number"] == 216 for group in groups)
    text = (Path(__file__).parents[3] / "tests" / "corpus" / "simple" / "Si.cif").read_text(
        encoding="utf-8"
    )
    lowered = xtalforge.lower_to_subgroup(text, 216, 0.5)
    assert lowered["hall"] == "F -4 2 3"


def test_standardize_reduce_wyckoff_and_zone() -> None:
    text = (Path(__file__).parents[3] / "tests" / "corpus" / "simple" / "Si.cif").read_text(
        encoding="utf-8"
    )
    standard = xtalforge.standardize(text, symprec=0.001)
    assert standard["n_sites"] == 1
    primitive = xtalforge.primitive_cell(text)
    assert primitive["p1"]["n_atoms"] == 2
    reduced = xtalforge.reduce_structure(text)
    assert reduced["determinant"] in (-1.0, 1.0)
    assignments = xtalforge.wyckoff_assignments(text, symprec=0.001)
    assert assignments
    assert assignments[0]["letter"] == "b"
    zone = xtalforge.brillouin_zone(text)
    assert len(zone["vertices"]) == 8
    ladder = xtalforge.tolerance_ladder(text)
    assert ladder
    diagnostics = xtalforge.diagnose(text, "227")
    assert diagnostics["max_displacement"] < 1e-6
    path = xtalforge.band_path(text)
    assert path["bravais"] == "cF"
    assert "L" in path["special_points"]
    changed = xtalforge.change_setting(text, 526)
    assert changed["hall"] == "-F 4vw 2vw 3"
    assert changed["p1"]["n_atoms"] == 8


def test_layer_group_surface() -> None:
    text = """data_layer
_cell_length_a 4
_cell_length_b 4
_cell_length_c 20
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
C1 C 0 0 0
C2 C 0.5 0.5 0
"""
    dataset = xtalforge.detect_layer_symmetry(text, 0.001)
    assert 1 <= dataset["number"] <= 80
    assert dataset["operations"]


def test_view_scene_and_stl() -> None:
    text = (Path(__file__).parents[3] / "tests" / "corpus" / "simple" / "Si.cif").read_text(
        encoding="utf-8"
    )
    scene = xtalforge.scene_summary(text, "space-fill", 1.0)
    assert scene["atoms"] == 8
    assert scene["cell_edges"] == 12
    stl = xtalforge.mesh_stl(text, "space-fill")
    assert stl.startswith("solid")
    assert "facet normal" in stl


def test_pxrd_surface() -> None:
    text = (Path(__file__).parents[3] / "tests" / "corpus" / "simple" / "Si.cif").read_text(
        encoding="utf-8"
    )
    pattern = xtalforge.pxrd(text, max_two_theta=60.0)
    assert pattern["reflections"]
    assert len(pattern["points"]) > 100
    assert all(b[0] <= b[1] for b in zip(
        [point[0] for point in pattern["points"]],
        [point[0] for point in pattern["points"][1:]],
    ))


def test_phase3_document_and_project_surface() -> None:
    text = (Path(__file__).parents[3] / "tests" / "corpus" / "simple" / "Si.cif").read_text(
        encoding="utf-8"
    )
    document = xtalforge.Document(text)
    assert document.can_undo is False
    document.move_sites([0], [0.01, 0.0, 0.0])
    assert document.can_undo is True
    assert document.dirty is True
    assert document.history == ["Move site"]
    document.undo()
    assert document.can_redo is True
    document.redo()
    assert document.can_redo is False

    archive = xtalforge.project_pack(text)
    assert archive[:2] == b"PK"
    payload = xtalforge.project_read(archive)
    assert payload["formula"] == "Si"

    schema = xtalforge.param_schema("symmetry")
    assert schema["properties"]["symprec"]["type"] == "number"


def test_workspace_surface(tmp_path: Path) -> None:
    import json

    root = tmp_path / "workspace"
    xtalforge.workspace_open(str(root))
    record = {
        "id": "run-1",
        "module": "demo",
        "input": "Si.cif",
        "status": "Done",
        "outputs": [],
        "data": {},
        "created": 1,
    }
    path = xtalforge.workspace_write_run(str(root), json.dumps(record))
    assert Path(path).exists()
    assert xtalforge.workspace_runs(str(root)) == ["run-1"]
    assert json.loads(xtalforge.workspace_read_run(str(root), "run-1"))["module"] == "demo"
    archive = xtalforge.workspace_archive(json.dumps([record]))
    assert xtalforge.workspace_read_archive(archive) == ["run-1"]


def test_force_field_surface() -> None:
    text = (Path(__file__).parents[3] / "tests" / "corpus" / "simple" / "Si.cif").read_text(
        encoding="utf-8"
    )
    energy = xtalforge.ff_energy(text)
    assert isinstance(energy, float)
    types = xtalforge.ff_types(text)
    assert types[0]["type"] == "Si3"
    vibrations = xtalforge.ff_vibrations(text)
    assert len(vibrations["frequencies"]) == 24
    thermo = xtalforge.ff_thermo(text, temperature=300.0)
    assert isinstance(thermo["gibbs"], float)
    optimized = xtalforge.ff_optimize(text, max_steps=2)
    assert isinstance(optimized["energy"], float)
    salt = (Path(__file__).parents[3] / "tests" / "corpus" / "simple" / "NaCl.cif").read_text(
        encoding="utf-8"
    )
    charges = xtalforge.ff_charges(salt, "qeq")
    assert abs(charges["total"]) < 1e-6


def test_analysis_and_hydrogens_surface() -> None:
    text = (Path(__file__).parents[3] / "tests" / "corpus" / "simple" / "Si.cif").read_text(
        encoding="utf-8"
    )
    metrics = xtalforge.porosity(text, spacing=0.5)
    assert metrics["cell_volume"] > 0.0
    assert metrics["surface_area"] >= 0.0
    names = xtalforge.rcsr_names()
    assert "pcu" in names

    cgd = """CRYSTAL
NAME pcu
GROUP P1
NODE A 6 0 0 0
NODE B 6 0.5 0 0
EDGE 0 0 0 0.5 0 0
END
"""
    invariants = xtalforge.net_invariants(cgd)
    assert invariants[0]["name"] == "pcu"
    assert invariants[0]["nodes"] == 2
    overlay = xtalforge.net_overlay(cgd)
    assert len(overlay["sites"]) == 2

    methane = "1\nmethane\nC 0 0 0\n"
    filled = xtalforge.add_hydrogens(methane)
    assert filled["n_sites"] == 5


def test_prepare_surface() -> None:
    text = (Path(__file__).parents[3] / "tests" / "corpus" / "simple" / "Si.cif").read_text(
        encoding="utf-8"
    )
    report = xtalforge.ff_prepare(
        text, minimize=False, charges=True, ewald=False, vibrations=False
    )
    assert report["energy"] is not None
    assert abs(report["total_charge"]) < 1e-6
    assert report["steps"][0]["name"]


def test_mof_surface() -> None:
    topologies = xtalforge.mof_topologies()
    assert any(item["name"] == "pcu" for item in topologies)
    blocks = xtalforge.mof_blocks()
    assert any(item["name"] == "Zn4O" for item in blocks)
    built = xtalforge.build_mof("pcu", "Zn4O", "benzene")
    assert built["n_sites"] > 0
    assert isinstance(built["warnings"], list)


def test_extra_builders_surface() -> None:
    c20 = xtalforge.build_c20()
    assert c20["n_sites"] == 20
    twist = xtalforge.twist(2, 2, 0.0, 2)
    assert twist["n_sites"] == 16
    ribbon = xtalforge.ribbon(2, 3, 12.0)
    assert ribbon["n_sites"] > 0
    ztc = xtalforge.build_ztc(samples=8, cell=12.0)
    assert ztc["n_sites"] > 0
    schwarz = xtalforge.build_schwarzite(samples=8, cell=12.0)
    assert schwarz["n_sites"] > 0
    polymer = xtalforge.build_polymer("C=C", 3)
    assert polymer["n_sites"] > 0
    assert polymer["steps"]
