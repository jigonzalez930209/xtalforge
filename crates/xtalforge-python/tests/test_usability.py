"""Tests for the Phase 11 usability surface."""

import os
import tempfile

import xtalforge

CIF = """data_Si
_cell_length_a 5.430941
_cell_length_b 5.430941
_cell_length_c 5.430941
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'F d -3 m'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Si1 Si 0 0 0
"""


def test_read_write_roundtrip_and_warnings() -> None:
    with tempfile.TemporaryDirectory() as directory:
        cif_path = os.path.join(directory, "si.cif")
        with open(cif_path, "w", encoding="utf-8") as handle:
            handle.write(CIF)
        structure = xtalforge.read(cif_path)
        assert structure["formula"] == "Si8"
        assert isinstance(structure["warnings"], list)
        poscar_path = os.path.join(directory, "si.vasp")
        written = xtalforge.write(poscar_path, CIF, "cif")
        assert written["bytes"] > 0
        rewritten = xtalforge.read(poscar_path)
        assert rewritten["formula"] == "Si8"


def test_molecule_and_sketch() -> None:
    molecule = xtalforge.molecule_from_smiles("c1ccccc1")
    assert molecule["formula"] == "C6"
    assert xtalforge.sketch_roundtrip("c1ccccc1")
    names = xtalforge.fragment_names()
    assert "Water" in names
    fragment = xtalforge.fragment_structure("Water")
    assert fragment["formula"] == "O"


def test_substitute_and_guest() -> None:
    replaced = xtalforge.substitute_structure(CIF, 0, "Ge")
    assert replaced["sites"][0]["element"] == "Ge"
    host = xtalforge.place_guest_structure(
        CIF,
        [("O", [0.0, 0.0, 0.0])],
        offset=(2.7, 2.7, 2.7),
        clearance=0.3,
    )
    assert any(site["element"] == "O" for site in host["sites"])


def test_error_category_and_threads() -> None:
    assert xtalforge.error_category("invalid_input") == "input"
    assert xtalforge.error_category("limit_exceeded") == "limit"
    assert xtalforge.set_threads(1) is True


def test_object_api_numpy_and_batch() -> None:
    structure = xtalforge.Structure(CIF)
    assert structure.formula == "Si8"
    assert len(structure.positions) == 8
    assert len(structure.frac_positions) == 8
    assert len(structure.cell) == 3
    assert len(structure.numbers) == 8
    rebuilt = xtalforge.Structure.from_arrays(
        structure.cell, [int(z) for z in structure.numbers], frac=structure.frac_positions
    )
    assert rebuilt.formula == "Si8"
    moved = structure.set_positions(structure.frac_positions)
    assert moved.n_sites == 8
    dataset = xtalforge.detect_symmetry(structure)
    assert isinstance(dataset, xtalforge.SymmetryDataset)
    assert dataset["number"] == 227
    supercell = xtalforge.supercell(structure, 2, 1, 1)
    assert supercell["formula"] == "Si16"
    batch = __import__("xtalforge.batch", fromlist=["detect_symmetry", "pxrd"])
    assert len(batch.detect_symmetry([CIF])) == 1
    assert len(batch.pxrd([CIF])) == 1
