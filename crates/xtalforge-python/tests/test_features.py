"""Smoke tests for the Phase 9 Python surface."""

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


def test_matcher_and_miller_indices() -> None:
    result = xtalforge.match_structures(CIF, CIF)
    assert result["matched"] is True
    assert result["rms"] < 1e-6
    indices = xtalforge.unique_miller_indices(CIF, max_index=1)
    assert any(
        tuple(index["hkl"]) == (1, 0, 0) and index["multiplicity"] == 6
        for index in indices
    )


def test_deform_and_rattle_are_deterministic() -> None:
    scaled = xtalforge.deform(CIF, volume=200.0)
    assert abs(scaled["volume"] - 200.0) < 1e-6
    first = xtalforge.rattle(CIF, amplitude=0.01, seed=7)
    second = xtalforge.rattle(CIF, amplitude=0.01, seed=7)
    assert first["max_displacement"] == second["max_displacement"]


def test_local_environment_and_rdf() -> None:
    report = xtalforge.voronoi(CIF, cutoff=4.0)
    assert len(report["cells"]) == 8
    assert abs(report["total_volume"] - 160.1862572980187) < 0.1
    rdf = xtalforge.radial_distribution(CIF, rmax=4.0, bin=0.05)
    assert len(rdf["r"]) == 80
    assert len(rdf["partial"]) >= 1


def test_eos_and_hull() -> None:
    fit = xtalforge.eos_fit(
        [(18.0, 1.0), (19.0, 0.4), (20.0, 0.0), (21.0, 0.5), (22.0, 1.6)]
    )
    assert 18.0 < fit["v0"] < 22.0
    hull = xtalforge.convex_hull(
        [
            ("A", {"A": 1.0}, 0.0),
            ("B", {"B": 1.0}, 0.0),
            ("AB", {"A": 1.0, "B": 1.0}, -1.0),
        ]
    )
    assert hull["points"][2]["stable"] is True


def test_edits_and_defects() -> None:
    edited = xtalforge.edit_structure(
        CIF, species=[("Si", "Ge")], fraction=0.5, seed=3
    )
    assert sorted(edited["formula"].split()) in (["Ge4", "Si4"], ["Si4", "Ge4"])
    defects = xtalforge.vacancies(CIF, min_image=8.0)
    assert len(defects) == 1
    assert len(defects[0]["sites"]) == 56
