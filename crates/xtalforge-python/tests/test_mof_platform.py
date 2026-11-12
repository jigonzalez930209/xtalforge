"""Smoke tests for the Phase 10 Python surface."""

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


def test_mof_check_and_pores() -> None:
    report = xtalforge.mof_check(CIF)
    assert isinstance(report["findings"], list)
    pores = xtalforge.mof_pores(CIF, spacing=1.5)
    assert pores["descriptors"]["lcd"] >= 0.0
    assert 0 <= pores["descriptors"]["channel_dimension"] <= 3


def test_mof_deconstruction_and_identity() -> None:
    deconstruction = xtalforge.mof_deconstruct(CIF)
    assert isinstance(deconstruction["units"], list)
    identity = xtalforge.mof_mofid(CIF)
    assert identity["mofid"].startswith("MOFid:")
    assert len(identity["mofkey"]) == 64


def test_mof_adsorption_exports() -> None:
    guests = xtalforge.mof_guests()
    assert {guest["name"] for guest in guests} >= {"CO2", "N2", "CH4", "H2", "H2O"}
    widom = xtalforge.mof_widom(CIF, samples=20)
    assert widom["henry"] >= 0.0
    isotherm = xtalforge.mof_gcmc(CIF, pressures=[1.0], steps=10)
    assert len(isotherm["points"]) == 1
    raspa = xtalforge.mof_raspa(CIF)
    assert "MonteCarlo" in raspa["simulation_input"]


def test_mof_modification_exports() -> None:
    defects = xtalforge.mof_defects(CIF, linker_fraction=0.0)
    assert defects["removed_linkers"] == 0
    layers = xtalforge.mof_stack_layers(CIF, layers=2)
    assert len(layers["sites"]) == 16
    mixed = xtalforge.mof_mixed(CIF, [("Si", "Ge")], ordered=False, fraction=1.0)
    assert all(site["element"] == "Ge" for site in mixed["sites"])
