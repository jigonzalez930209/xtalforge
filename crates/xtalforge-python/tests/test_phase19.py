"""Phase 19 Python surfaces: magnetism, disorder, diffraction."""

from pathlib import Path

import pytest

xtalforge = pytest.importorskip("xtalforge")

ROOT = Path(__file__).resolve().parents[3]
NACL = (ROOT / "tests" / "corpus" / "simple" / "NaCl.cif").read_text()
MCIF = """
data_mn_test
_cell_length_a 5.0
_cell_length_b 4.0
_cell_length_c 8.0
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_space_group_name_H-M 'P 1'
_space_group_name_Hall 'P 1'
loop_
_space_group_symop_magn_operation_xyz
'x,y,z,+1'
'-x,-y,z,-1'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Mn1 Mn 0.0 0.0 0.0 1.0
Mn2 Mn 0.5 0.5 0.0 0.5
Fe2 Fe 0.5 0.5 0.0 0.5
O1 O 0.25 0.25 0.25 1.0
loop_
_atom_site_moment_label
_atom_site_moment_crystalaxis_x
_atom_site_moment_crystalaxis_y
_atom_site_moment_crystalaxis_z
Mn1 0.0 0.0 4.0
Mn2 0.0 0.0 -4.0
"""


def test_magnetic_orderings_and_disorder():
    orderings = xtalforge.magnetic_orderings(MCIF, supercell=1, max_orderings=4)
    assert orderings
    assert all("magmom" in ordering for ordering in orderings)
    disorder = xtalforge.disorder_groups(MCIF)
    assert disorder["groups"]
    derivatives = xtalforge.derivative_structures(MCIF, max_index=2, max_structures=2)
    assert derivatives


def test_diffraction_patterns():
    neutron = xtalforge.neutron_powder(NACL, wavelength=1.54, step=0.1)
    assert neutron["reflections"]
    assert neutron["points"]
    electron = xtalforge.electron_diffraction(NACL, zone=(0, 0, 1), size=64)
    assert electron["spots"]
    magnetic = xtalforge.magnetic_powder(MCIF, wavelength=1.54, step=0.2, allow_missing_moments=True)
    assert "reflections" in magnetic


def test_tensor_averages():
    assert xtalforge.raman_powder_average([[[1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]]])[0] > 0
    assert xtalforge.born_diagonal_average([[[2.0, 0, 0], [0, 2.0, 0], [0, 0, 2.0]]])[0] == 2.0
