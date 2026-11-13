"""OBJ and glTF mesh exports (Task 14.4.4)."""

import json

import xtalforge

CIF = """
data_Si
_cell_length_a 5.43
_cell_length_b 5.43
_cell_length_c 5.43
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Si1 Si 0 0 0
"""


def test_mesh_obj_and_gltf() -> None:
    obj = xtalforge.mesh_obj(CIF)
    assert "\nv " in obj
    assert "\nvn " in obj
    assert "\nf " in obj
    gltf = json.loads(xtalforge.mesh_gltf(CIF))
    assert gltf["asset"]["version"] == "2.0"
    assert gltf["meshes"]
    assert gltf["accessors"]
