"""Phonopy and spglib-compatible outputs for the Python adapters (Task 13.1.3)."""

from __future__ import annotations

import json

from ._interop_json import interop as _interop
from ._xtalforge import Structure, detect_symmetry, element_info


def structure_to_phonopy_atoms(structure_json):
    """Builds a phonopy ``PhonopyAtoms`` object from a P1 structure document."""
    from phonopy.structure.atoms import PhonopyAtoms

    document = json.loads(structure_json)["structure"]
    sites = document["sites"]
    arguments = {
        "symbols": [site["element"] for site in sites],
        "cell": document["lattice"],
        "scaled_positions": [site["frac"] for site in sites],
    }
    moments = _interop(document).get("magnetic_moments")
    if moments is not None and len(moments) == len(sites) and any(moments):
        arguments["magnetic_moments"] = moments
    return PhonopyAtoms(**arguments)


def get_symmetry_dataset(structure, symprec=1e-5):
    """Returns a spglib-shaped symmetry dataset for a structure or its text.

    The ``number``, ``hall_number``, ``international``, ``hall``, ``choice``,
    ``transformation_matrix``, ``origin_shift``, ``rotations``,
    ``translations``, ``wyckoffs``, ``site_symmetry_symbols``,
    ``equivalent_atoms``, ``crystallographic_orbits``, ``primitive_lattice``,
    ``mapping_to_primitive``, ``std_lattice``, ``std_types``,
    ``std_positions``, ``std_rotation_matrix`` and ``pointgroup`` keys follow
    spglib. ``std_mapping_to_primitive`` is not tracked by the core dataset
    (the primitive standardized cell has its own origin) and is always empty.
    """
    text = structure.to_json() if isinstance(structure, Structure) else structure
    data = detect_symmetry(text, symprec=symprec)
    operations = data.get("operations") or []
    return {
        "number": data.get("number"),
        "hall_number": data.get("hall_number"),
        "international": str(data.get("hm_symbol", "")).replace(" ", ""),
        "hall": data.get("hall_symbol", ""),
        "choice": data.get("setting", ""),
        "transformation_matrix": data.get("std_linear"),
        "origin_shift": data.get("std_origin_shift"),
        "rotations": [operation["rotation"] for operation in operations],
        "translations": [operation["translation"] for operation in operations],
        "wyckoffs": data.get("wyckoffs"),
        "site_symmetry_symbols": data.get("site_symmetry_symbols"),
        "crystallographic_orbits": data.get("orbits"),
        "equivalent_atoms": data.get("orbits"),
        "primitive_lattice": data.get("prim_std_lattice"),
        "mapping_to_primitive": data.get("mapping_to_primitive"),
        "std_lattice": data.get("std_lattice"),
        "std_types": [
            int(element_info(element)["z"]) for element in data.get("std_elements", [])
        ],
        "std_positions": data.get("std_positions"),
        "std_rotation_matrix": data.get("std_rotation_matrix"),
        "std_mapping_to_primitive": [],
        "pointgroup": data.get("point_group", ""),
    }
