"""Phonopy and spglib-compatible outputs (Task 13.1.3)."""

from ._object_api import Structure

def structure_to_phonopy_atoms(structure_json: str) -> object: ...
def get_symmetry_dataset(
    structure: Structure | str, symprec: float = 1e-5
) -> dict[str, object]: ...
