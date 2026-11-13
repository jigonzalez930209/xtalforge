"""pymatgen SpacegroupAnalyzer parity (Task 13.1.2)."""

import pytest

import xtalforge

SILICON = [[5.430941, 0.0, 0.0], [0.0, 5.430941, 0.0], [0.0, 0.0, 5.430941]]


def test_pymatgen_spacegroup_number_matches():
    pytest.importorskip("pymatgen.core")
    from pymatgen.core import Structure as PmgStructure
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    source = PmgStructure(SILICON, ["Na", "K"], [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
    structure = xtalforge.Structure.from_pymatgen(source)
    expected = SpacegroupAnalyzer(source).get_space_group_number()
    assert xtalforge.get_symmetry_dataset(structure)["number"] == expected
