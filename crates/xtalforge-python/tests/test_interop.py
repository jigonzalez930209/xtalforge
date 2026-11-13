"""Integration tests for the Phase 13.1 Python object adapters."""

import json
import subprocess
import sys
import types

import numpy as np
import pytest

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

SILICON = [[5.430941, 0.0, 0.0], [0.0, 5.430941, 0.0], [0.0, 0.0, 5.430941]]


def document(structure):
    return json.loads(structure.to_json())["structure"]


def interop(structure):
    return document(structure)["metadata"]["interop"]


class FakeAseAtoms:
    def __init__(self, symbols, positions, cell, pbc, arrays=None, constraints=None):
        self.symbols = symbols
        self.positions = positions
        self.cell = cell
        self.pbc = pbc
        self.arrays = arrays or {}
        self.constraints = constraints or []
        self.info = {"tag": "fake"}
        self.calc = None

    def get_chemical_symbols(self):
        return self.symbols

    def get_positions(self):
        return self.positions

    def get_cell(self):
        return self.cell

    def get_pbc(self):
        return self.pbc


class FakeFixAtoms:
    def __init__(self, indices):
        self.indices = indices


class FakePmgLattice:
    def __init__(self, matrix):
        self.matrix = matrix


class FakePmgSite:
    def __init__(self, species, frac_coords=None, coords=None, properties=None, label=None):
        self.species = species
        self.frac_coords = frac_coords
        self.coords = coords
        self.properties = properties or {}
        self.label = label
        self.occupancy = sum(species.values())


class FakePmgStructure:
    def __init__(self, lattice, sites, pbc=(True, True, True)):
        self.lattice = FakePmgLattice(lattice)
        self.sites = sites
        self.pbc = pbc


class FakePmgMolecule:
    def __init__(self, sites):
        self.sites = sites


class FakeRdkitAtom:
    def __init__(self, symbol, charge=0):
        self.symbol = symbol
        self.charge = charge

    def GetSymbol(self):
        return self.symbol

    def GetFormalCharge(self):
        return self.charge


class FakePoint:
    def __init__(self, x, y, z):
        self.x, self.y, self.z = x, y, z


class FakeRdkitBond:
    def __init__(self, first, second, order):
        self.first, self.second, self.order = first, second, order

    def GetBeginAtomIdx(self):
        return self.first

    def GetEndAtomIdx(self):
        return self.second

    def GetBondTypeAsDouble(self):
        return self.order


class FakeRdkitMol:
    def __init__(self, atoms, bonds, positions):
        self.atoms, self.bonds, self.positions = atoms, bonds, positions

    def GetNumAtoms(self):
        return len(self.atoms)

    def GetAtomWithIdx(self, index):
        return self.atoms[index]

    def GetConformer(self):
        return self

    def GetAtomPosition(self, index):
        return FakePoint(*self.positions[index])

    def GetBonds(self):
        return self.bonds


class FakePhonopyAtoms:
    def __init__(self, **kwargs):
        self.kwargs = kwargs


class FakeRxAtom:
    def __init__(self, symbol):
        self.symbol = symbol
        self.charge = 0
        self.aromatic = False

    def SetFormalCharge(self, charge):
        self.charge = charge

    def SetIsAromatic(self, value):
        self.aromatic = value


class FakeRxConformer:
    def __init__(self, count):
        self.points = [None] * count

    def SetAtomPosition(self, index, point):
        self.points[index] = point


class FakeRxMol:
    def __init__(self):
        self.atoms, self.bonds, self.conformer = [], [], None

    def AddAtom(self, atom):
        self.atoms.append(atom)
        return len(self.atoms) - 1

    def AddBond(self, first, second, order):
        self.bonds.append((first, second, order))

    def GetMol(self):
        return self

    def GetAtomWithIdx(self, index):
        return self.atoms[index]

    def AddConformer(self, conformer, assignId=True):
        self.conformer = conformer


class FakeRxStructure:
    def __init__(self, lattice, species, coords, **kwargs):
        self.lattice, self.species, self.coords, self.kwargs = (
            lattice,
            species,
            coords,
            kwargs,
        )
        self.site_properties = kwargs.get("site_properties")


class FakeRxMolecule:
    def __init__(self, species, coords, site_properties=None):
        self.species, self.coords = species, coords
        self.site_properties = site_properties


class FakeSpecies:
    def __init__(self, element, charge=0):
        self.element = element
        self.oxi_state = charge


def test_from_ase_fake_atoms_keeps_metadata():
    atoms = FakeAseAtoms(
        ["Si", "Si"],
        [[0.0, 0.0, 0.0], [2.715, 2.715, 2.715]],
        SILICON,
        [True, True, True],
        arrays={
            "tags": [1, 2],
            "initial_charges": [1.5, -1.5],
            "initial_magnetic_moments": [0.5, 0.0],
            "occupancy": [1.0, 0.5],
        },
        constraints=[FakeFixAtoms([1])],
    )
    structure = xtalforge.Structure.from_ase(atoms)
    assert structure.n_sites == 2
    assert [site.element for site in structure.sites] == ["Si", "Si"]
    data = interop(structure)
    assert data["tags"] == [1, 2]
    assert data["charges"] == [1.5, -1.5]
    assert data["magnetic_moments"] == [0.5, 0.0]
    assert data["occupancy"] == [1.0, 0.5]
    assert data["constraints"] == [{"class": "FakeFixAtoms", "indices": [1]}]
    assert data["info"] == {"tag": "fake"}
    assert document(structure)["sites"][1]["occupancy"] == 0.5


def test_ase_real_roundtrip():
    ase = pytest.importorskip("ase")
    constraints = pytest.importorskip("ase.constraints")
    atoms = ase.Atoms(
        "Si2",
        positions=[[0.0, 0.0, 0.0], [2.715, 2.715, 2.715]],
        cell=SILICON,
        pbc=True,
    )
    atoms.set_tags([1, 2])
    atoms.set_initial_charges([1.5, -1.5])
    atoms.set_initial_magnetic_moments([0.5, 0.0])
    atoms.set_constraint(constraints.FixAtoms(indices=[1]))
    from ase.calculators.singlepoint import SinglePointCalculator

    atoms.calc = SinglePointCalculator(atoms, energy=-1.25, forces=np.zeros((2, 3)))
    structure = xtalforge.Structure.from_ase(atoms)
    back = structure.to_ase()
    assert back.get_chemical_symbols() == ["Si", "Si"]
    np.testing.assert_allclose(back.positions, atoms.positions, atol=1e-10)
    np.testing.assert_allclose(back.cell, atoms.cell, atol=1e-10)
    assert list(back.pbc) == [True, True, True]
    assert list(back.get_tags()) == [1, 2]
    np.testing.assert_allclose(back.get_initial_charges(), [1.5, -1.5])
    np.testing.assert_allclose(back.get_initial_magnetic_moments(), [0.5, 0.0])
    assert list(back.constraints[0].index) == [1]
    assert back.get_potential_energy() == pytest.approx(-1.25)
    np.testing.assert_allclose(back.get_forces(), np.zeros((2, 3)))


def test_ase_molecule_roundtrip():
    ase = pytest.importorskip("ase")
    water = ase.Atoms(
        "OH2", positions=[[0.0, 0.0, 0.0], [0.957, 0.0, 0.0], [-0.24, 0.927, 0.0]]
    )
    structure = xtalforge.Structure.from_ase(water)
    assert structure.is_p1
    assert interop(structure)["molecule"] is True
    back = structure.to_ase()
    assert not back.pbc.any()
    np.testing.assert_allclose(back.positions, water.positions, atol=1e-10)


def test_ase_calculator():
    ase = pytest.importorskip("ase")
    atoms = ase.Atoms("H2", positions=[[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]], cell=[5, 5, 5])
    atoms.calc = xtalforge.AseCalculator()
    energy = atoms.get_potential_energy()
    forces = atoms.get_forces()
    assert energy == pytest.approx(xtalforge.Structure.from_ase(atoms).ff_energy())
    assert forces.shape == (2, 3)
    assert np.isfinite(forces).all()
    assert xtalforge.adapters.UffAseCalculator is xtalforge.adapters.AseCalculator


def test_from_pymatgen_fake_keeps_disorder():
    structure = FakePmgStructure(
        SILICON,
        [
            FakePmgSite({"Na": 0.5, "K": 0.5}, frac_coords=[0.0, 0.0, 0.0]),
            FakePmgSite(
                {"Si": 1.0},
                frac_coords=[0.5, 0.5, 0.5],
                properties={"xtalforge_tag": 3},
            ),
        ],
    )
    converted = xtalforge.Structure.from_pymatgen(structure)
    sites = document(converted)["sites"]
    assert [site["element"] for site in sites] == ["Na", "K", "Si"]
    assert [site["occupancy"] for site in sites] == [0.5, 0.5, 1.0]
    assert sites[0]["disorder_group"] == sites[1]["disorder_group"] == 1
    assert interop(converted)["tags"] == [0.0, 0.0, 3.0]


def test_from_pymatgen_fake_molecule():
    molecule = FakePmgMolecule(
        [
            FakePmgSite({"O": 1.0}, coords=[0.0, 0.0, 0.0]),
            FakePmgSite({"H": 1.0}, coords=[0.957, 0.0, 0.0]),
        ]
    )
    converted = xtalforge.Structure.from_pymatgen(molecule)
    assert converted.formula == "O H"
    assert interop(converted)["molecule"] is True


def test_to_pymatgen_fake_module(monkeypatch):
    core = types.ModuleType("pymatgen.core")
    core.Lattice = FakePmgLattice
    core.Species = FakeSpecies
    core.Structure = FakeRxStructure
    core.Molecule = FakeRxMolecule
    package = types.ModuleType("pymatgen")
    package.core = core
    monkeypatch.setitem(sys.modules, "pymatgen", package)
    monkeypatch.setitem(sys.modules, "pymatgen.core", core)
    structure = xtalforge.Structure.from_arrays(
        SILICON, [11, 19], frac=[[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
    )
    converted = structure.to_pymatgen()
    assert isinstance(converted, FakeRxStructure)
    assert converted.species == [{"Na": 1.0}, {"K": 1.0}]
    np.testing.assert_allclose(converted.coords, [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
    disordered = xtalforge.Structure.from_pymatgen(
        FakePmgStructure(
            SILICON,
            [FakePmgSite({"Na": 0.5, "K": 0.5}, frac_coords=[0.0, 0.0, 0.0])],
        )
    ).to_pymatgen()
    assert isinstance(disordered, FakeRxStructure)
    assert disordered.species == [{"Na": 0.5, "K": 0.5}]


def test_pymatgen_real_roundtrip():
    core = pytest.importorskip("pymatgen.core")
    structure = xtalforge.Structure.from_arrays(
        SILICON, [11, 19], frac=[[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
    )
    converted = structure.to_pymatgen()
    assert isinstance(converted, core.Structure)
    back = xtalforge.Structure.from_pymatgen(converted)
    assert back.formula == structure.formula


def test_to_phonopy_atoms_fake_module(monkeypatch):
    package = types.ModuleType("phonopy")
    structure_module = types.ModuleType("phonopy.structure")
    atoms_module = types.ModuleType("phonopy.structure.atoms")
    atoms_module.PhonopyAtoms = FakePhonopyAtoms
    package.structure = structure_module
    structure_module.atoms = atoms_module
    monkeypatch.setitem(sys.modules, "phonopy", package)
    monkeypatch.setitem(sys.modules, "phonopy.structure", structure_module)
    monkeypatch.setitem(sys.modules, "phonopy.structure.atoms", atoms_module)
    converted = xtalforge.Structure.from_arrays(
        SILICON, [14, 14], frac=[[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
    ).to_phonopy_atoms()
    assert isinstance(converted, FakePhonopyAtoms)
    assert converted.kwargs["symbols"] == ["Si", "Si"]
    assert converted.kwargs["cell"] == SILICON


def test_to_phonopy_atoms_real():
    pytest.importorskip("phonopy")
    structure = xtalforge.Structure(CIF)
    converted = structure.to_phonopy_atoms()
    np.testing.assert_allclose(converted.cell, structure.cell, atol=1e-10)
    np.testing.assert_allclose(
        converted.scaled_positions, structure.frac_positions, atol=1e-10
    )
    assert converted.symbols == ["Si"] * 8


def test_get_symmetry_dataset_shape():
    dataset = xtalforge.get_symmetry_dataset(xtalforge.Structure(CIF), symprec=1e-5)
    assert dataset["number"] == 227
    assert xtalforge.get_symmetry_dataset(CIF)["number"] == 227
    assert dataset["hall_number"] == 525
    assert dataset["international"].replace(" ", "") == "Fd-3m"
    assert dataset["pointgroup"] == "m-3m"
    assert len(dataset["rotations"]) == len(dataset["translations"])
    assert len(dataset["equivalent_atoms"]) == 8
    assert len(dataset["std_types"]) == len(dataset["std_positions"])
    assert dataset["std_mapping_to_primitive"] == []


def test_get_symmetry_dataset_matches_spglib():
    spglib = pytest.importorskip("spglib")
    structure = xtalforge.Structure(CIF)
    dataset = xtalforge.get_symmetry_dataset(structure, symprec=1e-5)
    reference = spglib.get_symmetry_dataset(
        (structure.cell, structure.frac_positions, structure.numbers.astype(int)),
        symprec=1e-5,
    )
    assert dataset["number"] == reference.number
    assert dataset["wyckoffs"] == list(reference.wyckoffs)


def test_from_rdkit_fake_mol():
    mol = FakeRdkitMol(
        [FakeRdkitAtom("C"), FakeRdkitAtom("O", charge=-1)],
        [FakeRdkitBond(0, 1, 1.0)],
        [[0.0, 0.0, 0.0], [1.43, 0.0, 0.0]],
    )
    structure = xtalforge.Structure.from_rdkit(mol)
    data = document(structure)
    assert [site["element"] for site in data["sites"]] == ["C", "O"]
    assert data["sites"][1]["charge"] == -1.0
    assert data["explicit_bonds"] == [
        {"i": 0, "j": 1, "image": [0, 0, 0], "order": 1.0, "remove": False}
    ]
    assert interop(structure)["molecule"] is True


def test_to_rdkit_fake_module(monkeypatch):
    types_module = types.SimpleNamespace(
        SINGLE="single", DOUBLE="double", TRIPLE="triple", AROMATIC="aromatic"
    )
    chem = types.ModuleType("rdkit.Chem")
    chem.RWMol = FakeRxMol
    chem.Atom = FakeRxAtom
    chem.Conformer = FakeRxConformer
    chem.BondType = types_module
    chem.SanitizeMol = lambda mol: mol
    geometry = types.ModuleType("rdkit.Geometry")
    geometry.Point3D = lambda x, y, z: (x, y, z)
    package = types.ModuleType("rdkit")
    package.Chem = chem
    package.Geometry = geometry
    monkeypatch.setitem(sys.modules, "rdkit", package)
    monkeypatch.setitem(sys.modules, "rdkit.Chem", chem)
    monkeypatch.setitem(sys.modules, "rdkit.Geometry", geometry)
    structure = xtalforge.Structure.from_rdkit(
        FakeRdkitMol(
            [FakeRdkitAtom("C"), FakeRdkitAtom("O", charge=-1)],
            [FakeRdkitBond(0, 1, 2.0)],
            [[0.0, 0.0, 0.0], [1.43, 0.0, 0.0]],
        )
    )
    converted = structure.to_rdkit()
    assert isinstance(converted, FakeRxMol)
    assert [atom.symbol for atom in converted.atoms] == ["C", "O"]
    assert converted.atoms[1].charge == -1
    assert converted.bonds == [(0, 1, "double")]


def test_rdkit_real_roundtrip():
    pytest.importorskip("rdkit")
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.AddHs(Chem.MolFromSmiles("CO"))
    assert AllChem.EmbedMolecule(mol, randomSeed=7) == 0
    structure = xtalforge.Structure.from_rdkit(mol)
    back = structure.to_rdkit()
    assert back.GetNumAtoms() == mol.GetNumAtoms()
    assert len(back.GetBonds()) == len(mol.GetBonds())


def test_package_import_does_not_import_optional_dependencies():
    code = (
        "import sys, xtalforge; "
        "print([name for name in ('ase', 'pymatgen', 'phonopy', 'rdkit') "
        "if name in sys.modules])"
    )
    output = subprocess.check_output([sys.executable, "-c", code], text=True)
    assert output.strip() == "[]"
