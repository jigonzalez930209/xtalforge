"""Adapters to ASE, pymatgen, phonopy, spglib and RDKit (Task 13.1).

The optional dependencies are imported only inside the functions that need
them, so ``import xtalforge`` never imports ASE, pymatgen, phonopy, spglib or
RDKit.

Round trips keep the data the core model has no field for (pbc, tags, initial
charges, magnetic moments, occupancies, ASE constraints, site properties and
the last calculator results) in the core structure metadata under the
``interop`` key. Arrays are stored in site order; a P1 structure is the
canonical exchange object.

The ASE calculator is exposed as ``xtalforge.AseCalculator``; this module also
offers the ``UffAseCalculator`` alias.
"""

from __future__ import annotations

import json

from ._interop_json import (
    build_document as _build_document,
    cart as _cart,
    filled as _filled,
    interop as _interop,
    jsonable as _jsonable,
)
from ._rdkit_adapter import (
    rdkit_to_structure_json as rdkit_to_structure_json,
    structure_to_rdkit as structure_to_rdkit,
)
from ._spglib_adapter import (
    get_symmetry_dataset as get_symmetry_dataset,
    structure_to_phonopy_atoms as structure_to_phonopy_atoms,
)
from ._xtalforge import InputError, Structure


# --------------------------------------------------------------------- ASE


def _ase_cell(atoms):
    cell = atoms.get_cell() if hasattr(atoms, "get_cell") else atoms.cell
    rows = [[float(value) for value in row] for row in cell]
    if len(rows) != 3:
        raise InputError("invalid_input", "ASE Atoms cell must have three vectors")
    return rows


def _ase_pbc(atoms):
    pbc = getattr(atoms, "pbc", None)
    if pbc is None:
        pbc = atoms.get_pbc() if hasattr(atoms, "get_pbc") else [True, True, True]
    if not hasattr(pbc, "__len__"):
        return [bool(pbc)] * 3
    return [bool(value) for value in pbc]


def _cell_volume(cell):
    (ax, ay, az), (bx, by, bz), (cx, cy, cz) = cell
    return abs(
        ax * (by * cz - bz * cy)
        - ay * (bx * cz - bz * cx)
        + az * (bx * cy - by * cx)
    )


def _optional_array(arrays, names):
    for name in names:
        if name in arrays:
            return [float(value) for value in arrays[name]]
    return None


def _aligned(values, count):
    """A metadata array only when its length matches the atom count."""
    return values if values is not None and len(values) == count else None


def _ase_constraint(constraint):
    item = {"class": type(constraint).__name__}
    indices = getattr(constraint, "indices", None)
    if indices is None:
        indices = getattr(constraint, "index", None)
    if indices is not None:
        item["indices"] = [int(value) for value in indices]
    for key in ("a", "b", "mask"):
        value = getattr(constraint, key, None)
        if value is not None:
            item[key] = _jsonable(value)
    return item


def _ase_calculator_results(calculator):
    if calculator is None:
        return None
    results = _jsonable(getattr(calculator, "results", {}))
    if not isinstance(results, dict) or not results:
        return None
    return {
        "name": str(getattr(calculator, "name", type(calculator).__name__)),
        "results": results,
    }


def ase_to_structure_json(atoms):
    """Returns a structure document for an ASE ``Atoms`` object."""
    symbols = [str(symbol) for symbol in atoms.get_chemical_symbols()]
    positions = [[float(value) for value in point] for point in atoms.get_positions()]
    cell = _ase_cell(atoms)
    pbc = _ase_pbc(atoms)
    arrays = getattr(atoms, "arrays", None) or {}
    charges = _optional_array(arrays, ("initial_charges", "charges"))
    moments = _optional_array(
        arrays, ("initial_magmoms", "initial_magnetic_moments", "magnetic_moments")
    )
    occupancy = _optional_array(arrays, ("occupancy", "occupancies"))
    tags = [int(value) for value in arrays["tags"]] if "tags" in arrays else None
    labels = (
        [str(value) for value in arrays["labels"]]
        if "labels" in arrays
        else [f"{symbol}{index + 1}" for index, symbol in enumerate(symbols)]
    )
    constraints = [_ase_constraint(item) for item in getattr(atoms, "constraints", [])]
    info = _jsonable(getattr(atoms, "info", {}))
    finite = not any(pbc) or _cell_volume(cell) < 1.0e-8
    interop = {
        "source": "ase",
        "pbc": pbc,
        "molecule": finite,
        "labels": labels,
        "tags": tags,
        "charges": charges,
        "magnetic_moments": moments,
        "occupancy": occupancy,
        "constraints": constraints,
        "calculator": _ase_calculator_results(getattr(atoms, "calc", None)),
        "info": info if isinstance(info, dict) else {},
    }
    entries = [
        {
            "element": symbol,
            "coord": list(positions[index]),
            "occupancy": 1.0 if occupancy is None else occupancy[index],
            "charge": None if charges is None or charges[index] == 0.0 else charges[index],
            "label": labels[index],
            "group": None,
        }
        for index, symbol in enumerate(symbols)
    ]
    document = _build_document("ase", entries, not finite, cell, interop, cartesian=True)
    return json.dumps(document)


def _ase_constraints(items):
    from ase.constraints import FixAtoms

    constraints = []
    for item in items:
        if item.get("class") == "FixAtoms" and "indices" in item:
            constraints.append(FixAtoms(indices=[int(index) for index in item["indices"]]))
    return constraints


def structure_to_ase(structure_json):
    """Builds an ASE ``Atoms`` object from a P1 structure document.

    ``FixAtoms`` constraints are rebuilt; other constraint classes stay in the
    ``interop`` metadata because their constructor arguments are not portable.
    """
    import numpy
    from ase import Atoms
    from ase.calculators.singlepoint import SinglePointCalculator

    document = json.loads(structure_json)["structure"]
    sites = document["sites"]
    interop = _interop(document)
    lattice = document["lattice"]
    symbols = [site["element"] for site in sites]
    positions = [_cart(site["frac"], lattice) for site in sites]
    if interop.get("molecule"):
        origin = interop.get("cart_origin") or [0.0, 0.0, 0.0]
        positions = [
            [point[axis] - float(origin[axis]) for axis in range(3)] for point in positions
        ]
        cell = [[0.0, 0.0, 0.0] for _ in range(3)]
        pbc = [False, False, False]
    else:
        cell = lattice
        pbc = interop.get("pbc") or [True, True, True]
    atoms = Atoms(symbols=symbols, positions=positions, cell=cell, pbc=pbc)

    count = len(sites)
    charges = _aligned(interop.get("charges"), count)
    if charges is None:
        charges = [site.get("charge") for site in sites]
    if charges is not None and any(value not in (None, 0, 0.0) for value in charges):
        atoms.set_initial_charges([float(value or 0.0) for value in charges])
    moments = _aligned(interop.get("magnetic_moments"), count)
    if moments is not None and any(value not in (None, 0, 0.0) for value in moments):
        atoms.set_initial_magnetic_moments([float(value or 0.0) for value in moments])
    tags = _aligned(interop.get("tags"), count)
    if tags is not None:
        atoms.set_tags([int(value) for value in tags])
    occupancy = _aligned(interop.get("occupancy"), count)
    if occupancy is not None:
        atoms.set_array("occupancy", numpy.asarray(occupancy, dtype=float))
    labels = _aligned(interop.get("labels"), count)
    if labels is not None:
        atoms.set_array("labels", numpy.asarray(labels))
    atoms.set_constraint(_ase_constraints(interop.get("constraints") or []))
    info = interop.get("info")
    if isinstance(info, dict):
        atoms.info.update(info)
    calculator = interop.get("calculator")
    if isinstance(calculator, dict) and calculator.get("results"):
        results = dict(calculator["results"])
        for key in ("forces", "stress"):
            if key in results:
                results[key] = numpy.asarray(results[key], dtype=float)
        atoms.calc = SinglePointCalculator(atoms, **results)
    return atoms


# ---------------------------------------------------------------- pymatgen


def _pymatgen_species(site):
    species = getattr(site, "species", None)
    entries = []
    if hasattr(species, "items"):
        for key, value in species.items():
            element = getattr(key, "element", key)
            if not isinstance(element, str):
                element = getattr(key, "symbol", str(element))
            entries.append((str(element), float(value), getattr(key, "oxi_state", None)))
    if not entries:
        element = getattr(site, "species_string", None) or "X"
        occupancy = getattr(site, "occupancy", 1.0)
        entries.append(
            (str(element), 1.0 if occupancy is None else float(occupancy), None)
        )
    return entries


def pymatgen_to_structure_json(value):
    """Returns a structure document for a pymatgen ``Structure`` or ``Molecule``."""
    sites = list(getattr(value, "sites", ()))
    properties = [dict(getattr(site, "properties", {}) or {}) for site in sites]
    lattice = getattr(value, "lattice", None)
    if lattice is not None:
        matrix = [[float(component) for component in row] for row in lattice.matrix]
        coords = [[float(component) for component in site.frac_coords] for site in sites]
        pbc = None
        global_properties = getattr(value, "properties", None)
        if isinstance(global_properties, dict):
            pbc = global_properties.get("xtalforge_pbc")
        if pbc is None:
            pbc = getattr(value, "pbc", None)
        if pbc is None:
            pbc = [True, True, True]
        periodic = any(bool(flag) for flag in pbc)
    else:
        matrix = None
        coords = [[float(component) for component in site.coords] for site in sites]
        pbc = [False, False, False]
        periodic = False
    entries = []
    entry_properties = []
    for index, site in enumerate(sites):
        species = _pymatgen_species(site)
        group = index + 1 if len(species) > 1 else None
        props = properties[index]
        for element, occupancy, charge in species:
            formal = charge if charge is not None else props.get("xtalforge_charge")
            entry_properties.append(props)
            entries.append(
                {
                    "element": element,
                    "coord": coords[index],
                    "occupancy": occupancy,
                    "charge": None if formal in (None, 0, 0.0) else float(formal),
                    "label": props.get("xtalforge_label")
                    or getattr(site, "label", None)
                    or f"{element}{index + 1}",
                    "group": group,
                    "tag": props.get("xtalforge_tag"),
                    "moment": props.get("xtalforge_magnetic_moment"),
                }
            )
    interop = {
        "source": "pymatgen",
        "pbc": [bool(flag) for flag in pbc],
        "molecule": not periodic,
        "site_properties": _jsonable(entry_properties),
        "charges": _filled([entry["charge"] for entry in entries]),
        "tags": _filled([entry.get("tag") for entry in entries]),
        "magnetic_moments": _filled([entry.get("moment") for entry in entries]),
    }
    document = _build_document("pymatgen", entries, periodic, matrix, interop)
    return json.dumps(document)


def _pymatgen_groups(sites):
    order, groups = [], {}
    for index, site in enumerate(sites):
        group = site.get("disorder_group")
        key = ("group", group) if group is not None else ("site", index)
        if key not in groups:
            groups[key] = []
            order.append(key)
        groups[key].append((index, site))
    return [groups[key] for key in order]


def _pymatgen_key(species, element, charge):
    if charge not in (None, 0, 0.0) and float(charge).is_integer():
        return species(element, int(charge))
    return element


def structure_to_pymatgen(structure_json):
    """Builds a pymatgen ``Structure``, or a ``Molecule`` for finite systems."""
    from pymatgen.core import Lattice as PmgLattice
    from pymatgen.core import Molecule, Species, Structure as PmgStructure

    document = json.loads(structure_json)["structure"]
    sites = document["sites"]
    interop = _interop(document)
    lattice = document["lattice"]
    stored = interop.get("site_properties") or []
    tags = interop.get("tags")
    moments = interop.get("magnetic_moments")
    origin = interop.get("cart_origin") or [0.0, 0.0, 0.0]
    species_list, coords, properties = [], [], []
    for group in _pymatgen_groups(sites):
        index, first = group[0]
        entries = {}
        for _, site in group:
            key = _pymatgen_key(Species, site["element"], site.get("charge"))
            occupancy = site.get("occupancy")
            entries[key] = entries.get(key, 0.0) + (
                1.0 if occupancy is None else float(occupancy)
            )
        site_properties = dict(first.get("props") or {})
        if len(stored) == len(sites) and isinstance(stored[index], dict):
            site_properties = {**stored[index], **site_properties}
        if first.get("charge") not in (None, 0, 0.0):
            site_properties["xtalforge_charge"] = float(first["charge"])
        if tags is not None and len(tags) == len(sites) and tags[index]:
            site_properties["xtalforge_tag"] = int(tags[index])
        if moments is not None and len(moments) == len(sites) and moments[index]:
            site_properties["xtalforge_magnetic_moment"] = float(moments[index])
        site_properties["xtalforge_label"] = first["label"]
        species_list.append(entries)
        properties.append(site_properties)
        point = _cart(first["frac"], lattice)
        if interop.get("molecule"):
            coords.append([point[axis] - float(origin[axis]) for axis in range(3)])
        else:
            coords.append(first["frac"])
    keys = {key for prop in properties for key in prop}
    site_properties = {
        key: [prop.get(key) for prop in properties] for key in sorted(keys)
    }
    if interop.get("molecule"):
        return Molecule(
            species_list, coords, site_properties=site_properties or None
        )
    flags = [bool(flag) for flag in (interop.get("pbc") or [True, True, True])]
    arguments = {"properties": {"xtalforge_pbc": flags}} if not all(flags) else {}
    return PmgStructure(
        PmgLattice(lattice),
        species_list,
        coords,
        coords_are_cartesian=False,
        site_properties=site_properties or None,
        validate_proximity=False,
        **arguments,
    )


# --------------------------------------------------------------- calculator


_ASE_CALCULATOR = None


def _ase_calculator_class():
    global _ASE_CALCULATOR
    if _ASE_CALCULATOR is not None:
        return _ASE_CALCULATOR
    import numpy
    from ase.calculators.calculator import Calculator

    class AseCalculator(Calculator):
        """ASE calculator backed by the native UFF implementation.

        Exposed as ``xtalforge.AseCalculator``; ``UffAseCalculator`` is an
        alias in ``xtalforge.adapters``. Energy is in eV and forces in
        eV/Angstrom. Each evaluation converts the ``Atoms`` object to a P1
        ``Structure`` and calls the native UFF model.
        """

        name = "xtalforge-uff"
        implemented_properties = ["energy", "forces"]

        def calculate(self, atoms=None, properties=("energy",), system_changes=None):
            if atoms is None:
                atoms = self.atoms
            if atoms is None:
                raise InputError("invalid_input", "AseCalculator needs an Atoms object")
            structure = Structure._from_json(ase_to_structure_json(atoms))
            energy, forces = structure._ff_energy_forces()
            self.results = {
                "energy": float(energy),
                "forces": numpy.asarray(forces, dtype=float).reshape(-1, 3),
            }

    _ASE_CALCULATOR = AseCalculator
    AseCalculator.__qualname__ = "AseCalculator"
    return _ASE_CALCULATOR


def __getattr__(name):
    if name in ("AseCalculator", "UffAseCalculator"):
        return _ase_calculator_class()
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
