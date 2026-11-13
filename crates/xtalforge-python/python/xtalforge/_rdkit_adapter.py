"""RDKit molecule exchange for the xtalforge Python adapters (Task 13.1.4)."""

from __future__ import annotations

import json

from ._interop_json import build_document as _build_document
from ._interop_json import cart as _cart
from ._interop_json import filled as _filled
from ._interop_json import interop as _interop
from ._xtalforge import InputError


def rdkit_to_structure_json(mol):
    """Returns a structure document for an RDKit ``Mol`` with a conformer."""
    count = int(mol.GetNumAtoms())
    try:
        conformer = mol.GetConformer()
    except (ValueError, RuntimeError) as error:
        raise InputError(
            "invalid_input", "RDKit molecule needs a conformer with 3D coordinates"
        ) from error
    symbols, positions, charges, labels = [], [], [], []
    for index in range(count):
        atom = mol.GetAtomWithIdx(index)
        symbol = str(atom.GetSymbol())
        point = conformer.GetAtomPosition(index)
        positions.append([float(point.x), float(point.y), float(point.z)])
        symbols.append(symbol)
        charges.append(float(atom.GetFormalCharge()))
        labels.append(f"{symbol}{index + 1}")
    bonds = []
    for bond in mol.GetBonds():
        first = int(bond.GetBeginAtomIdx())
        second = int(bond.GetEndAtomIdx())
        order = float(bond.GetBondTypeAsDouble())
        bonds.append(
            {
                "i": min(first, second),
                "j": max(first, second),
                "image": [0, 0, 0],
                "order": order,
                "remove": False,
            }
        )
    entries = [
        {
            "element": symbols[index],
            "coord": positions[index],
            "occupancy": 1.0,
            "charge": None if charges[index] == 0.0 else charges[index],
            "label": labels[index],
            "group": None,
        }
        for index in range(count)
    ]
    interop = {
        "source": "rdkit",
        "pbc": [False, False, False],
        "molecule": True,
        "charges": _filled([entry["charge"] for entry in entries]),
        "labels": labels,
    }
    document = _build_document("rdkit", entries, False, None, interop)
    document["structure"]["explicit_bonds"] = bonds
    return json.dumps(document)


def _bond_type(chem, order):
    types = chem.BondType
    if abs(order - 1.5) < 1.0e-6:
        return types.AROMATIC
    if abs(order - 2.0) < 1.0e-6:
        return types.DOUBLE
    if abs(order - 3.0) < 1.0e-6:
        return types.TRIPLE
    return types.SINGLE


def structure_to_rdkit(structure_json):
    """Builds an RDKit ``Mol`` with a conformer from a P1 structure document."""
    from rdkit import Chem
    from rdkit.Geometry import Point3D

    document = json.loads(structure_json)["structure"]
    sites = document["sites"]
    origin = _interop(document).get("cart_origin") or [0.0, 0.0, 0.0]
    positions = [
        [
            _cart(site["frac"], document["lattice"])[axis] - float(origin[axis])
            for axis in range(3)
        ]
        for site in sites
    ]
    editable = Chem.RWMol()
    indices = []
    for site in sites:
        atom = Chem.Atom(site["element"])
        charge = site.get("charge")
        if charge not in (None, 0, 0.0):
            atom.SetFormalCharge(int(round(float(charge))))
        indices.append(editable.AddAtom(atom))
    aromatic = set()
    for bond in document.get("explicit_bonds") or []:
        if bond.get("remove"):
            continue
        order = float(bond.get("order", 1.0))
        first = indices[int(bond["i"])]
        second = indices[int(bond["j"])]
        editable.AddBond(first, second, _bond_type(Chem, order))
        if abs(order - 1.5) < 1.0e-6:
            aromatic.update((first, second))
    mol = editable.GetMol()
    for index in aromatic:
        mol.GetAtomWithIdx(index).SetIsAromatic(True)
    conformer = Chem.Conformer(len(sites))
    for index, point in enumerate(positions):
        conformer.SetAtomPosition(index, Point3D(point[0], point[1], point[2]))
    mol.AddConformer(conformer, assignId=True)
    Chem.SanitizeMol(mol)
    return mol
