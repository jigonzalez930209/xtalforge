"""Shared JSON helpers for the optional Python adapters (Task 13.1)."""

from __future__ import annotations

import json

from ._xtalforge import Structure, element_info

INTEROP = "interop"
MOLECULE_PADDING = 10.0


def jsonable(value):
    """Best-effort conversion of NumPy scalars and arrays to JSON values."""
    if value is None or isinstance(value, (bool, int, float, str)):
        return value
    if isinstance(value, dict):
        return {str(key): jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(item) for item in value]
    tolist = getattr(value, "tolist", None)
    if callable(tolist):
        return jsonable(tolist())
    item = getattr(value, "item", None)
    if callable(item):
        try:
            return jsonable(item())
        except (TypeError, ValueError):
            return None
    return None


def document(structure):
    """The bare structure payload of a ``Structure`` object."""
    return json.loads(structure.to_json())["structure"]


def interop(document_data):
    metadata = document_data.get("metadata") or {}
    value = metadata.get(INTEROP)
    return value if isinstance(value, dict) else {}


def filled(values):
    """A float list with zeros for missing entries, or ``None`` if all missing."""
    if all(value is None for value in values):
        return None
    return [0.0 if value is None else float(value) for value in values]


def cart(frac, lattice):
    return [
        sum(float(frac[axis]) * float(lattice[axis][column]) for axis in range(3))
        for column in range(3)
    ]


def build_document(name, entries, periodic, lattice, interop_data, cartesian=False):
    """Builds a versioned structure document from uniform per-atom entries."""
    numbers = [int(element_info(entry["element"])["z"]) for entry in entries]
    if periodic:
        points = [entry["coord"] for entry in entries]
        probe = (
            Structure.from_arrays(lattice, numbers, cart=points)
            if cartesian
            else Structure.from_arrays(lattice, numbers, frac=points)
        )
        document_data = document(probe)
    else:
        points = [(entry["element"], list(entry["coord"])) for entry in entries]
        probe = Structure._molecule(name, points, MOLECULE_PADDING)
        document_data = document(probe)
        origin = [
            float(probe.positions[0][axis] - entries[0]["coord"][axis]) for axis in range(3)
        ]
        interop_data["cart_origin"] = origin
        interop_data["padding"] = MOLECULE_PADDING
    for site, entry in zip(document_data["sites"], entries):
        site["label"] = entry["label"]
        site["occupancy"] = min(max(float(entry["occupancy"]), 0.0), 1.0)
        site["charge"] = entry["charge"]
        site["disorder_group"] = entry["group"]
    document_data["metadata"][INTEROP] = dict(interop_data, version=1)
    return {"schema_version": 1, "structure": document_data}
