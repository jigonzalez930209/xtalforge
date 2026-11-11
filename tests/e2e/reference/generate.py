#!/usr/bin/env python3
"""Generate the independent reference data for the end-to-end battery.

Every number comes from a standard library, never from xtalforge:

- `gemmi` reads the CIF, expands the symmetry and supplies the cell.
- `spglib` detects the space group and the primitive cell.
- `ase` provides the lattice matrix, the periodic neighbour list and the
  nanotube builder.
- NumPy computes the Bragg angles from the reciprocal metric tensor.

The output (`reference.json`) is committed. The Playwright suite loads the
same sample files in the browser and compares the WASM output against it.

Usage:
    python tests/e2e/reference/generate.py [--check-samples]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter
from itertools import product
from pathlib import Path

import gemmi
import numpy as np
import spglib
from ase import Atoms
from ase.build import nanotube
from ase.geometry import cellpar_to_cell
from ase.neighborlist import neighbor_list

ROOT = Path(__file__).resolve().parents[3]
SAMPLES = ["Si", "NaCl", "diamond", "graphene", "TiO2", "MFI", "MOF-5", "Fe"]
CORPUS = {
    "Si": "simple/Si.cif",
    "NaCl": "simple/NaCl.cif",
    "diamond": "simple/diamond.cif",
    "graphene": "simple/graphene.cif",
    "TiO2": "simple/TiO2.cif",
    "MFI": "simple/MFI.cif",
    "MOF-5": "prepared/MOF-5.cif",
    "Fe": "simple/Fe.cif",
}
CU_KA = 1.54184
MIN_TWO_THETA = 5.0
MAX_TWO_THETA = 60.0
SYMPREC = 1e-3


def element_of(site: gemmi.SmallStructure.Site) -> str:
    return gemmi.Element(site.type_symbol).name


def expanded_sites(structure: gemmi.SmallStructure) -> list[gemmi.SmallStructure.Site]:
    return structure.get_all_unit_cell_sites()


def cell_matrix(structure: gemmi.SmallStructure) -> np.ndarray:
    cell = structure.cell
    return cellpar_to_cell([cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma])


def spglib_cell(
    structure: gemmi.SmallStructure,
) -> tuple[np.ndarray, np.ndarray, list[int]]:
    sites = expanded_sites(structure)
    lattice = cell_matrix(structure)
    positions = np.array(
        [[site.fract.x, site.fract.y, site.fract.z] for site in sites], dtype=float
    )
    numbers = [gemmi.Element(site.type_symbol).atomic_number for site in sites]
    return lattice, positions, numbers


def symmetry_reference(structure: gemmi.SmallStructure) -> dict[str, object]:
    dataset = spglib.get_symmetry_dataset(spglib_cell(structure), symprec=SYMPREC)
    orbits = Counter(int(orbit) for orbit in dataset.equivalent_atoms)
    return {
        "spacegroup": int(dataset.number),
        "wyckoff_multiplicities": sorted(orbits.values()),
        "wyckoff_letters": dict(sorted(Counter(str(letter) for letter in dataset.wyckoffs).items())),
    }


def primitive_volume(structure: gemmi.SmallStructure) -> float:
    primitive = spglib.standardize_cell(
        spglib_cell(structure),
        to_primitive=True,
        no_idealize=True,
        symprec=SYMPREC,
    )
    return abs(float(np.linalg.det(primitive[0])))


def cell_contents(structure: gemmi.SmallStructure) -> dict[str, int]:
    counts: Counter[str] = Counter()
    for site in expanded_sites(structure):
        counts[element_of(site)] += 1
    return dict(sorted(counts.items()))


def density(structure: gemmi.SmallStructure) -> float:
    mass = sum(
        gemmi.Element(site.type_symbol).weight * site.occ
        for site in expanded_sites(structure)
    )
    return mass * 1.660_539_066_60 / structure.cell.volume


def ase_atoms(structure: gemmi.SmallStructure) -> Atoms:
    sites = expanded_sites(structure)
    return Atoms(
        [element_of(site) for site in sites],
        scaled_positions=[[site.fract.x, site.fract.y, site.fract.z] for site in sites],
        cell=cell_matrix(structure),
        pbc=True,
    )


def neighbour_reference(structure: gemmi.SmallStructure) -> dict[str, object]:
    atoms = ase_atoms(structure)
    distances = neighbor_list("d", atoms, 4.0)
    values = np.asarray(distances)
    values = values[values > 1e-6]
    shortest = float(values.min()) if values.size else 0.0
    unique = sorted({round(float(value), 3) for value in values})[:8]
    return {"shortest_distance": shortest, "distances": unique}


def reciprocal_vectors(lattice: np.ndarray) -> np.ndarray:
    return np.linalg.inv(lattice).T


def pxrd_peaks(lattice: np.ndarray, sites) -> list[dict[str, object]]:
    reciprocal = reciprocal_vectors(lattice)
    wavelength = CU_KA
    d_min = wavelength / (2.0 * math.sin(math.radians(MAX_TWO_THETA / 2.0)))
    limits = [
        int(math.ceil(float(np.linalg.norm(vector)) / d_min)) + 1 for vector in lattice
    ]
    peaks: dict[float, dict[str, object]] = {}
    for h, k, l in product(
        range(-limits[0], limits[0] + 1),
        range(-limits[1], limits[1] + 1),
        range(-limits[2], limits[2] + 1),
    ):
        if h == 0 and k == 0 and l == 0:
            continue
        g = np.array([h, k, l], dtype=float) @ reciprocal
        norm = float(np.linalg.norm(g))
        if norm <= 0.0:
            continue
        d = 1.0 / norm
        if wavelength >= 2.0 * d:
            continue
        two_theta = math.degrees(2.0 * math.asin(wavelength / (2.0 * d)))
        if not MIN_TWO_THETA <= two_theta <= MAX_TWO_THETA:
            continue
        key = round(d, 5)
        entry = peaks.setdefault(
            key,
            {
                "d": d,
                "two_theta": two_theta,
                "multiplicity": 0,
                "hkl": [h, k, l],
            },
        )
        entry["multiplicity"] = int(entry["multiplicity"]) + 1
        if (h, k, l) < tuple(entry["hkl"]):
            entry["hkl"] = [h, k, l]
    # Integrated intensity |F|^2 m LP with the gemmi IT92 form factors, as an
    # independent check of the relative intensities.
    for entry in peaks.values():
        h, k, l = entry["hkl"]
        s2 = 1.0 / (4.0 * float(entry["d"]) ** 2)
        factor = sum(
            site.occ
            * site.element.it92.calculate_sf(s2)
            * np.exp(2j * math.pi * (h * site.fract.x + k * site.fract.y + l * site.fract.z))
            for site in sites
        )
        theta = math.radians(float(entry["two_theta"]) / 2.0)
        lp = (1.0 + math.cos(2.0 * theta) ** 2) / (math.sin(theta) ** 2 * math.cos(theta))
        entry["intensity"] = float(abs(factor) ** 2) * int(entry["multiplicity"]) * lp
    strongest = max((float(entry["intensity"]) for entry in peaks.values()), default=1.0) or 1.0
    for entry in peaks.values():
        entry["relative_intensity"] = float(entry["intensity"]) / strongest
    return sorted(peaks.values(), key=lambda item: float(item["two_theta"]))


def sample_reference(name: str, path: Path) -> dict[str, object]:
    structure = gemmi.read_small_structure(str(path))
    lattice = cell_matrix(structure)
    sites = expanded_sites(structure)
    cell = structure.cell
    symmetry = symmetry_reference(structure)
    reference = {
        "name": name,
        "cell_parameters": [
            cell.a,
            cell.b,
            cell.c,
            cell.alpha,
            cell.beta,
            cell.gamma,
        ],
        "cell_matrix": lattice.tolist(),
        "volume": cell.volume,
        "density": density(structure),
        "formula": cell_contents(structure),
        "n_atoms": len(sites),
        "spacegroup": symmetry["spacegroup"],
        "wyckoff_multiplicities": symmetry["wyckoff_multiplicities"],
        "wyckoff_letters": symmetry["wyckoff_letters"],
        "primitive_volume": primitive_volume(structure),
        "atom_positions": [
            [
                element_of(site),
                float(site.fract.x),
                float(site.fract.y),
                float(site.fract.z),
            ]
            for site in sites
        ],
        "neighbours": neighbour_reference(structure),
        "pxrd": pxrd_peaks(lattice, sites),
    }
    return reference


def builder_reference() -> dict[str, object]:
    graphene_bond = 1.42
    ase_tube = nanotube(5, 5, length=1, bond=graphene_bond, vacuum=0.0)
    centre = ase_tube.positions[:, :2].mean(axis=0)
    tube_radius = float(
        np.linalg.norm(ase_tube.positions[:, :2] - centre, axis=1).mean()
    )
    return {
        "graphene": {
            "n_atoms": 32,
            "bond_length": graphene_bond,
            "cell_a": graphene_bond * math.sqrt(3.0) * 4,
        },
        "nanotube_5_5": {
            "n_atoms": int(len(ase_tube)),
            "radius": tube_radius,
            "bond_length": graphene_bond,
        },
        "nanotube_5_5_length_2": {
            "n_atoms": int(len(ase_tube)) * 2,
            "radius": tube_radius,
        },
        "c60": {"n_atoms": 60, "n_bonds": 90},
        "c20": {"n_atoms": 20, "n_bonds": 30},
    }


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_samples() -> None:
    for name in SAMPLES:
        corpus = ROOT / "tests" / "corpus" / CORPUS[name]
        docs = ROOT / "docs" / "public" / "samples" / f"{name}.cif"
        if digest(corpus) != digest(docs):
            raise SystemExit(
                f"{docs} differs from {corpus}; refresh the playground sample copy"
            )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-samples", action="store_true")
    parser.add_argument(
        "--out",
        default=str(ROOT / "tests" / "e2e" / "reference" / "reference.json"),
    )
    args = parser.parse_args()
    if args.check_samples:
        check_samples()
    reference = {
        "generated_by": "tests/e2e/reference/generate.py (gemmi, spglib, ASE, NumPy)",
        "constants": {
            "cu_ka_wavelength": CU_KA,
            "two_theta_range": [MIN_TWO_THETA, MAX_TWO_THETA],
            "symprec": SYMPREC,
        },
        "samples": {
            name: sample_reference(name, ROOT / "tests" / "corpus" / CORPUS[name])
            for name in SAMPLES
        },
        "builders": builder_reference(),
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(reference, indent=2, sort_keys=True) + "\n")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
