#!/usr/bin/env python3
"""Generate Phase 1 golden data from Crystal-Builder.

The generator is deliberately separate from the Rust tests: it runs only in a
virtual environment with Crystal-Builder installed, records the exact versions
and options, and writes deterministic JSON. CI never regenerates golden data
(Task 0.9.5).
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import platform
import subprocess
import sys
from pathlib import Path
from typing import Any

import gemmi
import numpy as np
import scipy
import spglib
from xtal.core import bonding, groups, measure, p1, rings, slab, supercell
from xtal.io.cif_reader import read_cif

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CORPUS = ROOT / "tests" / "corpus"
DEFAULT_OUT = ROOT / "tests" / "golden" / "phase-1"
FEATURES = ["cell", "sites", "p1", "supercell", "slab", "bonds", "rings", "groups", "measures"]
LATER_PHASES = [
    "space-group-dataset",
    "wyckoff-letters",
    "k-path",
    "pxrd",
    "uff",
    "charges",
    "porosity",
    "topology",
]


def attempt(function: Any) -> Any:
    """Run one feature and record its failure instead of aborting the corpus."""
    try:
        return function()
    except Exception as error:
        return {"error": {"type": type(error).__name__, "message": str(error)}}


def jsonable(value: Any) -> Any:
    """Convert numpy scalars and arrays to JSON values."""
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, tuple):
        return [jsonable(item) for item in value]
    if isinstance(value, list):
        return [jsonable(item) for item in value]
    if isinstance(value, dict):
        return {str(key): jsonable(item) for key, item in value.items()}
    return value


def versions(crystal_builder: Path) -> dict[str, Any]:
    try:
        commit = subprocess.run(
            ["git", "-C", str(crystal_builder), "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        commit = "unknown"
    try:
        package_version = importlib.metadata.version("crystal-builder")
    except importlib.metadata.PackageNotFoundError:
        package_version = "unknown"
    return {
        "crystal_builder": package_version,
        "crystal_builder_commit": commit,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "gemmi": gemmi.__version__,
        "spglib": spglib.__version__,
    }


def cell_data(structure: Any) -> dict[str, Any]:
    lattice = structure.lattice
    return {
        "parameters": [float(value) for value in lattice.parameters],
        "volume": float(lattice.volume),
        "matrix": [[float(value) for value in row] for row in lattice.matrix],
    }


def site_data(structure: Any) -> list[dict[str, Any]]:
    out = []
    for site in structure.sites:
        out.append(
            {
                "label": site.label,
                "element": site.element,
                "frac": [float(value) for value in site.frac],
                "occupancy": float(site.occupancy),
                "charge": jsonable(site.charge),
                "wyckoff": site.wyckoff,
            }
        )
    return out


def p1_data(structure: Any) -> dict[str, Any]:
    cell = p1.expand(structure)
    return {
        "n_atoms": int(cell.n_atoms),
        "lattice": cell_data(cell),
        "elements": list(cell.elements),
        "labels": list(cell.labels),
        "frac": [[float(value) for value in row] for row in cell.frac],
        "occupancy": [float(value) for value in cell.occupancy],
        "multiplicities": [int(cell.multiplicity(index)) for index in range(structure.n_sites)],
        "site_idx": [int(value) for value in cell.site_idx],
    }


def supercell_data(structure: Any, max_atoms: int) -> dict[str, Any] | None:
    for repeat in (2, 1):
        if structure.n_sites * (repeat**3) > max_atoms:
            continue
        cell = supercell.supercell(structure, repeat, repeat, repeat)
        return {
            "repeat": [repeat, repeat, repeat],
            "n_sites": int(cell.n_sites),
            "lattice": cell_data(cell),
        }
    return None


def slab_data(structure: Any, max_atoms: int) -> dict[str, Any] | None:
    if structure.n_sites > max_atoms:
        return None
    result = slab.make_slab(structure, (1, 0, 0), layers=1, vacuum=10.0)
    return {
        "hkl": [1, 0, 0],
        "layers": 1,
        "vacuum": 10.0,
        "thickness": float(result.thickness),
        "n_sites": int(result.structure.n_sites),
        "lattice": cell_data(result.structure),
    }


def bonds_data(structure: Any, max_atoms: int) -> list[dict[str, Any]] | None:
    if structure.n_sites > max_atoms:
        return None
    graph = bonding.graph(structure)
    return [
        {
            "i": int(bond.i),
            "j": int(bond.j),
            "image": [int(value) for value in bond.image],
            "distance": float(bond.distance),
            "order": float(bond.order),
            "stated": bool(bond.stated),
        }
        for bond in graph.bonds
    ]


def rings_data(structure: Any, max_atoms: int, limit: int = 2000) -> dict[str, Any] | None:
    if structure.n_sites > max_atoms:
        return None
    found = rings.rings_of(structure, max_size=8)
    paths = [
        [[int(index), [int(value) for value in image]] for index, image in ring]
        for ring in found[:limit]
    ]
    return {"count": len(found), "truncated": len(found) > limit, "paths": paths}


def groups_data(structure: Any, max_atoms: int) -> dict[str, Any] | None:
    if structure.n_sites > max_atoms:
        return None
    found = {}
    for kind in sorted(groups.KINDS):
        matches = groups.find(structure, kind)
        found[kind] = {
            "count": len(matches),
            "atoms": [list(match.atoms) for match in matches[:20]],
        }
    return found


def measures_data(structure: Any) -> dict[str, Any]:
    cell = p1.expand(structure)
    limit = min(cell.n_atoms, 5)
    distances = []
    for i in range(limit):
        for j in range(i + 1, limit):
            distances.append(
                {
                    "i": i,
                    "j": j,
                    "distance": float(measure.distance(cell, cell.lattice, i, j)),
                }
            )
    angles = []
    for i in range(min(cell.n_atoms, 3)):
        for j in range(min(cell.n_atoms, 3)):
            for k in range(j + 1, min(cell.n_atoms, 3)):
                if len({i, j, k}) < 3:
                    continue
                angles.append(
                    {
                        "i": i,
                        "j": j,
                        "k": k,
                        "angle": float(measure.angle(cell, cell.lattice, i, j, k)),
                    }
                )
    return {"distances": distances, "angles": angles}


def dump(path: Path, corpus: Path, options: argparse.Namespace) -> dict[str, Any]:
    relative = path.relative_to(corpus).as_posix()
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    metadata = {
        "source": relative,
        "sha256": digest,
        "features": FEATURES,
        "skipped_features": LATER_PHASES,
        "options": {
            "max_atoms": options.max_atoms,
            "slab_max_atoms": options.slab_max_atoms,
            "symprec": options.symprec,
        },
        "versions": versions(options.crystal_builder),
    }
    try:
        structure = read_cif(path)
    except Exception as error:  # The malformed edge case is expected to land here.
        message = str(error).replace(str(path), relative)
        return {**metadata, "error": {"type": type(error).__name__, "message": message}}

    data = {
        **metadata,
        "formula": structure.composition(),
        "space_group": {
            "hall": structure.space_group.hall,
            "number": int(structure.space_group.number),
            "hm": structure.space_group.hm,
        },
        "n_sites": int(structure.n_sites),
        "cell": cell_data(structure),
        "sites": site_data(structure),
        "p1": attempt(lambda: p1_data(structure)),
        "supercell": attempt(lambda: supercell_data(structure, options.max_atoms)),
        "slab": attempt(lambda: slab_data(structure, options.slab_max_atoms)),
        "bonds": attempt(lambda: bonds_data(structure, options.max_atoms)),
        "rings": attempt(lambda: rings_data(structure, options.max_atoms)),
        "groups": attempt(lambda: groups_data(structure, options.max_atoms)),
        "measures": attempt(lambda: measures_data(structure)),
    }
    return jsonable(data)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corpus", type=Path, default=DEFAULT_CORPUS)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument(
        "--crystal-builder",
        type=Path,
        default=Path("/home/lestad/github/Crystal-Builder"),
    )
    parser.add_argument("--max-atoms", type=int, default=600)
    parser.add_argument("--slab-max-atoms", type=int, default=200)
    parser.add_argument("--symprec", type=float, default=1e-3)
    args = parser.parse_args()

    files = sorted(
        path
        for path in args.corpus.rglob("*.cif")
        if path.is_file()
    )
    if not files:
        print(f"no CIF files under {args.corpus}", file=sys.stderr)
        return 1

    args.out.mkdir(parents=True, exist_ok=True)
    for path in files:
        relative = path.relative_to(args.corpus).with_suffix(".json")
        target = args.out / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        payload = dump(path, args.corpus, args)
        target.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        try:
            shown = target.relative_to(ROOT)
        except ValueError:
            shown = target
        print(f"wrote {shown}")

    try:
        shown_out = args.out.relative_to(ROOT)
    except ValueError:
        shown_out = args.out
    print(f"generated {len(files)} golden files under {shown_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
