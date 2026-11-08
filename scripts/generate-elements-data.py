#!/usr/bin/env python3
"""Generate the embedded element bundle from Crystal-Builder.

Run in the project virtual environment with Crystal-Builder installed:

    .venv/bin/python scripts/generate-elements-data.py

The output is reproducible except for the compressed bytes' timestamp, which
is fixed to zero. The Rust test verifies the decompressed payload checksum.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import subprocess
from pathlib import Path

from xtal.core import pyykko
from xtal.core.elements import (
    ELECTRONEGATIVITY,
    VALENCE,
    _NAMES,
    element,
    symbol_from_z,
)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "crates" / "xtalforge-core" / "data" / "elements.json.gz"
CRYSTAL_BUILDER = Path("/home/lestad/github/Crystal-Builder")


def commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "-C", str(CRYSTAL_BUILDER), "rev-parse", "HEAD"],
            text=True,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def row(z: int, symbol: str) -> dict[str, object]:
    record = element(symbol)
    single = pyykko.SINGLE[symbol] / 100.0
    double = pyykko.DOUBLE[symbol]
    return {
        "z": z,
        "symbol": symbol,
        "name": _NAMES[symbol],
        "mass": record.mass,
        "covalent_radius": record.covalent_radius,
        "vdw_radius": record.vdw_radius,
        "color": list(record.color),
        "is_metal": bool(record.is_metal),
        "valence": VALENCE.get(symbol),
        "electronegativity": ELECTRONEGATIVITY.get(symbol),
        "pyykko_single": single,
        "pyykko_double": None if double is None else double / 100.0,
        "pyykko_aromatic": pyykko.aromatic_radius(symbol),
    }


def main() -> None:
    elements = [row(z, symbol_from_z(z)) for z in range(1, 119)]
    deuterium = element("D")
    elements.append(
        {
            "z": 1,
            "symbol": "D",
            "name": "Deuterium",
            "mass": deuterium.mass,
            "covalent_radius": deuterium.covalent_radius,
            "vdw_radius": deuterium.vdw_radius,
            "color": list(deuterium.color),
            "is_metal": bool(deuterium.is_metal),
            "valence": VALENCE.get("D"),
            "electronegativity": ELECTRONEGATIVITY.get("D"),
            "pyykko_single": pyykko.SINGLE["H"] / 100.0,
            "pyykko_double": None,
            "pyykko_aromatic": pyykko.aromatic_radius("H"),
        }
    )
    dummy = element("X")
    elements.append(
        {
            "z": 0,
            "symbol": "X",
            "name": "Dummy",
            "mass": dummy.mass,
            "covalent_radius": dummy.covalent_radius,
            "vdw_radius": dummy.vdw_radius,
            "color": list(dummy.color),
            "is_metal": bool(dummy.is_metal),
            "valence": None,
            "electronegativity": None,
            "pyykko_single": None,
            "pyykko_double": None,
            "pyykko_aromatic": None,
        }
    )
    data = {
        "schema_version": 1,
        "source": (
            f"Crystal-Builder {commit()} core/elements.py, core/pyykko.py, "
            "ff/isotopes.py, ff/vibrations.py"
        ),
        "elements": elements,
        "isotopes": [{"symbol": "D", "element": "H", "mass": deuterium.mass}],
    }
    payload = json.dumps(data, indent=2, sort_keys=True).encode() + b"\n"
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with gzip.GzipFile(filename="", mode="wb", fileobj=open(OUT, "wb"), mtime=0) as stream:
        stream.write(payload)
    print(f"wrote {OUT.relative_to(ROOT)}")
    print(f"sha256 {hashlib.sha256(payload).hexdigest()}")


if __name__ == "__main__":
    main()
