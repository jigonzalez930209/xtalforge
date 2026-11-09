#!/usr/bin/env python3
"""Generates the Rietveld golden pattern (Task 7.8.3).

Reproduces the Crystal-Builder `_rutile_xy` fixture: rutile in `P4_2/mnm`
with the oxygen at (0.29, 0.29, 0), Cu K alpha 1 + K alpha 2 from rietx's
Bragg-Brentano preset, a flat 100-count background and Poisson noise with
seed 0 on a 0.02 degree grid from 20 to 80 degrees.

Crystal-Builder and rietx must be importable:

    PYTHONPATH=/path/to/Crystal-Builder python scripts/golden/rietveld.py

The commit that adds the file records the generator; the golden test reads the
committed `.xy`, so the Rust test suite does not need Python.
"""

import numpy as np

from xtal.core.lattice import Lattice
from xtal.core.structure import Structure
from xtal.powder import bridge
from xtal.powder.data import Radiation

X_TRUE = 0.29
CELL = (4.5940, 4.5940, 2.9590, 90.0, 90.0, 90.0)
TWO_THETA = np.arange(20.0, 80.0, 0.02)

TARGET = (
    __file__.rsplit("scripts/golden/rietveld.py", 1)[0]
    + "tests/golden/phase-7/rietveld/rutile.xy"
)


def main() -> None:
    structure = Structure.from_arrays(
        Lattice.from_parameters(*CELL),
        ["Ti", "O"],
        [[0.0, 0.0, 0.0], [X_TRUE, X_TRUE, 0.0]],
        space_group="P4_2/mnm",
    )
    y = bridge.predict(structure, Radiation("cu"), TWO_THETA)
    counts = np.random.default_rng(0).poisson(y / y.max() * 5000 + 100)
    with open(TARGET, "w", encoding="utf-8") as handle:
        handle.write("# synthetic rutile, Cu Ka1+Ka2, Poisson noise, seed 0\n")
        for angle, count in zip(TWO_THETA, counts):
            handle.write(f"{angle:.2f} {int(count)}\n")
    print(f"wrote {TARGET} ({len(counts)} points, max {counts.max()})")


if __name__ == "__main__":
    main()
