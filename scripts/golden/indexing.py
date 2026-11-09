#!/usr/bin/env python3
"""Generates the powder-indexing golden data (Task 7.5).

Peak positions are derived from the cell with Bragg's law; each unique
reflection up to 100 degrees two-theta is kept and a small deterministic
jitter is added so the figures of merit are finite. The M20 and F_N values are
computed by a second, independent implementation of the definitions in
``docs/algorithms/indexing.md``.

Run with a Python that has numpy:

    python scripts/golden/indexing.py
"""

import json
import math
from pathlib import Path

import numpy as np

WAVELENGTH = 1.540598
TOLERANCE = 0.03
MAX_TWO_THETA = 100.0
MAX_INDEX = 10
MAX_FOM_INDEX = 8

CELLS = [
    ("Si", "cubic", (5.430941, 5.430941, 5.430941, 90.0, 90.0, 90.0)),
    ("TiO2", "tetragonal", (4.5937, 4.5937, 2.9587, 90.0, 90.0, 90.0)),
    ("ZnO", "hexagonal", (3.2417, 3.2417, 5.1876, 90.0, 90.0, 120.0)),
    ("MFI", "orthorhombic", (20.0511, 19.8757, 13.36823, 90.0, 90.0, 90.0)),
    ("monoclinic", "monoclinic", (5.0, 6.0, 7.0, 90.0, 100.0, 90.0)),
]


def lattice_matrix(cell):
    a, b, c, alpha, beta, gamma = cell
    alpha, beta, gamma = map(math.radians, (alpha, beta, gamma))
    return np.array(
        [
            [a, 0.0, 0.0],
            [b * math.cos(gamma), b * math.sin(gamma), 0.0],
            [
                c * math.cos(beta),
                c * (math.cos(alpha) - math.cos(beta) * math.cos(gamma))
                / math.sin(gamma),
                c
                * math.sqrt(
                    1.0
                    - math.cos(alpha) ** 2
                    - math.cos(beta) ** 2
                    - math.cos(gamma) ** 2
                    + 2.0 * math.cos(alpha) * math.cos(beta) * math.cos(gamma)
                )
                / math.sin(gamma),
            ],
        ]
    )


def reciprocal_squared(cell, hkl):
    """Q = |A^-1 v|^2 for a direct matrix whose rows are the lattice vectors."""
    inverse = np.linalg.inv(lattice_matrix(cell))
    return float(np.linalg.norm(inverse @ np.array(hkl, dtype=float)) ** 2)


def two_theta(q):
    sine = min(1.0, math.sqrt(q) * WAVELENGTH / 2.0)
    return 2.0 * math.degrees(math.asin(sine))


def peak_list(cell):
    positions = []
    for h in range(MAX_INDEX + 1):
        for k in range(MAX_INDEX + 1):
            for l in range(MAX_INDEX + 1):
                if h == 0 and k == 0 and l == 0:
                    continue
                value = two_theta(reciprocal_squared(cell, (h, k, l)))
                if 5.0 < value <= MAX_TWO_THETA:
                    positions.append(value)
    positions.sort()
    unique = []
    for value in positions:
        if not unique or abs(value - unique[-1]) >= 0.01:
            unique.append(value)
    seed = 0x12345678
    peaks = []
    for value in unique[:60]:
        seed = (seed * 6364136223846793005 + 1442695040888963407) % (1 << 64)
        jitter = ((seed >> 11) / float(1 << 53) - 0.5) * 0.006
        peaks.append(value + jitter)
    return peaks


def predicted(cell, max_q):
    limit = MAX_FOM_INDEX
    lines = []
    for h in range(-limit, limit + 1):
        for k in range(-limit, limit + 1):
            for l in range(-limit, limit + 1):
                if h < 0 or (h == 0 and k < 0) or (h == 0 and k == 0 and l <= 0):
                    continue
                q = reciprocal_squared(cell, (h, k, l))
                if q <= max_q * (1.0 + 1e-9):
                    lines.append((q, (h, k, l)))
    lines.sort()
    return lines


def nearest(lines, q):
    best = None
    for line in lines:
        delta = abs(line[0] - q)
        if best is None or delta < best[0]:
            best = (delta, line)
    return best[1] if best else None


def figures(cell, peaks):
    observed = sorted((2.0 * math.sin(math.radians(value / 2.0)) / WAVELENGTH) ** 2 for value in peaks)
    lines = predicted(cell, observed[-1])
    window = min(20, len(observed))
    delta_q = 0.0
    delta_theta = 0.0
    matched = 0
    for q in observed[:window]:
        line = nearest(lines, q)
        if line is None:
            continue
        calculated = two_theta(line[0])
        target = two_theta(q)
        if abs(calculated - target) > TOLERANCE:
            continue
        delta_q += abs(line[0] - q)
        delta_theta += abs(calculated - target)
        matched += 1
    q20 = observed[window - 1]
    n20 = sum(1 for line in lines if line[0] <= q20)
    mean_q = max(delta_q / matched, 1e-12)
    mean_theta = max(delta_theta / matched, 1e-12)
    return {
        "m20": q20 / (2.0 * mean_q * n20),
        "f_n": len(observed) / (mean_theta * len(lines)),
        "n20": n20,
        "n_possible": len(lines),
        "mean_delta_q": mean_q,
        "mean_delta_two_theta": mean_theta,
    }


def main():
    entries = []
    for name, system, cell in CELLS:
        peaks = peak_list(cell)
        entries.append(
            {
                "name": name,
                "system": system,
                "cell": dict(
                    zip(
                        ("a", "b", "c", "alpha", "beta", "gamma"),
                        [float(value) for value in cell],
                    )
                ),
                "wavelength": WAVELENGTH,
                "tolerance": TOLERANCE,
                "peaks": [round(value, 8) for value in peaks],
                "figures": figures(cell, peaks),
            }
        )
    target = (
        Path(__file__).resolve().parents[2]
        / "tests/golden/phase-7/indexing/cells.json"
    )
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps({"cells": entries, "max_index": MAX_INDEX}, indent=1) + "\n")
    print(f"wrote {target}")


if __name__ == "__main__":
    main()
