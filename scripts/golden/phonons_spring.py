#!/usr/bin/env python3
"""Generates the phonopy reference for the spring-model phonon test.

The model is a simple cubic lattice with one carbon atom, nearest-neighbour
springs of 10 eV/Angstrom^2 and equilibrium length 4 Angstrom. Phonopy reads
the analytic supercell force constants and reports the frequencies on a
2x2x2 Monkhorst-Pack mesh. Run with a virtual environment that has phonopy:

    python scripts/golden/phonons_spring.py
"""
import json
from pathlib import Path

import numpy as np
import phonopy
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms

K = 10.0
A = 4.0
CELL = PhonopyAtoms(symbols=["C"], cell=np.eye(3) * A, scaled_positions=[[0, 0, 0]])
PHONOPY = Phonopy(CELL, supercell_matrix=np.diag([3, 3, 3]))
POSITIONS = PHONOPY.supercell.scaled_positions
LATTICE = np.array(PHONOPY.supercell.cell)
N = len(POSITIONS)
FC = np.zeros((N, N, 3, 3))
for i in range(N):
    for j in range(N):
        if i == j:
            continue
        delta = POSITIONS[j] - POSITIONS[i]
        delta -= np.round(delta)
        cart = delta @ LATTICE
        distance = np.linalg.norm(cart)
        if abs(distance - A) < 1e-6:
            unit = cart / distance
            FC[i, j] = -K * np.outer(unit, unit)
for i in range(N):
    FC[i, i] = -FC[i].sum(axis=0)
PHONOPY.force_constants = FC
POINTS = [
    [x, y, z]
    for x in (-0.25, 0.25)
    for y in (-0.25, 0.25)
    for z in (-0.25, 0.25)
]
FREQUENCIES = [
    [float(value) * 33.35641 for value in PHONOPY.get_frequencies(q)]
    for q in POINTS
]
OUT = {
    "model": "simple cubic nearest-neighbour springs",
    "spring_ev_per_a2": K,
    "equilibrium_a": A,
    "supercell": [3, 3, 3],
    "q_points": POINTS,
    "frequencies_cm1": FREQUENCIES,
    "versions": {"phonopy": phonopy.__version__},
}
TARGET = Path(__file__).resolve().parents[2] / "tests/golden/phase-4/phonons/spring-sc.json"
TARGET.parent.mkdir(parents=True, exist_ok=True)
TARGET.write_text(json.dumps(OUT, indent=1) + "\n")
print(f"wrote {TARGET}")
