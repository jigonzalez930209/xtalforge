#!/usr/bin/env python3
"""Generates the Pawley golden data (Task 7.6.2).

A cubic cell with one point scatterer allows every reflection, so the peak
list is a property of the cell alone. The pattern is a sum of area-normalised
Gaussians with known integrated intensities, a linear background and a zero
shift; the Rust Pawley fit must recover all three.

    python scripts/golden/pawley.py
"""

import json
import math
from pathlib import Path

WAVELENGTH = 1.54184
FWHM = 0.1
ZERO = 0.02
CELL = (4.0, 4.0, 4.0, 90.0, 90.0, 90.0)
RANGE = (10.0, 60.0)
STEP = 0.02
MAX_INDEX = 4


def two_theta(hkl):
    h, k, l = hkl
    a = CELL[0]
    q = (h * h + k * k + l * l) / (a * a)
    sine = math.sqrt(q) * WAVELENGTH / 2.0
    if sine >= 1.0:
        return None
    return 2.0 * math.degrees(math.asin(sine))


def intensity(hkl):
    h, k, l = hkl
    return 1000.0 / (1.0 + h * h + k * k + l * l)


def reflections():
    """Every signed hkl: a P1 cell has no Friedel merging."""
    found = []
    for h in range(-MAX_INDEX, MAX_INDEX + 1):
        for k in range(-MAX_INDEX, MAX_INDEX + 1):
            for l in range(-MAX_INDEX, MAX_INDEX + 1):
                if h == 0 and k == 0 and l == 0:
                    continue
                value = two_theta((h, k, l))
                if value is None or value < RANGE[0] or value > RANGE[1]:
                    continue
                found.append(
                    {
                        "hkl": [h, k, l],
                        "d": CELL[0] / math.sqrt(h * h + k * k + l * l),
                        "two_theta": value,
                        "intensity": intensity((h, k, l)),
                    }
                )
    found.sort(key=lambda item: (item["two_theta"], item["hkl"]))
    return found


def main():
    reflections_list = reflections()
    normalisation = 2.0 / FWHM * math.sqrt(math.log(2.0) / math.pi)
    count = int((RANGE[1] - RANGE[0]) / STEP) + 1
    peaks = []
    for index in range(count):
        value = RANGE[0] + index * STEP
        total = 5.0 + 0.02 * (value - RANGE[0])
        for reflection in reflections_list:
            offset = value - (reflection["two_theta"] + ZERO)
            scaled = offset / FWHM
            total += reflection["intensity"] * normalisation * math.exp(
                -4.0 * math.log(2.0) * scaled * scaled
            )
        peaks.append([round(value, 6), round(total, 8)])
    golden = {
        "cell": list(CELL),
        "wavelength": WAVELENGTH,
        "fwhm": FWHM,
        "zero_shift": ZERO,
        "background": [5.0, 0.02],
        "reflections": reflections_list,
        "peaks": peaks,
    }
    target = (
        Path(__file__).resolve().parents[2]
        / "tests/golden/phase-7/pawley/one-atom-cubic.json"
    )
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(golden, indent=1) + "\n")
    print(f"wrote {target} ({len(reflections_list)} reflections, {len(peaks)} points)")


if __name__ == "__main__":
    main()
