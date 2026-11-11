"""Compares the CLI, Python and WASM surface outputs field by field.

Usage: python tests/surfaces/compare.py CLI.json PY.json WASM.json
Exits non-zero when any field differs beyond a 1e-9 relative tolerance.
"""

import json
import math
import sys


def close(first, second):
    if isinstance(first, list) and isinstance(second, list):
        return len(first) == len(second) and all(close(a, b) for a, b in zip(first, second))
    if isinstance(first, (int, float)) and isinstance(second, (int, float)):
        return math.isclose(first, second, rel_tol=1e-9, abs_tol=1e-9)
    return first == second


surfaces = {name: json.load(open(path)) for name, path in zip(["cli", "python", "wasm"], sys.argv[1:4])}
differences = []
fields = 0
for sample, values in surfaces["cli"].items():
    keys = set().union(*(surfaces[name][sample].keys() for name in surfaces))
    for key in sorted(keys):
        fields += 1
        reference = values.get(key)
        for name in ("python", "wasm"):
            if not close(reference, surfaces[name][sample].get(key)):
                differences.append((sample, key, name))
for difference in differences:
    print("DIFF", *difference)
print(f"{fields} fields compared across CLI, Python and WASM; {len(differences)} differences")
sys.exit(1 if differences else 0)
