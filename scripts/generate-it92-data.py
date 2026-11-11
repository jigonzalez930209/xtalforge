"""Writes crates/xtalforge-analysis/data/it92.json, the neutral-atom X-ray
form-factor coefficients of International Tables for Crystallography Vol. C
(1992), Table 6.1.1.4, as tabulated by gemmi (MPL-2.0).

f(s) = sum_i a_i exp(-b_i s^2) + c, with s = sin(theta)/lambda in 1/Angstrom.

Usage: .venv/bin/python scripts/generate-it92-data.py
"""
import json
import pathlib

import gemmi

rows = {}
for z in range(1, 119):
    element = gemmi.Element(z)
    coefficients = element.it92
    if coefficients is None:
        continue
    rows[element.name] = [round(v, 6) for v in list(coefficients.a) + list(coefficients.b) + [coefficients.c]]

target = pathlib.Path(__file__).resolve().parent.parent / "crates/xtalforge-analysis/data/it92.json"
target.write_text(
    json.dumps(
        {
            "source": "International Tables for Crystallography Vol. C (1992), Table 6.1.1.4; via gemmi "
            + gemmi.__version__,
            "layout": "a1 a2 a3 a4 b1 b2 b3 b4 c",
            "coefficients": rows,
        },
        separators=(",", ":"),
    )
    + "\n"
)
print(f"wrote {len(rows)} elements to {target}")
