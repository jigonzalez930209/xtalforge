"""Powder surface smoke tests."""

import math

import xtalforge


def test_powder_index_surface() -> None:
    centers = [19.32, 27.46, 30.21, 36.11, 39.22, 41.27, 44.08, 54.37, 56.68]
    lines = []
    for step in range(3000):
        two_theta = 15.0 + step * 0.02
        intensity = sum(
            math.exp(-(((two_theta - center) / 0.12) ** 2)) / (1.0 + index)
            for index, center in enumerate(centers)
        )
        lines.append(f"{two_theta:.3f} {intensity:.6f}")
    result = xtalforge.powder_index("\n".join(lines))
    assert result["peaks"]
    best = result["candidates"][0]
    assert best["system"] == "tetragonal"
    assert abs(best["a"] - 4.5937) < 0.05
    assert abs(best["c"] - 2.9587) < 0.05
    assert best["m20"] > 0.0


def test_powder_pawley_surface() -> None:
    import json
    from pathlib import Path

    golden = json.loads(
        (Path(__file__).parents[3] / "tests" / "golden" / "phase-7" / "pawley" / "one-atom-cubic.json").read_text(
            encoding="utf-8"
        )
    )
    cell = golden["cell"]
    cif = f"""data_one
_cell_length_a {cell[0]}
_cell_length_b {cell[1]}
_cell_length_c {cell[2]}
_cell_angle_alpha {cell[3]}
_cell_angle_beta {cell[4]}
_cell_angle_gamma {cell[5]}
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
C1 C 0 0 0 1
"""
    pattern = "\n".join(f"{peak[0]} {peak[1]}" for peak in golden["peaks"])
    result = xtalforge.powder_pawley(
        cif,
        pattern,
        refine_cell=True,
        refine_zero=True,
        profile="gaussian",
    )
    assert abs(result["cell"][0] - cell[0]) < 5e-3
    assert abs(result["zero_shift"] - golden["zero_shift"]) < 5e-3
    assert result["rwp"] < 0.02
    assert len(result["reflections"]) == len(golden["reflections"])
    assert all(reflection["esd"] > 0.0 for reflection in result["reflections"])


def test_powder_refine_surface() -> None:
    from pathlib import Path

    root = Path(__file__).parents[3]
    cif = (root / "tests" / "corpus" / "simple" / "TiO2.cif").read_text(
        encoding="utf-8"
    )
    pattern = (
        root / "tests" / "golden" / "phase-7" / "rietveld" / "rutile.xy"
    ).read_text(encoding="utf-8")
    result = xtalforge.powder_refine(
        cif,
        pattern,
        refine_positions=True,
        refine_zero=True,
        refine_cell=True,
        refine_width=True,
        doublet=True,
    )
    assert result["stages"] == [
        "background_scale",
        "zero_cell",
        "width",
        "structure",
    ]
    assert result["rwp"] < 0.3
    assert result["parameter_esds"]
    assert any(name == "zero" for name, _ in result["parameter_esds"])
    assert all(esd >= 0.0 for _, esd in result["parameter_esds"])


def test_powder_automatic_leading_cell_table() -> None:
    from pathlib import Path

    root = Path(__file__).parents[3]
    cif = (root / "tests" / "corpus" / "simple" / "TiO2.cif").read_text(
        encoding="utf-8"
    )
    pattern = (
        root / "tests" / "golden" / "phase-7" / "rietveld" / "rutile.xy"
    ).read_text(encoding="utf-8")
    result = xtalforge.powder_automatic(cif, pattern)
    assert result["steps"][:3] == ["peaks", "index", "leading_cells"]
    rows = result["leading_cells"]
    assert rows
    assert any(row["pawley_rwp"] is not None for row in rows)
    fitted = [row["pawley_rwp"] for row in rows if row["pawley_rwp"] is not None]
    assert all(first <= second * 1.01 for first, second in zip(fitted, fitted[1:]))
    numbers = {
        group["number"]
        for group in rows[0]["space_groups"]
    }
    assert 136 in numbers
    assert all(row["space_groups"] for row in rows)
    numbers = {
        group["number"]
        for row in rows
        for group in row["space_groups"]
    }
    assert numbers
    assert rows[0]["class_symbol"]
    assert rows[0]["class_rank"] >= 1
    assert isinstance(rows[0]["delta_bic"], float)
    first = rows[0]["space_groups"][0]
    assert first["symbol"]
    assert first["centering"] in "PABCIRF"
    assert first["order"] > 0
