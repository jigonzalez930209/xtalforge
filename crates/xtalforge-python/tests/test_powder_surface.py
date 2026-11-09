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
