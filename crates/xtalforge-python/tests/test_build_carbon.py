"""Disordered-carbon recipe surface smoke test."""

import xtalforge


def test_build_carbon_recipe_returns_the_report() -> None:
    result = xtalforge.build_carbon(repeat=(1, 1, 1), seed=0)
    assert len(result["sites"]) > 300
    assert abs(result["density"] - 0.42) < 0.03
    assert result["pieces"] == 1
    assert result["periodicity"] == 3
    assert 1.25 <= result["mean_bond"] <= 1.56
    assert result["census"]
