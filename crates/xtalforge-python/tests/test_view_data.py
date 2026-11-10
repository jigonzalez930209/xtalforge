"""Integration tests for the full scene-geometry surface."""

from pathlib import Path

import xtalforge

SI = (Path(__file__).parents[3] / "tests" / "corpus" / "simple" / "Si.cif").read_text(
    encoding="utf-8"
)


def test_scene_data_returns_the_full_geometry() -> None:
    scene = xtalforge.scene_data(SI, "ball-stick", 1.0)
    assert scene["atoms"] and scene["bonds"] and scene["cell_edges"]
    assert len(scene["atoms"]) == 8
    assert len(scene["cell_edges"]) == 12
    first = scene["atoms"][0]
    assert first["element"] == "Si"
    assert len(first["position"]) == 3
    assert first["radius"] > 0
    assert len(first["color"]) == 3
    bond = scene["bonds"][0]
    assert first["index"] in (bond["i"], bond["j"])
    assert scene["settings"]["style"] == "ball-stick"


def test_scene_data_rejects_an_unknown_style() -> None:
    try:
        xtalforge.scene_data(SI, "dotted")
    except xtalforge.InputError as error:
        assert "draw style" in str(error)
    else:
        raise AssertionError("an unknown style must fail")
