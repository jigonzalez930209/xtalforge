"""Integration tests for the Python hello-world surface."""

import json
from pathlib import Path

import pytest
import xtalforge

VECTORS = json.loads(
    (Path(__file__).parents[3] / "tests" / "vectors" / "hello_world.json").read_text(
        encoding="utf-8"
    )
)


def test_version_matches_the_package() -> None:
    assert xtalforge.version() == xtalforge.__version__
    assert isinstance(xtalforge.version(), str)
    assert xtalforge.version()


def test_cell_volume_shared_vectors() -> None:
    for case in VECTORS["cell_volume"]:
        got = xtalforge.cell_volume(
            case["a"],
            case["b"],
            case["c"],
            case["alpha"],
            case["beta"],
            case["gamma"],
        )
        assert abs(got - case["expected"]) <= case["tolerance"]


def test_keyword_arguments() -> None:
    got = xtalforge.cell_volume(
        a=1.0, b=1.0, c=1.0, alpha=90.0, beta=90.0, gamma=90.0
    )
    assert abs(got - 1.0) < 1e-12


def test_element_info_shared_vectors() -> None:
    vectors = json.loads(
        (Path(__file__).parents[3] / "tests" / "vectors" / "elements.json").read_text(
            encoding="utf-8"
        )
    )
    for case in vectors["elements"]:
        info = xtalforge.element_info(case["input"])
        assert info["symbol"] == case["symbol"]
        assert info["z"] == case["z"]
        assert abs(info["mass"] - case["mass"]) <= case["tolerance"]


def test_element_info_accepts_an_atomic_number() -> None:
    info = xtalforge.element_info(26)
    assert info["symbol"] == "Fe"


def test_lattice_shared_vectors() -> None:
    lattice_vectors = json.loads(
        (Path(__file__).parents[3] / "tests" / "vectors" / "lattice.json").read_text(
            encoding="utf-8"
        )
    )
    for case in lattice_vectors["cases"]:
        info = xtalforge.lattice_info(
            case["a"],
            case["b"],
            case["c"],
            case["alpha"],
            case["beta"],
            case["gamma"],
        )
        assert abs(info["volume"] - case["volume"]) <= 1e-9
        result = xtalforge.lattice_minimum_image(
            case["a"],
            case["b"],
            case["c"],
            case["alpha"],
            case["beta"],
            case["gamma"],
            case["minimum_image"]["delta"],
        )
        assert abs(result["distance"] - case["minimum_image"]["distance"]) <= 1e-9


def test_invalid_element_raises_input_error() -> None:
    with pytest.raises(xtalforge.InputError) as caught:
        xtalforge.element_info("Zz")
    assert caught.value.args[0] == "invalid_input"


def test_invalid_cell_raises_input_error() -> None:
    with pytest.raises(xtalforge.InputError) as caught:
        xtalforge.cell_volume(0.0, 1.0, 1.0, 90.0, 90.0, 90.0)
    assert caught.value.args[0] == "invalid_input"
    assert "invalid input" in caught.value.args[1]
    assert isinstance(caught.value, xtalforge.XtalforgeError)
