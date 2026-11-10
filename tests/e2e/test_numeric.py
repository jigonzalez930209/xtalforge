"""Numeric end-to-end checks against gemmi, spglib and ASE references."""

from __future__ import annotations

import math

import numpy as np
import pytest

from support import (
    REFERENCE,
    origin_vectors,
    corners_match,
    formula_counts,
    operation,
    open_sample,
    periodic_match,
    projection_error,
    pxrd,
    reciprocal_spacing,
)

SAMPLES = ["Si", "NaCl", "diamond", "graphene", "TiO2", "MFI", "MOF-5", "Fe"]


@pytest.mark.parametrize("sample", SAMPLES)
def test_sample_matches_standard_library(page, base_url, sample):
    state = open_sample(page, base_url, sample)
    reference = REFERENCE["samples"][sample]
    summary = state["summary"]
    assert summary is not None, state.get("error")
    assert formula_counts(summary["formula"]) == reference["formula"]
    assert summary["volume"] == pytest.approx(reference["volume"], rel=1e-8)
    assert summary["density"] == pytest.approx(reference["density"], rel=1e-6)
    assert summary["n_atoms"] == reference["n_atoms"]
    assert state["symmetry"]["number"] == reference["spacegroup"]
    if sample != "MOF-5":
        # MOF-5 is stored in its primitive rhombohedral cell while the detected
        # group is the conventional Fm-3m setting; the Wyckoff multiplicity is
        # reported in the conventional cell, so the two orbit lists differ by
        # the fourfold F centring.
        assert sorted(state["wyckoffOrbits"]) == reference["wyckoff_multiplicities"]
    assert corners_match(state, reference["cell_matrix"]) < 2e-3

    _, worst = periodic_match(
        state["atomPositions"], reference["atom_positions"], reference["cell_matrix"]
    )
    assert worst < 5e-3, f"{sample}: atom positions differ by {worst} A"

    assert projection_error(page, state, sample) < 2e-3
    if sample != "Fe":
        shortest = min(state["bondDistances"])
        assert shortest == pytest.approx(
            reference["neighbours"]["shortest_distance"], abs=5e-3
        )


def test_supercell_matches_reference(page, base_url):
    state = open_sample(page, base_url, "Si")
    reference = REFERENCE["samples"]["Si"]
    result = operation(page, "supercell", na=2, nb=2, nc=2)
    assert result["summary"]["n_atoms"] == 8 * reference["n_atoms"]
    assert result["summary"]["volume"] == pytest.approx(reference["volume"] * 8, rel=1e-8)
    assert len(result["atomPositions"]) == 8 * reference["n_atoms"]
    assert formula_counts(result["summary"]["formula"]) == {
        "Si": 8 * reference["formula"]["Si"]
    }
    assert result["scene"]["atoms"] == 8 * reference["n_atoms"]


def test_slab_geometry_matches_reference(page, base_url):
    state = open_sample(page, base_url, "Si")
    reference = REFERENCE["samples"]["Si"]
    result = operation(page, "slab", h=1, k=1, l=1, layers=2, vacuum=10.0)
    assert result["summary"]["n_atoms"] == 2 * reference["n_atoms"]
    spacing = reciprocal_spacing(reference["cell_matrix"], (1, 1, 1))
    thickness = 2 * spacing
    height = thickness + 10.0
    lengths = [float(np.linalg.norm(vector)) for vector in origin_vectors(result)]
    assert min(abs(length - height) for length in lengths) < 2e-3

    area = reference["volume"] / spacing
    assert result["summary"]["volume"] == pytest.approx(area * height, rel=1e-6)

    normal = min(origin_vectors(result), key=lambda vector: abs(np.linalg.norm(vector) - height))
    positions = np.array([atom["position"] for atom in result["atomPositions"]])
    fractional_z = positions @ normal / float(normal @ normal)
    floor = min(0.5, 10.0 / 2.0) / height
    assert fractional_z.min() >= floor - 1e-6
    assert fractional_z.max() <= floor + thickness / height + 1e-6


def test_symmetry_transforms_match_reference(page, base_url):
    open_sample(page, base_url, "Si")
    reference = REFERENCE["samples"]["Si"]

    primitive = operation(page, "primitive", tolerance=0.001)
    assert primitive["summary"]["volume"] == pytest.approx(
        reference["primitive_volume"], rel=1e-6
    )
    assert primitive["summary"]["n_atoms"] == 2

    open_sample(page, base_url, "Si")
    standard = operation(page, "standardize", tolerance=0.001)
    assert standard["summary"]["volume"] == pytest.approx(reference["volume"], rel=1e-6)
    assert standard["summary"]["n_atoms"] == reference["n_atoms"]

    open_sample(page, base_url, "Si")
    expanded = operation(page, "p1")
    assert expanded["summary"]["n_atoms"] == reference["n_atoms"]
    assert expanded["summary"]["n_sites"] == reference["n_atoms"]

    open_sample(page, base_url, "Si")
    reduced = operation(page, "reduce")
    assert reduced["summary"]["volume"] == pytest.approx(reference["volume"], rel=1e-6)

    open_sample(page, base_url, "Si")
    asymmetric = operation(page, "asymmetrize", tolerance=0.001)
    assert asymmetric["summary"]["n_sites"] == 1
    assert asymmetric["summary"]["n_atoms"] == reference["n_atoms"]


def test_builders_match_reference(page, base_url):
    open_sample(page, base_url, "Si")
    builders = REFERENCE["builders"]

    graphene = operation(page, "graphene", nx=4, ny=4)
    assert graphene["summary"]["n_atoms"] == builders["graphene"]["n_atoms"]
    assert min(graphene["bondDistances"]) == pytest.approx(
        builders["graphene"]["bond_length"], abs=1e-3
    )
    lengths = [float(np.linalg.norm(vector)) for vector in origin_vectors(graphene)]
    assert min(abs(length - builders["graphene"]["cell_a"]) for length in lengths) < 1e-3

    tube = operation(page, "nanotube", n=5, m=5, length=2)
    assert tube["summary"]["n_atoms"] == builders["nanotube_5_5_length_2"]["n_atoms"]
    positions = np.array([atom["position"] for atom in tube["atomPositions"]])
    radial = positions[:, [0, 2]] - positions[:, [0, 2]].mean(axis=0)
    radius = float(np.linalg.norm(radial, axis=1).mean())
    assert radius == pytest.approx(builders["nanotube_5_5"]["radius"], abs=2e-3)
    assert min(tube["bondDistances"]) == pytest.approx(1.41, abs=0.01)

    c60 = operation(page, "c60")
    assert c60["summary"]["n_atoms"] == builders["c60"]["n_atoms"]
    assert c60["bonds"] == builders["c60"]["n_bonds"]

    c20 = operation(page, "c20")
    assert c20["summary"]["n_atoms"] == builders["c20"]["n_atoms"]
    assert c20["bonds"] == builders["c20"]["n_bonds"]

    mof = operation(page, "mof", topology="pcu", node="Zn4O", linker="benzene")
    assert mof["summary"] is not None, mof.get("error")
    assert mof["summary"]["n_atoms"] > 100
    assert mof["summary"]["volume"] > 1000.0


@pytest.mark.parametrize(
    ("sample", "counts"),
    [
        ("Si", (24, 14, 36)),
        ("NaCl", (24, 14, 36)),
        ("Fe", (14, 12, 24)),
        ("graphene", (12, 8, 18)),
        ("MFI", (8, 6, 12)),
    ],
)
def test_brillouin_zone_matches_the_bravais_lattice(page, base_url, sample, counts):
    open_sample(page, base_url, sample)
    zone = page.evaluate("window.__xf.brillouinZone()")
    assert (len(zone["vertices"]), len(zone["faces"]), len(zone["edges"])) == counts


@pytest.mark.parametrize("sample", ["Si", "NaCl", "TiO2", "Fe", "graphene"])
def test_pxrd_positions_match_bragg_reference(page, base_url, sample):
    open_sample(page, base_url, sample)
    reference = REFERENCE["samples"][sample]["pxrd"]
    pattern = pxrd(page)
    reflections = pattern["reflections"]
    assert reflections, f"{sample}: no reflections"
    for reflection in reflections:
        nearest = min(
            reference,
            key=lambda peak: abs(peak["two_theta"] - reflection["two_theta"]),
        )
        assert abs(nearest["two_theta"] - reflection["two_theta"]) < 0.02
        assert abs(nearest["d"] - reflection["d"]) < 1e-3


def test_invalid_input_reports_the_stable_code(page, base_url):
    open_sample(page, base_url, "Si")
    page.evaluate("window.__xf.loadText('not a structure', 'broken')")
    page.wait_for_function("document.querySelector('.error')", timeout=10_000)
    message = page.locator(".error").inner_text()
    assert ":" in message
    assert "invalid" in message.lower() or "detect" in message.lower()


def test_error_path_keeps_the_last_good_scene(page, base_url):
    open_sample(page, base_url, "Si")
    page.evaluate("window.__xf.loadText('data_x\\n_cell_length_a -1', 'broken')")
    page.wait_for_timeout(200)
    state = page.evaluate("window.__xf.state()")
    assert state["summary"] is None
    assert page.locator(".error").count() == 1
    assert math.isfinite(state.get("bonds", 0))
