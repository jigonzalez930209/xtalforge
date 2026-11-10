"""Deterministic 3D screenshots compared with committed golden images."""

from __future__ import annotations

import pytest

from support import (
    MISMATCH_LIMIT,
    compare_image,
    open_sample,
    operation,
    screenshot,
)

SAMPLES = ["Si", "NaCl", "diamond", "graphene", "TiO2", "MFI", "MOF-5", "Fe"]


@pytest.mark.parametrize("sample", SAMPLES)
def test_sample_renders_like_golden(page, base_url, sample):
    open_sample(page, base_url, sample)
    path = screenshot(page, f"sample-{sample}")
    mismatch = compare_image(path, f"sample-{sample}")
    assert mismatch <= MISMATCH_LIMIT, f"{sample}: {mismatch:.3%} pixels differ"


def test_supercell_renders_like_golden(page, base_url):
    open_sample(page, base_url, "Si")
    operation(page, "supercell", na=2, nb=2, nc=2)
    path = screenshot(page, "operation-supercell")
    assert compare_image(path, "operation-supercell") <= MISMATCH_LIMIT


def test_slab_renders_like_golden(page, base_url):
    open_sample(page, base_url, "Si")
    operation(page, "slab", h=1, k=1, l=1, layers=2, vacuum=10.0)
    path = screenshot(page, "operation-slab")
    assert compare_image(path, "operation-slab") <= MISMATCH_LIMIT


def test_graphene_renders_like_golden(page, base_url):
    open_sample(page, base_url, "Si")
    operation(page, "graphene", nx=4, ny=4)
    path = screenshot(page, "operation-graphene")
    assert compare_image(path, "operation-graphene") <= MISMATCH_LIMIT


def test_nanotube_renders_like_golden(page, base_url):
    open_sample(page, base_url, "Si")
    operation(page, "nanotube", n=5, m=5, length=2)
    path = screenshot(page, "operation-nanotube")
    assert compare_image(path, "operation-nanotube") <= MISMATCH_LIMIT


def test_c60_renders_like_golden(page, base_url):
    open_sample(page, base_url, "Si")
    operation(page, "c60")
    path = screenshot(page, "operation-c60")
    assert compare_image(path, "operation-c60") <= MISMATCH_LIMIT


def test_mof_renders_like_golden(page, base_url):
    open_sample(page, base_url, "Si")
    operation(page, "mof", topology="pcu", node="Zn4O", linker="benzene")
    path = screenshot(page, "operation-mof")
    assert compare_image(path, "operation-mof") <= MISMATCH_LIMIT


def test_space_fill_style_renders_like_golden(page, base_url):
    open_sample(page, base_url, "NaCl")
    page.evaluate("window.__xf.setStyle('space-fill', 1.0)")
    path = screenshot(page, "style-space-fill")
    assert compare_image(path, "style-space-fill") <= MISMATCH_LIMIT


def test_screenshots_are_deterministic(page, base_url):
    open_sample(page, base_url, "Si")
    first = screenshot(page, "deterministic-a")
    second = screenshot(page, "deterministic-b")
    assert compare_image(first, "sample-Si") == 0.0
    assert compare_image(second, "sample-Si") == 0.0
