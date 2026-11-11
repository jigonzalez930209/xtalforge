"""Browser checks for the interactive playground panels."""

from __future__ import annotations

from pathlib import Path

import pytest

from support import ARTIFACTS, operation, open_sample, state


def test_tabs_and_brillouin_panel(page, base_url):
    open_sample(page, base_url, "Si")
    for name in ["Symmetry", "Analysis", "Powder", "Force field", "Export", "Structure"]:
        page.get_by_role("button", name=name, exact=True).click()
        assert page.locator(".tab-section").is_visible()
    page.get_by_role("button", name="Symmetry", exact=True).click()
    page.locator("summary", has_text="Brillouin zone").click()
    page.get_by_role("button", name="Show Brillouin zone").click()
    current = state(page)
    assert current["brillouin"] is True
    assert current["zone"] == {"vertices": 24, "faces": 14, "edges": 36}
    page.get_by_role("button", name="Back to the structure").click()
    assert state(page)["brillouin"] is False


def test_atom_selection_reports_measures(page, base_url):
    open_sample(page, base_url, "Si")
    page.evaluate("window.__xf.selectAtoms([0, 1])")
    distance = state(page)["measure"]
    assert distance["kind"] == "distance (Å)"
    assert distance["value"] > 0
    page.evaluate("window.__xf.selectAtoms([0, 1, 2])")
    assert state(page)["measure"]["kind"] == "angle (deg)"
    page.evaluate("window.__xf.selectAtoms([0, 1, 2, 3])")
    assert state(page)["measure"]["kind"] == "torsion (deg)"
    assert len(state(page)["selection"]) == 4


def test_new_builders_from_the_playground(page, base_url):
    open_sample(page, base_url, "Si")
    ribbon = operation(page, "ribbon", width=4, length=4, vacuum=12)
    assert ribbon["summary"]["n_atoms"] == 32
    polymer = operation(page, "polymer", smiles="C=C", monomers=4, density=1.0)
    assert polymer["summary"]["n_atoms"] == 8
    schwarzite = operation(page, "schwarzite", samples=10, cell=12, isovalue=0)
    assert schwarzite["summary"]["n_atoms"] > 0


def test_powder_plot_and_uff_energy(page, base_url):
    open_sample(page, base_url, "Si")
    page.get_by_role("button", name="Powder", exact=True).click()
    radiation = page.locator(".tab-section select").first
    assert radiation.locator("option").count() == 7
    assert radiation.input_value() == "cu-ka"
    page.get_by_role("button", name="Calculate", exact=True).click()
    assert page.locator(".pxrd-plot").count() == 1
    assert page.locator(".tab-section tbody tr").count() >= 3
    page.get_by_role("button", name="Force field", exact=True).click()
    page.get_by_role("button", name="UFF energy", exact=True).click()
    page.wait_for_function(
        "document.querySelector('.tab-section .mono')?.textContent?.includes('Energy =')",
        timeout=30_000,
    )


def test_export_downloads_poscar(page, base_url):
    open_sample(page, base_url, "Si")
    page.get_by_role("button", name="Export", exact=True).click()
    ARTIFACTS.mkdir(parents=True, exist_ok=True)
    target = ARTIFACTS / "exported.poscar"
    with page.expect_download() as download:
        page.get_by_role("button", name="POSCAR", exact=True).click()
    download.value.save_as(str(target))
    text = Path(target).read_text()
    assert len(text.splitlines()) >= 8
    assert "Si" in text


def test_fullscreen_toggle(page, base_url):
    open_sample(page, base_url, "Si")
    page.get_by_role("button", name="Fullscreen", exact=True).click()
    page.wait_for_function("document.fullscreenElement !== null", timeout=10_000)
    assert page.evaluate("document.fullscreenElement?.className") == "viewer"
    page.get_by_role("button", name="Exit fullscreen", exact=True).click()
    page.wait_for_function("document.fullscreenElement === null", timeout=10_000)


def test_tour_steps_drive_the_viewer(page, base_url):
    open_sample(page, base_url, "Si")
    assert page.locator(".tour-steps li").count() == 16
    zone = page.evaluate("window.__xf.tourStep('zone')")
    assert zone["tab"] == "symmetry"
    assert zone["zone"] == {"vertices": 24, "faces": 14, "edges": 36}
    measure = page.evaluate("window.__xf.tourStep('measure')")
    assert measure["brillouin"] is False
    assert measure["measure"]["kind"] == "angle (deg)"
    assert abs(measure["measure"]["value"] - 109.4712) < 0.01
    supercell = page.evaluate("window.__xf.tourStep('supercell')")
    assert supercell["summary"]["n_atoms"] == 64
    twist = page.evaluate("window.__xf.tourStep('twist')")
    assert twist["summary"]["n_atoms"] == 2 * 2 * 28
    for step in ["pxrd", "nanotube", "c60", "mof", "carbon", "polymer", "uff", "export"]:
        current = page.evaluate(f"window.__xf.tourStep('{step}')")
        assert current["scene"]["atoms"] > 0, step
    assert page.locator(".tour-steps li.active").count() == 0
    page.get_by_role("button", name="Tour", exact=True).click()
    assert page.locator(".tour-steps li.active").count() == 1
