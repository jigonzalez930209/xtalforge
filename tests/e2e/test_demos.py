"""Smoke tests for the older documentation demos."""

from __future__ import annotations


def test_cell_volume_demo_runs_wasm(page, base_url):
    page.goto(f"{base_url}/demos/cell-volume", wait_until="load")
    page.wait_for_function(
        "document.querySelector('.status')?.textContent?.includes('WASM')",
        timeout=60_000,
    )
    page.wait_for_function(
        "document.querySelector('.result')?.textContent?.includes('48.827')",
        timeout=10_000,
    )


def test_cif_inspector_reads_silicon(page, base_url):
    page.goto(f"{base_url}/demos/cif-inspector", wait_until="load")
    page.wait_for_selector(".summary", timeout=60_000)
    summary = page.locator(".summary h3").inner_text()
    assert "Si8" in summary
    assert page.locator("details table tbody tr").count() > 0


def test_document_editor_undoes_a_move(page, base_url):
    page.goto(f"{base_url}/demos/editor", wait_until="load")
    page.wait_for_selector(".summary", timeout=60_000)
    page.get_by_role("button", name="Move atom 0").click()
    page.wait_for_function(
        "document.querySelector('.row button')?.disabled === false", timeout=10_000
    )
    page.get_by_role("button", name="Undo").click()
    assert page.locator(".summary").is_visible()
