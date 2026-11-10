"""Shared helpers for the xtalforge end-to-end battery."""

from __future__ import annotations

import json
import math
import os
import re
from pathlib import Path

import numpy as np
from PIL import Image
from scipy.optimize import linear_sum_assignment

ROOT = Path(__file__).resolve().parents[2]
REFERENCE_PATH = ROOT / "tests" / "e2e" / "reference" / "reference.json"
GOLDEN = ROOT / "tests" / "e2e" / "golden"
ARTIFACTS = ROOT / "tests" / "e2e" / "test-results"
FOV = math.radians(45.0)
AZIMUTH = 0.85
ELEVATION = 0.42
PIXEL_TOLERANCE = 12
MISMATCH_LIMIT = 0.01

with REFERENCE_PATH.open() as handle:
    REFERENCE = json.load(handle)


def open_sample(page, base: str, sample: str) -> dict:
    page.goto(f"{base}/playground?sample={sample}", wait_until="load")
    page.wait_for_function("window.__xf && window.__xf.ready === true", timeout=60_000)
    return state(page)


def state(page) -> dict:
    return page.evaluate("window.__xf.state()")


def operation(page, name: str, **params) -> dict:
    return page.evaluate(
        "([name, params]) => window.__xf.operation(name, params)",
        [name, params],
    )


def pxrd(page) -> dict:
    return page.evaluate("window.__xf.pxrd()")


def canvas_metrics(page) -> tuple[int, int]:
    return tuple(
        page.evaluate(
            "(() => { const c = document.querySelector("
            "'[data-testid=playground-canvas]'); return [c.width, c.height]; })()"
        )
    )


def formula_counts(formula: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    for element, number in re.findall(r"([A-Z][a-z]?)([0-9]*\.?[0-9]+)?", formula):
        if not element:
            continue
        counts[element] = counts.get(element, 0) + round(float(number or 1))
    return counts


def cell_corners(matrix: list[list[float]]) -> np.ndarray:
    m = np.asarray(matrix, dtype=float)
    return np.array(
        [x * m[0] + y * m[1] + z * m[2] for x in (0, 1) for y in (0, 1) for z in (0, 1)]
    )


def origin_vectors(state_data: dict) -> list[np.ndarray]:
    vectors: list[np.ndarray] = []
    origin = np.zeros(3)
    for start, end in state_data["cellEdges"]:
        if np.linalg.norm(np.asarray(start) - origin) > 1e-6:
            continue
        vector = np.asarray(end, dtype=float)
        if all(np.linalg.norm(vector - known) > 1e-6 for known in vectors):
            vectors.append(vector)
    return vectors


def corners_match(state_data: dict, matrix: list[list[float]], tolerance=2e-3) -> float:
    page_points = np.unique(
        np.round(np.asarray(state_data["cellEdges"], dtype=float).reshape(-1, 3), 4),
        axis=0,
    )
    reference = np.round(cell_corners(matrix), 4)
    worst = 0.0
    for point in page_points:
        distance = min(float(np.linalg.norm(point - candidate)) for candidate in reference)
        worst = max(worst, distance)
    return worst


def periodic_match(
    page_atoms: list[dict], reference_atoms: list[list], matrix: list[list[float]]
) -> tuple[list[dict], float]:
    """Match P1 atoms to the gemmi expansion modulo lattice translations."""
    m = np.asarray(matrix, dtype=float)
    shifts = np.array(
        [[i, j, k] for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)]
    )
    expected: list[dict] = [{} for _ in page_atoms]
    worst = 0.0
    elements = {atom["element"] for atom in page_atoms}
    for element in elements:
        page_indices = [
            index
            for index, atom in enumerate(page_atoms)
            if atom["element"] == element
        ]
        reference_indices = [
            index
            for index, atom in enumerate(reference_atoms)
            if atom[0] == element
        ]
        if len(page_indices) != len(reference_indices):
            raise AssertionError(
                f"{element}: {len(page_indices)} page atoms vs "
                f"{len(reference_indices)} reference atoms"
            )
        page_positions = np.array(
            [page_atoms[index]["position"] for index in page_indices], dtype=float
        )
        reference_fractional = np.array(
            [reference_atoms[index][1:] for index in reference_indices], dtype=float
        )
        cost = np.zeros((len(page_positions), len(reference_fractional)))
        best_shifts = np.zeros((len(page_positions), len(reference_fractional), 3))
        for i, position in enumerate(page_positions):
            for j, fractional in enumerate(reference_fractional):
                candidates = (fractional + shifts) @ m
                distances = np.linalg.norm(candidates - position, axis=1)
                best = int(np.argmin(distances))
                cost[i, j] = distances[best]
                best_shifts[i, j] = shifts[best]
        rows, columns = linear_sum_assignment(cost)
        for row, column in zip(rows, columns):
            index = page_indices[row]
            adjusted = (reference_fractional[column] + best_shifts[row, column]) @ m
            expected[index] = {"element": element, "position": adjusted}
            worst = max(worst, float(cost[row, column]))
    return expected, worst


def project_positions(
    positions: np.ndarray,
    corners: np.ndarray,
    width: int,
    height: int,
) -> np.ndarray:
    points = np.vstack([positions, corners])
    low = points.min(axis=0)
    high = points.max(axis=0)
    center = (low + high) / 2.0
    radius = float(np.linalg.norm(high - low) / 2.0)
    distance = radius / math.sin(FOV / 2.0) * 1.15
    direction = np.array(
        [
            math.cos(ELEVATION) * math.cos(AZIMUTH),
            math.cos(ELEVATION) * math.sin(AZIMUTH),
            math.sin(ELEVATION),
        ]
    )
    camera = center + direction * distance
    forward = center - camera
    forward /= np.linalg.norm(forward)
    right = np.cross(forward, np.array([0.0, 0.0, 1.0]))
    right /= np.linalg.norm(right)
    up = np.cross(right, forward)
    relative = positions - camera
    x_camera = relative @ right
    y_camera = relative @ up
    z_camera = relative @ forward
    tangent = math.tan(FOV / 2.0)
    aspect = width / height
    x = (x_camera / (z_camera * tangent * aspect) + 1.0) / 2.0
    y = (1.0 - y_camera / (z_camera * tangent)) / 2.0
    return np.column_stack([x, y])


def projection_error(page, state_data: dict, sample: str) -> float:
    reference = REFERENCE["samples"][sample]
    expected, _ = periodic_match(
        state_data["atomPositions"], reference["atom_positions"], reference["cell_matrix"]
    )
    positions = np.array([atom["position"] for atom in expected], dtype=float)
    width, height = canvas_metrics(page)
    projected = project_positions(
        positions, cell_corners(reference["cell_matrix"]), width, height
    )
    page_projected = np.array(
        [[atom["x"], atom["y"]] for atom in state_data["projected"]], dtype=float
    )
    return float(np.abs(projected - page_projected).max())


def compare_image(actual: Path, name: str, update: bool | None = None) -> float:
    if update is None:
        update = os.environ.get("XF_UPDATE_GOLDENS") == "1"
    GOLDEN.mkdir(parents=True, exist_ok=True)
    golden = GOLDEN / f"{name}.png"
    if update and (not golden.exists() or actual.read_bytes() != golden.read_bytes()):
        golden.write_bytes(actual.read_bytes())
    if not golden.exists():
        raise AssertionError(f"missing golden {golden}; run with XF_UPDATE_GOLDENS=1")
    first = np.asarray(Image.open(actual).convert("RGB"), dtype=int)
    second = np.asarray(Image.open(golden).convert("RGB"), dtype=int)
    if first.shape != second.shape:
        raise AssertionError(f"{name}: screenshot shape {first.shape} != {second.shape}")
    difference = np.abs(first - second).max(axis=2)
    return float((difference > PIXEL_TOLERANCE).mean())


def screenshot(page, name: str) -> Path:
    ARTIFACTS.mkdir(parents=True, exist_ok=True)
    page.evaluate("window.__xf.setVisibility({cell: true, axes: true, rotate: false})")
    page.evaluate("window.__xf.setView({azimuth: 0.85, elevation: 0.42, zoom: 1})")
    page.wait_for_timeout(100)
    path = ARTIFACTS / f"{name}.png"
    page.locator("[data-testid=playground-canvas]").screenshot(path=str(path))
    return path


def reciprocal_spacing(matrix: list[list[float]], hkl: tuple[int, int, int]) -> float:
    m = np.asarray(matrix, dtype=float)
    g = np.asarray(hkl, dtype=float) @ np.linalg.inv(m).T
    return 1.0 / float(np.linalg.norm(g))


def surface_area(state_data: dict) -> float:
    vectors = origin_vectors(state_data)
    if len(vectors) < 2:
        raise AssertionError("the scene does not expose two cell vectors")
    return float(np.linalg.norm(np.cross(vectors[0], vectors[1])))
