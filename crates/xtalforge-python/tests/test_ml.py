"""Machine-learning bridge and model golden data (Tasks 8.4.1 and 8.4.3)."""

import json
import os
from pathlib import Path

import numpy as np
import pytest
import xtalforge

FIXTURE = Path(__file__).parent / "fixtures" / "ml" / "mace_mp0_small.json"


class FixedCalculator:
    """A duck-typed ASE calculator with constant results."""

    name = "fixed"

    def get_potential_energy(self, atoms):
        return -1.25

    def get_forces(self, atoms):
        return [[0.1, 0.0, 0.0], [-0.1, 0.0, 0.0]]

    def get_stress(self, atoms, *args, **kwargs):
        return np.array([0.001, 0.0, 0.0, 0.0, 0.0, 0.0])


def structure(cell, atoms):
    lines = [
        "data_ml",
        f"_cell_length_a {cell}",
        f"_cell_length_b {cell}",
        f"_cell_length_c {cell}",
        "_cell_angle_alpha 90",
        "_cell_angle_beta 90",
        "_cell_angle_gamma 90",
        "_symmetry_space_group_name_H-M 'P 1'",
        "loop_",
        "_atom_site_label",
        "_atom_site_type_symbol",
        "_atom_site_fract_x",
        "_atom_site_fract_y",
        "_atom_site_fract_z",
    ]
    for index, (element, point) in enumerate(atoms):
        fraction = [value / cell for value in point]
        lines.append(
            f"{element}{index + 1} {element} {fraction[0]:.12f} {fraction[1]:.12f} {fraction[2]:.12f}"
        )
    return xtalforge.Structure("\n".join(lines) + "\n", "cif")


def test_the_bridge_converts_energy_forces_and_stress() -> None:
    ml = xtalforge.MlCalculator(FixedCalculator())
    assert ml.name() == "fixed"
    document = structure(10.0, [("H", [0.0, 0.0, 0.0]), ("H", [0.74, 0.0, 0.0])])
    assert ml.energy(document.to_json()) == -1.25
    result = json.loads(ml.energy_forces(document.to_json()))
    assert result["energy_ev"] == -1.25
    assert result["forces"] == [[0.1, 0.0, 0.0], [-0.1, 0.0, 0.0]]
    assert result["n_atoms"] == 2
    assert result["stress_gpa"][0][0] == pytest.approx(0.001 * 160.21766208)


def test_the_backend_report_lists_the_modules() -> None:
    report = json.loads(xtalforge.ml_backends())
    for module in ["torch", "ase", "mace", "orb_models", "mattersim"]:
        assert module in report
        assert "available" in report[module]
    assert report["ase"]["available"] is True


@pytest.mark.skipif(
    not os.environ.get("XTALFORGE_MACE_MODEL"),
    reason="set XTALFORGE_MACE_MODEL to the cached MACE model file",
)
def test_mace_matches_the_recorded_golden_data() -> None:
    from mace.calculators import MACECalculator

    golden = json.loads(FIXTURE.read_text())
    calculator = MACECalculator(
        model_paths=os.environ["XTALFORGE_MACE_MODEL"],
        device="cpu",
        default_dtype="float64",
    )
    ml = xtalforge.MlCalculator(calculator)
    cases = {
        "h2": [("H", [0.0, 0.0, 0.0]), ("H", [0.74, 0.0, 0.0])],
        "water": [
            ("O", [0.0, 0.0, 0.0]),
            ("H", [0.96, 0.0, 0.0]),
            ("H", [-0.24, 0.93, 0.0]),
        ],
        "methane": [
            ("C", [0.0, 0.0, 0.0]),
            ("H", [0.63, 0.63, 0.63]),
            ("H", [0.63, -0.63, -0.63]),
            ("H", [-0.63, 0.63, -0.63]),
            ("H", [-0.63, -0.63, 0.63]),
        ],
    }
    for name, atoms in cases.items():
        document = structure(10.0, atoms)
        result = json.loads(ml.energy_forces(document.to_json()))
        expected = golden["cases"][name]
        assert result["energy_ev"] == pytest.approx(expected["energy_ev"], abs=1e-6)
        for ours, theirs in zip(result["forces"], expected["forces"]):
            for component, reference in zip(ours, theirs):
                assert component == pytest.approx(reference, abs=1e-5)
        scale = 160.21766208
        stress = expected["stress_ev_per_a3"]
        expected_tensor = [
            [stress[0] * scale, stress[5] * scale, stress[4] * scale],
            [stress[5] * scale, stress[1] * scale, stress[3] * scale],
            [stress[4] * scale, stress[3] * scale, stress[2] * scale],
        ]
        for ours_row, expected_row in zip(result["stress_gpa"], expected_tensor):
            for component, reference in zip(ours_row, expected_row):
                assert component == pytest.approx(reference, abs=1e-6)
