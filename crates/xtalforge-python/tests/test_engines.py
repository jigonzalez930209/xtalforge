"""Simulation-engine Python surfaces (Task 18.5)."""

from pathlib import Path

import pytest

xtalforge = pytest.importorskip("xtalforge")

ROOT = Path(__file__).resolve().parents[3]
NACL = (ROOT / "tests" / "corpus" / "simple" / "NaCl.cif").read_text()
EIGENVAL = (
    ROOT
    / "crates"
    / "xtalforge-engines"
    / "tests"
    / "fixtures"
    / "engines"
    / "vasp"
    / "EIGENVAL"
).read_text()


def test_presets_and_failures():
    presets = xtalforge.engine_presets()
    assert len(presets) >= 20
    assert any(name == "vasp/relax" for name, _ in presets)
    codes = {code for code, _, _ in xtalforge.engine_failures()}
    assert "scf_not_converged" in codes


def test_input_set_and_parser():
    summary = xtalforge.engine_input_set(NACL, "vasp/relax")
    assert summary["preset"] == "vasp/relax"
    names = {item["name"] for item in summary["files"]}
    assert names == {"INCAR", "KPOINTS", "POSCAR", "POTCAR.spec"}
    parsed = xtalforge.engine_parse(EIGENVAL, "vasp-eigenval")
    assert parsed["bands"]["energies"]
    bands = xtalforge.engine_bands(EIGENVAL, "vasp-eigenval")
    assert bands["gap"]["energy"] >= 0.0
    script = xtalforge.engine_script("slurm", "xtalforge info a.cif")
    assert "#SBATCH" in script
