"""Benchmark cases shared by the comparative suite (Task 15.1.1).

Every case returns a tiny, JSON-serialisable output summary so the result file
stays deterministic while still proving that the timed call did real work.
Reference libraries are imported lazily; a missing package removes its cases
from the run instead of failing the suite.
"""

from __future__ import annotations

import gc
import io
import os
import re
import statistics
import subprocess
import time
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
CIF_READ = ROOT / "tests" / "corpus" / "prepared" / "UiO-66.cif"
MOF5 = ROOT / "tests" / "corpus" / "cod" / "MOF-5.cif"
SI = ROOT / "tests" / "corpus" / "simple" / "Si.cif"

CIF_READ_TEXT = CIF_READ.read_text(encoding="utf-8")
MOF5_TEXT = MOF5.read_text(encoding="utf-8")
SI_TEXT = SI.read_text(encoding="utf-8")

IMPORT_ERRORS: dict[str, str] = {}

try:
    import xtalforge as xf
except Exception as exc:  # pragma: no cover - reported by the CLI
    xf = None
    IMPORT_ERRORS["xtalforge"] = str(exc)

try:
    import gemmi
except Exception as exc:  # pragma: no cover
    gemmi = None
    IMPORT_ERRORS["gemmi"] = str(exc)

try:
    import spglib
except Exception as exc:  # pragma: no cover
    spglib = None
    IMPORT_ERRORS["spglib"] = str(exc)

try:
    from ase.io.cif import read_cif as ase_read_cif
    from ase.neighborlist import neighbor_list
except Exception as exc:  # pragma: no cover
    ase_read_cif = None
    neighbor_list = None
    IMPORT_ERRORS["ase"] = str(exc)

try:
    from pymatgen.analysis.diffraction.xrd import XRDCalculator
    from pymatgen.core import Structure as PmgStructure
    from pymatgen.io.cif import CifParser
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
except Exception as exc:  # pragma: no cover
    XRDCalculator = None
    PmgStructure = None
    CifParser = None
    SpacegroupAnalyzer = None
    IMPORT_ERRORS["pymatgen"] = str(exc)


TASK_INFO: dict[str, dict] = {
    "cif_read": {
        "title": "CIF read",
        "input": "tests/corpus/prepared/UiO-66.cif (114 P1 atoms)",
        "description": "Parse a P1 CIF into a usable structure object.",
        "timed": "file text in memory to a parsed structure; the atom count is returned",
        "references": ["gemmi", "ase", "pymatgen"],
        "unavailable": {"spglib": "spglib has no structure-file reader"},
    },
    "p1": {
        "title": "P1 expansion",
        "input": "tests/corpus/cod/MOF-5.cif (7 asymmetric sites)",
        "description": "Apply the space-group operations and return the 424-atom P1 cell.",
        "timed": "CIF text to the P1 cell; parsers are inside the timed region for every library",
        "references": ["gemmi", "ase", "pymatgen"],
        "unavailable": {"spglib": "spglib does not expand asymmetric units"},
    },
    "symmetry": {
        "title": "Symmetry detection",
        "input": "tests/corpus/cod/MOF-5.cif",
        "description": "Find the space-group number from the atomic coordinates.",
        "timed": "CIF text to the detected international number at symprec 0.001 A",
        "references": ["spglib", "pymatgen"],
        "unavailable": {
            "gemmi": "gemmi has no atom-position space-group finder",
            "ase": "ASE has no symmetry detector",
        },
    },
    "neighbors": {
        "title": "Periodic neighbors",
        "input": "tests/corpus/cod/MOF-5.cif, cutoff 3.0 A",
        "description": "List periodic atom pairs inside a cutoff.",
        "timed": "CIF text to the complete neighbor list",
        "references": ["ase", "pymatgen"],
        "unavailable": {
            "gemmi": "the gemmi NeighborSearch API targets macromolecular models",
            "spglib": "spglib has no neighbor search",
        },
    },
    "pxrd": {
        "title": "PXRD pattern",
        "input": "tests/corpus/cod/MOF-5.cif, Cu K-alpha, 2-theta 5-80 deg",
        "description": "Calculate the powder diffraction peak list.",
        "timed": "CIF text to the peak list",
        "references": ["pymatgen"],
        "unavailable": {
            "gemmi": "gemmi has no PXRD calculator",
            "ase": "ASE has no PXRD calculator",
            "spglib": "spglib has no PXRD calculator",
        },
    },
    "supercell": {
        "title": "Supercell",
        "input": "tests/corpus/simple/Si.cif, 50 x 50 x 50 (1,000,000 atoms)",
        "description": "Expand the eight-atom cubic cell into one million atoms.",
        "timed": "in-memory build only; xtalforge uses the core benchmark, references use their Python supercell operator",
        "references": ["ase", "pymatgen"],
        "unavailable": {
            "gemmi": "gemmi has no supercell builder",
            "spglib": "spglib has no supercell builder",
        },
    },
    "porosity": {
        "title": "Porosity",
        "input": "tests/corpus/cod/MOF-5.cif, 0.8 A grid, 1.2 A probe",
        "description": "Compute geometric pore volume, accessible volume and surface area.",
        "timed": "CIF text to the porosity report",
        "references": [],
        "unavailable": {
            "gemmi": "no geometric porosity implementation",
            "spglib": "no geometric porosity implementation",
            "ase": "no geometric porosity implementation",
            "pymatgen": "no geometric porosity implementation",
        },
    },
}

# Repetitions are deliberately small for the expensive cases. Warm-up calls are
# outside the timed samples and are never reported.
REPETITIONS = {
    "cif_read": 15,
    "p1": 9,
    "symmetry": 9,
    "neighbors": 9,
    "pxrd": 5,
    "supercell": 5,
    "porosity": 5,
}


@dataclass(frozen=True)
class Case:
    task: str
    library: str
    run: Callable[..., dict]
    repetitions: int
    warmups: int = 1
    external: bool = False


def _pymatgen_neighbors() -> dict:
    structure = _pmg_structure(MOF5_TEXT)
    return {
        "pairs": sum(
            len(structure.get_neighbors(site, 3.0)) for site in structure
        )
    }


def _gemmi_structure(text: str):
    """Parses CIF text with gemmi and returns its P1 site list."""
    document = gemmi.cif.read_string(text)
    structure = gemmi.make_small_structure_from_block(document.sole_block())
    return structure, structure.get_all_unit_cell_sites()


def _ase_atoms(text: str):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return ase_read_cif(io.StringIO(text))


def _pmg_structure(text: str):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return CifParser.from_str(text).parse_structures(primitive=False)[0]


def _symmetry_spglib() -> dict:
    structure, sites = _gemmi_structure(MOF5_TEXT)
    lattice = np.asarray(structure.cell.orth.mat.tolist(), dtype=float).T
    positions = np.asarray([tuple(site.fract) for site in sites], dtype=float)
    numbers = np.asarray(
        [site.element.atomic_number for site in sites], dtype=np.int32
    )
    dataset = spglib.get_symmetry_dataset(
        (lattice, positions, numbers), symprec=1e-3
    )
    return {"space_group": int(dataset.number)}


def _symmetry_pymatgen() -> dict:
    structure = PmgStructure.from_file(str(MOF5))
    analyzer = SpacegroupAnalyzer(structure, symprec=1e-3)
    return {"space_group": int(analyzer.get_space_group_number())}


def _cif_read_xtalforge() -> dict:
    return {"atoms": xf.read_cif(CIF_READ_TEXT)["n_atoms"]}


def _p1_xtalforge() -> dict:
    return {"atoms": xf.to_p1(MOF5_TEXT)["n_atoms"]}


def _symmetry_xtalforge() -> dict:
    return {"space_group": xf.detect_symmetry(MOF5_TEXT)["number"]}


def _neighbors_xtalforge() -> dict:
    return {"pairs": len(xf.neighbors(MOF5_TEXT, 3.0))}


def _pxrd_xtalforge() -> dict:
    return {"reflections": len(xf.pxrd(MOF5_TEXT)["reflections"])}


def _porosity_xtalforge() -> dict:
    report = xf.porosity(MOF5_TEXT, spacing=0.8, probe=1.2)
    return {"accessible_volume": round(float(report["accessible_volume"]), 6)}


def _cif_read_gemmi() -> dict:
    return {"atoms": len(_gemmi_structure(CIF_READ_TEXT)[1])}


def _p1_gemmi() -> dict:
    return {"atoms": len(_gemmi_structure(MOF5_TEXT)[1])}


def _cif_read_ase() -> dict:
    return {"atoms": len(_ase_atoms(CIF_READ_TEXT))}


def _p1_ase() -> dict:
    return {"atoms": len(_ase_atoms(MOF5_TEXT))}


def _neighbors_ase() -> dict:
    atoms = _ase_atoms(MOF5_TEXT)
    return {
        "pairs": len(
            neighbor_list("ijd", atoms, 3.0, self_interaction=False)[0]
        )
    }


def _supercell_ase() -> dict:
    return {"atoms": len(_ase_atoms(SI_TEXT) * (50, 50, 50))}


def _cif_read_pymatgen() -> dict:
    return {"atoms": len(_pmg_structure(CIF_READ_TEXT))}


def _p1_pymatgen() -> dict:
    return {"atoms": len(_pmg_structure(MOF5_TEXT))}


def _pxrd_pymatgen() -> dict:
    structure = _pmg_structure(MOF5_TEXT)
    pattern = XRDCalculator(wavelength=1.5406).get_pattern(
        structure, two_theta_range=(5, 80)
    )
    return {"peaks": len(pattern.x)}


def _supercell_pymatgen() -> dict:
    return {"atoms": len(_pmg_structure(SI_TEXT) * (50, 50, 50))}


def _supercell_xtalforge(repetitions: int = 3) -> dict:
    command = [
        "cargo",
        "bench",
        "-q",
        "-p",
        "xtalforge-core",
        "--bench",
        "supercell",
    ]
    environment = dict(os.environ)
    environment["XTALFORGE_BENCH_REPS"] = str(repetitions)
    completed = subprocess.run(
        command,
        cwd=ROOT,
        env=environment,
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(completed.stderr.strip() or "supercell benchmark failed")
    match = re.search(
        r"supercell atoms=(\d+) median_ms=([0-9.]+)", completed.stdout
    )
    samples = re.search(r"samples_ms=([0-9.,]+)", completed.stdout)
    if match is None or samples is None:
        raise RuntimeError("could not parse the supercell benchmark output")
    values = sorted(float(value) for value in samples.group(1).split(","))
    return {
        "median_ms": float(match.group(2)),
        "min_ms": values[0],
        "max_ms": values[-1],
        "samples_ms": values,
        "output": {"atoms": int(match.group(1))},
    }


def cases() -> list[Case]:
    """Builds the case table for the libraries that imported successfully."""
    result: list[Case] = []
    reps = REPETITIONS

    if xf is not None:
        result.append(
            Case("cif_read", "xtalforge", _cif_read_xtalforge, reps["cif_read"])
        )
        result.append(Case("p1", "xtalforge", _p1_xtalforge, reps["p1"]))
        result.append(
            Case("symmetry", "xtalforge", _symmetry_xtalforge, reps["symmetry"])
        )
        result.append(
            Case("neighbors", "xtalforge", _neighbors_xtalforge, reps["neighbors"])
        )
        result.append(Case("pxrd", "xtalforge", _pxrd_xtalforge, reps["pxrd"]))
        result.append(
            Case(
                "supercell",
                "xtalforge",
                _supercell_xtalforge,
                reps["supercell"],
                external=True,
            )
        )
        result.append(
            Case("porosity", "xtalforge", _porosity_xtalforge, reps["porosity"])
        )

    if gemmi is not None:
        result.append(
            Case("cif_read", "gemmi", _cif_read_gemmi, reps["cif_read"])
        )
        result.append(Case("p1", "gemmi", _p1_gemmi, reps["p1"]))

    if spglib is not None:
        result.append(
            Case("symmetry", "spglib", _symmetry_spglib, reps["symmetry"])
        )

    if ase_read_cif is not None:
        result.append(Case("cif_read", "ase", _cif_read_ase, reps["cif_read"]))
        result.append(Case("p1", "ase", _p1_ase, reps["p1"]))
        result.append(
            Case("neighbors", "ase", _neighbors_ase, reps["neighbors"])
        )
        result.append(Case("supercell", "ase", _supercell_ase, reps["supercell"]))

    if PmgStructure is not None:
        result.append(
            Case("cif_read", "pymatgen", _cif_read_pymatgen, reps["cif_read"])
        )
        result.append(Case("p1", "pymatgen", _p1_pymatgen, reps["p1"]))
        result.append(
            Case("symmetry", "pymatgen", _symmetry_pymatgen, reps["symmetry"])
        )
        result.append(
            Case("neighbors", "pymatgen", _pymatgen_neighbors, reps["neighbors"])
        )
        result.append(Case("pxrd", "pymatgen", _pxrd_pymatgen, reps["pxrd"]))
        result.append(
            Case(
                "supercell",
                "pymatgen",
                _supercell_pymatgen,
                reps["supercell"],
            )
        )

    return result


def _measure(case: Case) -> dict:
    if case.external:
        return case.run(case.repetitions)
    # A full collection immediately before a timed call makes the next
    # allocation pay for the collector's new generation, which inflated small
    # cases by up to 4x in early runs. Collect once, then keep the collector
    # out of the timed region.
    gc.collect()
    was_enabled = gc.isenabled()
    gc.disable()
    try:
        for _ in range(case.warmups):
            case.run()
        samples: list[float] = []
        output: dict = {}
        for _ in range(case.repetitions):
            start = time.perf_counter_ns()
            output = case.run()
            samples.append((time.perf_counter_ns() - start) / 1e6)
    finally:
        if was_enabled:
            gc.enable()
    return {
        "median_ms": round(statistics.median(samples), 4),
        "min_ms": round(min(samples), 4),
        "max_ms": round(max(samples), 4),
        "samples_ms": [round(value, 4) for value in samples],
        "output": output,
    }
