#!/usr/bin/env python3
"""Generates the disordered-carbon golden data (Task 6.8.5).

Runs the Crystal-Builder recipe builder for the diamond net at 1x1x1,
0.42 g/cm3, 38 % coverage, one layer and the default edge ratios, with
relaxation disabled so the comparison is about the topology and the
termination. The RNG and the remesh differ from the Rust port, so the test
compares the distributions and the density, not the exact atom list.

    PYTHONPATH=/path/to/Crystal-Builder python scripts/golden/carbon.py
"""

import json
from pathlib import Path

from xtal.carbon.build import Recipe, build

TARGET = (
    Path(__file__).resolve().parents[2]
    / "tests/golden/phase-6/carbon/dia.json"
)


def main() -> None:
    recipe = Recipe(repeat=(1, 1, 1), relax="none", seed=0)
    built = build(recipe)
    data = {
        "net": recipe.net,
        "repeat": list(recipe.repeat),
        "density": recipe.density,
        "radius_ratio": recipe.radius_ratio,
        "coverage": recipe.coverage,
        "seed": recipe.seed,
        "atoms": len(built.structure.sites),
        "built_density": built.density,
        "ratios": built.ratios,
        "edge_fraction": built.edge_fraction,
        "pieces": built.pieces,
        "periodicity": built.periodicity,
        "census": {str(size): count for size, count in built.census.items()},
        "closed_carbons": built.closed_carbons,
        "edge_length": built.edge_length,
        "euler": built.euler,
        "gauss_bonnet": built.gauss_bonnet,
    }
    TARGET.parent.mkdir(parents=True, exist_ok=True)
    TARGET.write_text(json.dumps(data, indent=1) + "\n")
    print(f"wrote {TARGET}: {data['atoms']} atoms, {data['built_density']:.4f} g/cm3")


if __name__ == "__main__":
    main()
