# End-to-end audit battery

This directory verifies the public surfaces against independent reference
implementations and against deterministic 3D frames.

- `reference/generate.py` reads the corpus with **gemmi**, detects space groups
  and primitive cells with **spglib**, derives the lattice, neighbour list and
  nanotube geometry with **ASE**, and computes Bragg angles with **NumPy**.
  It writes `reference/reference.json`, which is committed.
- `test_numeric.py` drives the WASM playground in Chromium and compares the
  numbers, P1 positions (modulo lattice translations), lattice corners,
  projected atom centres, supercells, slabs, symmetry transforms, builders
  and PXRD peak positions with that reference.
- `test_visual.py` compares deterministic canvas screenshots with the PNGs in
  `golden/`.
- `test_demos.py` smoke-tests the older cell-volume, CIF-inspector and editor
  demos.

## Requirements

The repository virtual environment needs the reference stack and Playwright:

```bash
.venv/bin/python -m pip install -r tests/e2e/requirements.txt
.venv/bin/python -m playwright install chromium   # only if Chrome is absent
```

## Run

```bash
scripts/e2e.sh
```

The script rebuilds the WASM package and the VitePress site, regenerates the
reference JSON, starts a preview server on a free port, and runs the tests.
Set `XF_UPDATE_GOLDENS=1` to accept new screenshots after reviewing the diff;
the previous goldens are overwritten, so review them with `git diff --stat`
and the image viewer before committing.

## Known gaps

The battery covers the structure, symmetry, builder, view and PXRD surfaces.
Force fields, porosity, MOF chemistry and the powder refinement stages still
rely on their Rust golden-data tests; extending the same scenario runner to
those surfaces is tracked in the roadmap. The CLI, Python and WASM smoke tests
run in their own crate test suites.
