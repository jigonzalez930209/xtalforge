# Generating structures and screening MOFs

This tutorial builds a small hypothetical-MOF screen in one workflow: enumerate
frameworks, relax the promising ones, describe them and keep the provenance.

## 1. Check the catalogue

```bash
xtalforge generate prototypes
xtalforge generate layer-groups | head
xtalforge generate rod-groups
```

The prototype list is the curated Strukturbericht/AFLAT label subset shipped
with the crate. The layer list is the 80 layer groups; the rod list is the
fourteen-group rod subset.

## 2. Build a reference

```bash
xtalforge generate prototype B1 --out nacl.cif
xtalforge generate identify nacl.cif --json
```

`identify` prints the label, signature and matcher RMS. Substituting an element
keeps the prototype:

```bash
xtalforge generate substitute nacl.cif --radius-tolerance 0.4 --json
```

## 3. Generate candidates

```bash
xtalforge generate crystal \
  --composition "Zn:2,O:2" --space-group 186 --seed 11 \
  --out zno.cif

xtalforge generate crystal-range \
  --groups 200-210 --composition "Na:4,Cl:4" --seed 5
```

The range command reports one JSON object per group with `status`, `attempt`
and `wyckoff`, or the rejection reason.

## 4. Relax and rank

```bash
xtalforge generate search \
  --composition "Na:4,Cl:4" --space-groups 225 \
  --structures 8 --seed 3 --steps 100 --report search.json

xtalforge generate evolve \
  --composition "Na:4,Cl:4" --space-groups 225 \
  --population 8 --generations 3 --seed 3
```

`search` prints the best energy per atom and writes the best structure with
`--out best.cif`; `evolve` adds the generation table to the report. Both use
UFF4MOF and the Phase 16 deduplicator.

## 5. Screen MOFs

```bash
xtalforge generate mof-dataset \
  --topologies pcu --limit 4 --seed 1 --no-relax \
  --out mof-dataset.json

xtalforge generate functionalize framework.cif \
  --groups amino,fluorine --fractions 0.5,1.0 \
  --out functionalization.json
```

The dataset records topology, node, linker, structure, optional energy,
framework descriptors, pore descriptors and MOFkey duplicates. A combination
that cannot be assembled is recorded with its error and the run continues.

## 6. Make the search resumable

```bash
xtalforge generate workflow \
  --kind random --seed-path seed.cif \
  --options '{"composition":{"Na":4,"Cl":4},"space_groups":[225],"n_structures":16,"seed":0}' \
  --checkpoint search.checkpoint.json --output search-results.json \
  > search.workflow.json

xtalforge run search.workflow.json
```

The workflow's `compute` step calls the agent `search` verb, which reads the
checkpoint, runs the next chunk and writes it back. The same operation is
available as `search_workflow` in Python and Rust.

## 7. From Python

```python
import xtalforge

crystal = xtalforge.generate_crystal({"Na": 4, "Cl": 4}, space_group=225, seed=7)
report = xtalforge.random_search({"Na": 4, "Cl": 4}, [225], structures=4, seed=1)
label = xtalforge.identify_prototype(open("nacl.cif").read())
print(label["label"], report["hits"][0]["energy_per_atom"])
```

## Limits to keep in mind

- A generated structure is a proposal; relax it before comparing energies.
- The prototype and rod tables are documented subsets.
- The generator caps a conventional cell at 400 atoms.
- WASM searches run synchronously; put them in a worker for long runs.
- Checkpoints need a filesystem, so they are Rust/CLI/Python only.
