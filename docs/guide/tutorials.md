# Tutorials

Eight end-to-end workflows. Commands use the shipped corpus; replace the paths
with your files. Every step is also available as an agent verb; the
[automation guide](./automation) shows the JSON form.

## 1. CIF to supercell, bonds and symmetry

```bash
xtalforge info tests/corpus/simple/NaCl.cif
xtalforge supercell tests/corpus/simple/NaCl.cif 2 2 2 --out nacl-222.cif
xtalforge bonds nacl-222.cif --json
xtalforge symmetry detect nacl-222.cif --symprec 1e-5 --json
```

`info` reports the formula, cell, density and charge. `supercell` writes the
expanded cell (the CLI works on the P1 expansion). `bonds` returns the
perceived graph with lattice images, and `symmetry detect` the space group,
Wyckoff letters and standardized cell. Check that the space group survived the
supercell; it must be the same group with the same origin.

## 2. Prepare a deposited CIF

```bash
xtalforge ff prepare tests/corpus/cod/MOF-5.cif --json
xtalforge info tests/corpus/cod/MOF-5.cif
xtalforge symmetry detect tests/corpus/cod/MOF-5.cif --symprec 0.001
```

The preparation pipeline relaxes the geometry, assigns charges, runs Ewald and
reports every step. Read the report before using the file: a large relaxation
or a charge warning means the deposited model needs a decision. `symprec 0.001`
is deliberately loose; report the value you used.

The interactive equivalent is `load`, `info`, `diagnostics`, `hydrogens`,
`wrap`, `symmetry` in the agent session.

## 3. UFF optimization and a relaxed scan

```bash
xtalforge ff energy water.xyz
xtalforge ff optimize water.xyz --max-steps 200 --out water-relaxed.xyz
xtalforge ff vibrations water-relaxed.xyz --json
```

`energy` reports eV and forces. `optimize` uses FIRE by default and writes the
relaxed geometry; check `converged` and `max_force`. `vibrations` gives the
Gamma-point modes and thermochemistry; imaginary modes below the documented
threshold are reported, not hidden.

A relaxed scan is a Rust or Python call (`xtalforge_ff::relaxed_scan`), because
it needs a constraint object rather than a CLI flag.

## 4. Porosity report

```bash
xtalforge analysis porosity tests/corpus/cod/MOF-5.cif --spacing 0.4 --probe 1.2 --json
```

Report the accessible volume, fraction, surface area and largest included
sphere diameter, and state the spacing and probe radius. Run a second spacing
(0.3) when the number matters; if the two differ by more than a few percent,
the grid is not converged.

## 5. PXRD against a measured pattern

```bash
xtalforge pxrd tests/corpus/simple/NaCl.cif --radiation cu-ka --max-two-theta 80 --out nacl.xy
xtalforge powder automatic measured.xy --json
```

The first command writes a two-column calculated pattern. The automatic run
finds peaks, indexes them, runs Le Bail/Pawley and ranks the leading cells and
extinction classes. Compare positions first, then intensities; a zero shift or
a different instrument profile changes the latter.

## 6. Powder refinement

```bash
xtalforge powder index measured.xy --json
xtalforge powder lebail nacl.cif measured.xy --json
xtalforge powder pawley nacl.cif measured.xy --json
xtalforge powder rietveld nacl.cif measured.xy --refine --json
```

Indexing proposes cells; Le Bail extracts intensities without a structural
model; Pawley refines the cell and profile; Rietveld refines the structure.
Report Rwp, the cell esds and every fixed parameter. The automatic command
does the sequence and applies the acceptance thresholds.

## 7. Build a MOF

```bash
xtalforge mof topologies
xtalforge mof blocks
xtalforge mof build pcu Zn4O benzene --out mof.cif
xtalforge analysis porosity mof.cif --spacing 0.5
xtalforge render mof.cif --out mof.png --format png
```

The crate ships a curated topology/block subset, not the full PORMAKE
database. The builder places blocks; it does not relax them. Check `info` and
the diagnostics, then measure. Do not assign a topology from the builder name
alone: run net identification and report ambiguity.

## 8. Build a carbon

```bash
xtalforge build graphene 4 4 --out graphene.cif
xtalforge build nanotube 10 10 --length 4 --out tube.cif
xtalforge build c60 --out c60.cif
xtalforge build carbon --net dia --repeat-x 2 --repeat-y 2 --repeat-z 2 --seed 7 --out carbon.cif
xtalforge info carbon.cif
```

The disordered-carbon recipe is seeded and deterministic. It reports density,
ring census, mean C-C and periodicity in `--json` mode. Run `info` and
diagnostics afterwards: an unconsumed edge appears as a dangling atom.

## Where to go next

- [Agent API](../api/agent): the same work as JSON verbs.
- [Algorithms](../algorithms/): what each number means and its tolerances.
- [Migration guide](./migration): the Crystal-Builder equivalent of each step.
