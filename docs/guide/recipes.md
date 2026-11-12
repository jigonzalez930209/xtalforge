# Recipes: how do I …

Thirty task-oriented recipes. Each row names the operation and the entry point
on the four surfaces; the CLI column is runnable as written, and the Python,
WASM and Rust names are the exported symbols from the parity matrix. Format
hints are optional when the content is a CIF, POSCAR, XYZ or structure JSON.

| # | How do I … | CLI | Python | WASM | Rust |
|---|---|---|---|---|---|
| 1 | read a structure | `xtalforge info a.cif` | `xtalforge.read_structure(text)` | `read_structure(text, hint)` | `xtalforge::read(path)` |
| 2 | write POSCAR | `xtalforge convert a.cif --format poscar --out a.vasp` | `xtalforge.write(path, text)` | `convert_structure(text, "", "poscar")` | `xtalforge::write(path, &s)` |
| 3 | detect the space group | `xtalforge symmetry detect a.cif` | `detect_symmetry(text)` | `detect_symmetry(text, 0.001)` | `symmetry::detect_symmetry` |
| 4 | expand to P1 | `xtalforge symmetry p1 a.cif` | `to_p1(text)` | `to_p1_info(text)` | `symmetry::to_p1` |
| 5 | make a 2x2x2 supercell | `xtalforge supercell a.cif 2 2 2` | `supercell(text, 2, 2, 2)` | `supercell_info(text, 2, 2, 2)` | `Structure::supercell` |
| 6 | cut a (1 1 1) slab | `xtalforge slab a.cif 1 1 1 --layers 3` | `slab(text, 1, 1, 1, 3, 15)` | `slab_info(text, 1, 1, 1, 3, 15)` | `core::make_slab` |
| 7 | compute bonds and rings | `xtalforge bonds a.cif` | `bonds(text)` | `bonds(text)` | `Structure::bond_graph` |
| 8 | measure a distance | `xtalforge distance a.cif 0 1` | `distance(text, 0, 1)` | `distance(text, 0, 1)` | `core::distance` |
| 9 | compute PXRD | `xtalforge pxrd a.cif` | `pxrd(text)` | `pxrd_info(text)` | `analysis::calculate_pxrd` |
| 10 | run UFF energy | `xtalforge ff energy a.cif` | `ff_energy(text)` | `ff_energy(text)` | `UffCalculator::energy` |
| 11 | optimise with UFF | `xtalforge ff optimize a.cif` | `ff_optimize(text)` | `ff_optimize_text(text)` | `ff::minimize` |
| 12 | compute charges | `xtalforge ff charges a.cif --method eeq` | `ff_charges(text, "eeq")` | `ff_charges_info(text, "eeq")` | `ff::eeq_charges` |
| 13 | compare two structures | `xtalforge compare a.cif b.cif` | `match_structures(a, b)` | `match_structures(a, b, "", "{}")` | `symmetry::match_structures` |
| 14 | group a list | `xtalforge group a.cif b.cif c.cif` | `group_structures([...])` | `group_structures([...], "", "{}")` | `symmetry::group_structures` |
| 15 | apply a strain | `xtalforge deform a.cif --voigt 0.01,0,0,0,0,0` | `deform(text, voigt=[...])` | `deform(text, "", options)` | `core::strain` |
| 16 | rattle with a seed | `xtalforge features rattle a.cif --seed 7` | `rattle(text, seed=7)` | `rattle(text, "", "{\"seed\":7}")` | `core::rattle` |
| 17 | interpolate two frames | `xtalforge features interpolate a.cif b.cif --images 5` | `interpolate(a, b, 5)` | `interpolate(a, b, 5, "", "{}")` | `symmetry::interpolate_matched` |
| 18 | make a cubic supercell | `xtalforge features optimal-supercell a.cif 8` | `optimal_supercell(text, 8)` | `optimal_supercell(text, 8, "")` | `core::optimal_supercell` |
| 19 | compute Voronoi cells | `xtalforge features voronoi a.cif` | `voronoi(text)` | `voronoi(text, "", 6.0)` | `core::voronoi` |
| 20 | get a coordination shell | `xtalforge features coordination a.cif 0` | `coordination(text, 0)` | `coordination(text, "", 0, "voronoi", 5.0, 0.1)` | `core::coordination` |
| 21 | compute oxidation states | `xtalforge features oxidation a.cif` | `oxidation_states(text)` | `oxidation_states(text, "", 5.0)` | `analysis::oxidation_states` |
| 22 | compute an RDF | `xtalforge features rdf a.cif` | `radial_distribution(text)` | `radial_distribution(text, "", 10, 0.05, 0)` | `core::radial_distribution` |
| 23 | list unique Miller indices | `xtalforge features miller a.cif` | `unique_miller_indices(text)` | `unique_miller_indices(text, "", "{}")` | `analysis::unique_miller_indices` |
| 24 | build a Wulff shape | `xtalforge features wulff a.cif --surface 1,0,0,1.0` | `wulff(text, [[1,0,0,1.0]])` | `wulff(text, "", options)` | `core::wulff` |
| 25 | find adsorption sites | `xtalforge features adsorption a.cif` | `adsorption_sites(text)` | `adsorption_sites(text, "", 0.05)` | `core::adsorption_sites` |
| 26 | check a MOF | `xtalforge mof check a.cif` | `mof_check(text)` | `mof_check(text, "")` | `mof::check_mof` |
| 27 | compute MOFid | `xtalforge mof mofid a.cif` | `mof_mofid(text)` | `mof_mofid(text, "")` | `mof::mofid` |
| 28 | compute PLD/LCD/ASA | `xtalforge mof pores a.cif` | `mof_pores(text)` | `mof_pores(text, "", "{}")` | `mof::pore_descriptors` |
| 29 | run Widom insertion | `xtalforge mof adsorption widom a.cif --guest CO2` | `mof_widom(text)` | `mof_widom(text, "", "{}")` | `mof::widom_insertion` |
| 30 | build a polymer | `xtalforge build polymer c1ccc(cc1)C=C 10` | `build_polymer("c1ccc(cc1)C=C", 10)` | `build_polymer_text(...)` | `build::build_polymer` |

## Pipelining

The CLI reads standard input with `-` and writes a structure document to
standard output with `--out -`:

```bash
xtalforge supercell a.cif 2 2 2 --out - \
  | xtalforge slab - 1 1 1 --layers 3 --out slab.cif
```

The input format is sniffed from the content, so no temporary files or format
flags are needed.

## Keyword conventions

The naming table in [`docs/api/naming.md`](../api/naming.md) lists the
canonical option names on every surface. When a name differs, the 1.x alias is
kept and documented; see [`deprecation.md`](./deprecation.md).
