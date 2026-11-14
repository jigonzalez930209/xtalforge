# MOF builder

`xtalforge-mof` follows the PORMAKE idea (MIT; see `THIRD_PARTY.md`): a net
supplies the topology, a node block sits on every vertex and a linker spans
every edge. The catalogue here is a curated subset, not the full PORMAKE data.

## Catalogues

`topologies()` reads six nets from `data/topologies.json`: `pcu` (6-connected),
`dia` (4), `srs` (3), `sql` (4), `hcb` (3) and `bcu` (8). Each entry carries
fractional vertices and edge pairs.

`blocks()` reads four nodes (`Zn4O`, `Cu2` paddlewheel, `Cr3O` trimer,
`benzene_tri`) and four linkers (`benzene`, `biphenyl`, `stilbene`,
`triphenylene`) from `data/blocks.json`. A block carries its own cell, atoms,
connection points and elements. Coordinates are illustrative geometries, not
optimised structures.

## Build

```rust
use xtalforge::mof::{MofRequest, build_mof};

let built = build_mof(&MofRequest::new("pcu", "Zn4O", "benzene"))?;
println!("{} atoms, {} warnings", built.structure.sites().len(), built.warnings.len());
```

The builder fits a uniform cell scale to the linker, places the node at every
vertex, orients it so its connection points follow the incident edges, places
and scales the linker along each edge, and returns a P1 structure. `warnings`
reports a degree mismatch, a linker outside the 0.8-1.25 scale range, a linker
with more than two connection points and atom pairs closer than 0.9 Angstrom.
A vertex degree that does not match the node's connection count is a warning,
not a silent error. Unknown names are typed errors. Two metal-free blocks that
both carry oxygen (the carboxylate `benzene_tri` node with a carboxylate
linker) are refused with `invalid_input`, because their oxygens would overlap.

`MofRequest`, `BuiltMof { structure, topology, node, linker, warnings }`,
`BlockKind` (`Node`/`Linker`) and the `find_topology`/`find_block` lookups are
the public types. `build_request_json` accepts the same request as JSON;
`xtalforge mof build` exposes it in the CLI, `build_mof` in Python and
`build_mof_info` in WASM.

## Advanced builders

The crate also exposes documented geometric helpers:

- `linkers_by_span(min, max)` lists linker blocks by connection-point distance
  and `invent_linker(span, aromatic)` builds a linear carbon chain.
- `build_rod_framework(repeats)` places a metal chain with carbon linkers and
  `build_layered_framework(layers)` stacks honeycomb sheets.
- `merge_frameworks(first, second, shift)` merges two structures at a
  fractional shift; `decorate_super_tetrahedron(structure, node)` adds a
  tetrahedral O shell.
- `partition_pore_space(frameworks, spacing, probe)` reports each framework's
  accessible volume.
- `node_oriented(block, target_direction)` is the deterministic orientation
  helper used by the builder.

## Deliberate limits

- Six topologies and eight blocks are always embedded. The `mof_data` (`mof-data` in the facade) feature
  adds the full PORMAKE 0.2.3 bundle (2,129 quotient nets and 867 blocks) from
  two lazy gzip files; Python and WASM leave it off for package size, recorded
  in the parity exceptions.
- The global RMSD/Kabsch locator is implemented (`align`); the linker/node
  choice is explicit in the request.
- `database_fetch`: `mof::fetch_core_mof` / `fetch_qmof` download an archive into a cache
  directory, verify its SHA-256 and reuse it offline. The CLI is
  `xtalforge mof database fetch`; Python is `mof_fetch_database`; WASM uses
  `fetch_data_bundle` and keeps the bytes in JavaScript.
- Rods, layers, merged nets and decoration are the simplified helpers above,
  not full PORMAKE net manipulations.
- Blocks are not relaxed; the builder is geometric.
- The RCSR catalogue is the 2019-06-01 snapshot. The PORMAKE bundle is
  embedded only behind `mof-data`; its license and conversion are recorded in
  `THIRD_PARTY.md`.


## MOF platform (Phase 10)

### Cleaning and validation

- `mof_check(text, format_hint=None)`: overlaps, unbound atoms and fragments,
  impossible valences, hypervalent carbon, missing hydrogens, disorder, metal
  coordination and charge balance. Each finding carries a severity, the P1
  atoms and a `fix` verb. In Rust the function is `check_mof`; the CLI is
  `xtalforge mof check`.
- `mof_solvent_remove(text, remove_ions=False)`: removes free solvent from
  the embedded H2O/DMF/EtOH table and small organic fragments; counter-ions
  are kept unless `remove_ions` is set. Rust `remove_solvent`, CLI
  `mof solvent-remove`.
- `mof_disorder(text)`: picks the lowest-energy ordered decoration from the
  partial-occupancy groups. Rust `resolve_disorder`, CLI `mof disorder-resolve`.
- `mof_charge_balance(text)`: metal states from BVS, deprotonated carboxylate
  charge and the neutral/charged verdict. Rust `charge_balance`, CLI
  `mof charge-balance`.

### Deconstruction and identity

- `mof_deconstruct(text)`: metal nodes (all-node rule), organic linkers,
  solvent units, connection points, the contracted net and its RCSR names;
  `interpenetration` counts the framework components. Rust `deconstruct`,
  CLI `mof deconstruct`.
- `mof_mofid(text)`: canonical SMILES of the distinct nodes and linkers plus
  the topology and interpenetration in a deterministic `MOFid`; `mofkey` is
  the SHA-256 of that string. Rust `mofid`, CLI `mof mofid`.
- `mof_open_metal_sites(text)`: metals below their typical coordination with
  a free direction and a probe clearance. Rust `open_metal_sites`, CLI
  `mof open-metal-sites`.
- `mof_descriptors(text)`: formula, density, metals, node nuclearity and
  connectivity, linker span and functional groups. Rust
  `framework_descriptors`, CLI `mof descriptors`.

### Porosity

- `mof_pores(text, spacing=0.8, probe=1.2)`: native PLD by periodic
  union-find percolation, LCD, accessible surface area in m2/g and m2/cm3,
  accessible volume, void fraction and channel dimensionality, plus the
  geometric pore-size distribution. Rust `pore_descriptors`, CLI `mof pores`.

### Adsorption

- `mof_guests()`: the embedded TraPPE CO2, N2, CH4, H2 and H2O models.
- `mof_forcefield(text)`: UFF host parameters, EQeq framework charges and
  Lorentz-Berthelot mixing.
- `mof_widom(text, guest="CO2", temperature=298, samples=1000, seed=0)`:
  deterministic Widom insertion, dimensionless Henry coefficient, mean
  interaction energy and heat of adsorption.
- `mof_gcmc(text, guest="CO2", temperature=298, pressures=[1,5], steps=1000,
  seed=0)`: single-component isotherm with ideal-gas chemical potential.
- `mof_energy_grid(text, guest="CO2", spacing=1.0)`: guest energy grid
  summary for binding-site inspection.
- `mof_raspa(text, guest="CO2", temperature=298, pressure=1)`: RASPA
  `simulation.input`, CIF, force-field and mixing files.
- Rust names are `guests`, `force_field`, `widom_insertion`,
  `gcmc_isotherm`, `energy_grid` and `raspa_input`; CLI commands live under
  `mof adsorption`.

### Modification and variants

- `mof_functionalize(text, group="amino", positions=[], fraction=1, seed=0)`:
  -NH2, -OH, -CH3 and -F at chosen or all eligible linker carbons.
- `mof_defects(text, linker_fraction=0, cluster_fraction=0, cap="formate",
  seed=0)`: missing linkers and clusters with formate, water or hydroxide
  capping.
- `mof_interpenetrate(text, folds=2)`: n-fold interpenetration.
- `mof_stack_layers(text, layers=2, offset=(0.5,0.5))`: AA/AB layered
  variants.
- `mof_mixed(text, replacements, ordered=True, fraction=1, seed=0)`:
  ordered symmetry-orbit or seeded random mixed-metal/linker frameworks.
- Rust names are `functionalize`, `missing_defects`, `interpenetrate`,
  `stack_layers` and `mixed_framework`; CLI commands live under `mof modify`.

### Databases and screening

- `mof_read_database(path, kind="core")`: CoRE MOF or QMOF directory reader
  with property CSV support and license/citation metadata. Rust
  `read_core_mof` / `read_qmof`, CLI `mof database read`.
- `mof_screen(path, guest="CO2", spacing=0.8, samples=100, limit=0)`: cleans,
  deconstructs and computes descriptors, pores and Henry coefficients for a
  database. Rust `screen_database`, CLI `mof database screen`.
- `mof_deduplicate(structures)`: MOFkey first, structure matcher for ties.
  Rust `deduplicate`, CLI `mof database deduplicate`.
