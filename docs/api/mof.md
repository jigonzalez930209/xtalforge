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
not a silent error. Unknown names are typed errors.

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

- Six topologies and eight blocks, not the full PORMAKE catalogue.
- No RMSD locator and no automatic SBU recognition; the linker/node choice is
  explicit in the request.
- Rods, layers, merged nets and decoration are the simplified helpers above,
  not full PORMAKE net manipulations.
- Blocks are not relaxed; the builder is geometric.
- The full MOF and RCSR bundles are external data (`data-YYYY.MM.DD`), not
  part of this crate.
