# Molecule and carbon builders

`xtalforge-build` builds molecules, carbon structures, polymers and MOF
frameworks. Every builder returns a typed error and a deterministic structure.

## Carbon

```rust
use xtalforge_build::{fullerene, graphene, nanotube};

let sheet = graphene(4, 4, 1, 3.35, 15.0)?;
let tube = nanotube(5, 5, 2, 10.0)?;
let cage = fullerene(60)?;
```

`graphene(nx, ny, layers, interlayer, vacuum)` builds an AA-stacked honeycomb
sheet in a hexagonal cell (C-C = 1.42 Angstrom). `nanotube(n, m, length,
vacuum)` rolls the honeycomb along the chiral vector using a direct lattice
scan bucketed into the tube cell; the armchair and zigzag cases are validated
by atom counts and cell dimensions. `fullerene(60)` builds C60 from the
embedded ASE-derived coordinates; `face_spiral(20 | 60)` and
`spiral_sequence(size)` cover the C20 dodecahedron and the C60 cage.

`build_net(&CgdEntry, "C")` turns a CGD net into a P1 structure: vertices
become sites and edges become topology bonds, kept separate from chemistry.

## SMILES, embedding and sketches

`parse_smiles` returns a `SmilesMolecule` with `SmilesAtom` and `SmilesBond`
entries: atoms, aromatic atoms, branches, ring closures, bond orders, charges
and `.`-separated components.
`embed_smiles(&SmilesMolecule, &EmbedOptions)` produces deterministic 3D
coordinates with ideal bond lengths and a spring relaxation; rings close
within the documented tolerance and no stereochemistry is claimed.
`coordination_geometry(element)` returns a `CoordinationGeometry` and
`apply_coordination` places metal neighbours on the octahedral, tetrahedral,
square-planar or trigonal-bipyramidal shape, returning a warning when the graph
degree differs.

`Sketch` is a 2D editable graph (`SketchAtom`, `SketchBond`, `add_atom`,
`add_bond`, `remove_atom`); `to_smiles` writes a deterministic SMILES and
`from_smiles` builds a sketch from a parsed molecule.
`RDKIT_BEHAVIOUR_INVENTORY` records the reference behaviours (aromaticity,
sanitization, conformers, stereo, fragments) and what this crate implements.

## Fragments, substitution and guests

`fragments()` returns the embedded library (water, methanol, acetone and
others) with name, SMILES, category and description; `find_fragment` looks one
up by name. `clearance(structure, point)` returns the minimum periodic
distance to an atom centre; `substitute(structure, site, element)` changes one
asymmetric site. `place_guest(host, guest, offset, clearance)` inserts a guest
molecule when every atom keeps the requested periodic clearance from the host.

`add_hydrogens(structure, &HydrogenOptions)` fills the missing valences from
the perceived bond graph with deterministic geometry (C 4, N 3, O 2, B/P/S 3
or 2). The CLI exposes it as `xtalforge hydrogens`.

## Twists, ribbons, spirals and ZTC

`commensurate_twist(n, m, theta_degrees, layers)` builds two graphene layers
with a relative twist. Zero and 60 degrees are exact; another angle is
approximated by the nearest of the two and recorded in the structure metadata
(`twist_requested`, `twist_approximation`). `nanoribbon(width, length, vacuum)`
adds vacuum to a periodic honeycomb strip.

`ZtcOptions { samples, cell, isovalue }` controls the level-set surfaces.
`build_ztc` samples the gyroid `sin x cos y + sin y cos z + sin z cos x` and
`build_schwarzite` the Schwarz P surface `cos x + cos y + cos z = 0`;
`ztc_surface_mesh` and `schwarzite_surface_mesh` return the periodic
`TriangleMesh`. The result is a geometric carbon model, not a relaxed
amorphous structure.

## Polymers

`build_polymer(smiles, monomers, &PolymerOptions)` parses a repeat unit,
embeds it, replicates it and packs the chain, recording every step in
`PolymerBuild::steps`. `pack_chains(&[Structure], density)` fills a cubic cell
at a target density with a deterministic grid and a clash check.
`polymer_protocol()` loads the documented equilibration protocol (minimise,
anneal, NVT, NPT, production) from `data/polymer_protocol.json`.

## MOF frameworks

`xtalforge-mof` assembles a framework from a topology, a node block and a
linker block. Topologies, blocks, rods, layers, merged nets, decoration,
invented linkers and the pore-space partition are documented in the
[MOF builder](./mof) page.

## Surfaces

| Operation | Rust | CLI | Python | WASM |
|---|---|---|---|---|
| Graphene | `graphene` | `build graphene` | `build_graphene` | `build_graphene_info` |
| Nanotube | `nanotube` | `build nanotube` | `build_nanotube` | `build_nanotube_info` |
| C60 | `fullerene(60)` | `build c60` | `build_c60` | `build_c60_info` |
| C20 | `face_spiral(20)` | `build c20` | `build_c20` | `build_c20_info` |
| Twist | `commensurate_twist` | `build twist` | `twist` | `twist_info` |
| Ribbon | `nanoribbon` | `build ribbon` | `ribbon` | `ribbon_info` |
| ZTC | `build_ztc` | `build ztc` | `build_ztc` | `ztc_info` |
| Schwarz P | `build_schwarzite` | `build schwarzite` | `build_schwarzite` | `schwarzite_info` |
| Polymer | `build_polymer` | `build polymer` | `build_polymer` | `polymer_info` |
| Hydrogens | `add_hydrogens` | `hydrogens` | `add_hydrogens` | `hydrogens_info` |
