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

`commensurate_twist(n, m, theta_degrees, layers)` builds a twisted graphene
multilayer in the smallest commensurate moire cell whose angle is within
0.05 degrees of `theta_degrees` (or the closest one), repeated `n x m` times in
plane. The metadata records `twist_requested`, the exact `twist_angle`, the
moire indices `twist_indices` and `twist_approximation` (equal to
`twist_angle`). `(1, 1, 21.79, 2)` gives the `(2, 1)` cell: 21.787 degrees and
28 atoms. `nanoribbon(width, length, vacuum)`
adds vacuum to a periodic honeycomb strip.

`ZtcOptions { samples, cell, isovalue }` controls the level-set surfaces.
`build_carbon(&CarbonRecipe)` builds a disordered carbon on a periodic net:
a convolution surface of Gaussians along the net edges, an isotropic
remesh to the 2.46 Angstrom target, a ribbon cut that follows the net, a
three-coordinate dual and the requested edge termination. It returns
`CarbonBuild { structure, recipe, census, density, ratios, edge_fraction,
pieces, periodicity, edge_length, closed_carbons, bare_edges, mean_bond,
notes }`. `CarbonSurface` exposes the remeshed sheet for tests. The CLI is
`build carbon`, Python `build_carbon`, WASM `build_carbon_recipe_info` (statistics) and `build_carbon_text` (CIF for a viewer).
See [Disordered carbons](../algorithms/carbon.md).

`build_ztc` samples the gyroid `sin x cos y + sin y cos z + sin z cos x` and
`build_schwarzite` the Schwarz P surface `cos x + cos y + cos z = 0`;
`ztc_surface_mesh` and `schwarzite_surface_mesh` return the periodic
`TriangleMesh`. The carbon structures are the relaxed dual of the remeshed
surface: one three-coordinate carbon per triangle, 1.42 Angstrom bonds, with
explicit bonds stored. The carbon builders sample their own ~1.2 Angstrom
grid; `samples` sets the resolution of the returned meshes. `build_ztc`
refuses cells below `MIN_GYROID_CELL` (12 Angstrom) and `build_schwarzite`
below `MIN_SCHWARZ_P_CELL` (8 Angstrom), where fewer than 90 % of the carbons
would be three-coordinate. `build_carbon` relaxes its dual
carbon graph to 1.42 Angstrom bonds and 120 degree angles before terminating
the edges, and drops a termination that would overlap another atom.

## Polymers

`build_polymer(smiles, monomers, &PolymerOptions)` parses a repeat unit,
embeds it, joins `monomers` copies head (first SMILES atom) to tail (last
atom) into a chain that is periodic along `a`, and sizes the square
cross-section for the target density, recording every step in
`PolymerBuild::steps`. `pack_chains(&[Structure], density)` fills a cubic cell
at a target density with a deterministic grid and a clash check; it does not
unwrap a chain longer than the cell.
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
| Disordered carbon | `build_carbon` | `build carbon` | `build_carbon` | `build_carbon_recipe_info` |
| Schwarz P | `build_schwarzite` | `build schwarzite` | `build_schwarzite` | `schwarzite_info` |
| Polymer | `build_polymer` | `build polymer` | `build_polymer` | `polymer_info` |
| Hydrogens | `add_hydrogens` | `hydrogens` | `add_hydrogens` | `hydrogens_info` |
