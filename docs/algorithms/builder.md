# Builders: molecules, polymers, carbon and nets

Reference: Crystal-Builder `xtal/build/` and sci-form `smiles.rs` (MIT). The
builders produce P1 structures with explicit bonds; they do not infer a space
group.

## SMILES and 3D embedding

`parse_smiles` builds the graph (atoms, bonds, aromatic flags, charges) and
`embed_smiles` places it in 3D. The embedding is a distance-geometry pipeline:
bond lengths and angles set the initial distances, a deterministic
pseudo-random seed breaks symmetry, and a refinement minimises the distance
error. The same SMILES and seed give the same coordinates; the seed is an
argument, never a hidden global.

Aromatic rings use the pyykko aromatic radius; the parser records the Kekulé
form so a writer can round-trip it. Unknown elements and malformed ring
closures are typed errors.

## Fragments and substitution

The fragment library is a name, SMILES, category and description. A fragment is
embedded on demand. `clearance` reports the shortest distance from a point to
the existing atoms, and `substitute` replaces a site while keeping the
occupancy, label and charge metadata. A substitution that would overlap an atom
is refused unless the caller passes the documented override.

## Hydrogens

`add_hydrogens` reads the perceived bond graph, compares each C/N/O/B/P/S atom
with its target valence and places the missing H along the free directions.
Directions are computed through the lattice minimum image, so a hydrogen on a
cell boundary is not duplicated. Element defaults are C-H 1.09, N-H 1.01 and
O-H 0.96 A; setting any other distance uses it for every element.

The function assumes the bond rules are right: a wrong perception gives wrong
hydrogens. The diagnostics page lists `missing_hydrogens` as the signal that
this step is needed.

## Polymers

`build_polymer(smiles, monomers, density)` embeds one repeat unit, polymerises
along the connection points and packs chains until the target density is
reached. Packing is deterministic for the same arguments. The result is P1 with
explicit bonds; the density, monomer count and SMILES are recorded in metadata.

## Carbon allotropes

- `graphene(nx, ny, layers, interlayer, vacuum)` is an AA-stacked honeycomb
  sheet. The cell is `nx a1` and `ny a2` with `a1 = (a, 0)`,
  `a2 = (a/2, a sqrt(3)/2)` and `a = 1.42 sqrt(3)`; the two atoms of primitive
  cell `(i, j)` sit at fractional `((i + u)/nx, (j + v)/ny)` with
  `(u, v) = (0, 0)` and `(1/3, 1/3)`. Every bond is 1.42 Angstrom.
- `nanotube(n, m, length, vacuum)` rolls the honeycomb along the chiral vector
  `C = n a1 + m a2`. The strip is `|C|` wide and `length` translation vectors
  `T = ((2m + n) a1 - (2n + m) a2) / d_R` long, with `d_R = gcd(2n + m, 2m + n)`,
  and holds `4 (n^2 + n m + m^2) / d_R` atoms per `T`. Each atom goes to
  `(R sin t, y, R cos t)` with `R = |C| / 2 pi` and `t = 2 pi (r . C) / |C|^2`,
  in a cell `2R + vacuum` wide on x and z with the axis along y; the bonds are
  chords, about 1 % shorter than 1.42 Angstrom for a narrow tube.
- `fullerene(60)` uses the built-in C60 coordinates; `face_spiral(20)` gives
  C20 by a documented face traversal.
- `nanoribbon(width, length, vacuum)` is a honeycomb strip periodic along `a`;
  the `b` vector grows by `vacuum` while the atoms keep their Cartesian
  positions.
- `build_ztc`/`build_schwarzite` sample a periodic trigonometric level set and
  extract the isosurface with the analysis crate.
- `build_carbon` is the disordered recipe documented in [carbon](./carbon).

## Nets and guests

`build_net` turns a CGD entry into a P1 structure whose sites are vertices and
whose topology bonds are edges. `place_guest` puts a molecule into a pore with
a clearance check. Interpenetration is a separate command over a document.

## Validation

- `crates/xtalforge-build/tests/` covers SMILES round-trips, embedding
  determinism, fragment lookup, hydrogen counts, polymer density, C20/C60
  vertex and face counts, graphene/nanotube ring census, twist/ribbon
  periodicity and the carbon golden data.
- The C60 coordinates are a recorded data file, not a generated geometry; the
  test checks 60 vertices and 90 edges.
- The carbon recipe is validated against the Crystal-Builder distributions in
  `tests/golden/phase-6/carbon/`, with the documented remesh/ribbon
  deviations.

## Differences from Crystal-Builder

- The embedder is deterministic from a seed; the reference may use a random
  initialisation.
- A substitution overlap is refused rather than warned.
- The polymer packer reports the achieved density; it does not silently
  continue when the target is out of range.
