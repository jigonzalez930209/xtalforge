# Builders: molecules, polymers, carbon and nets

Reference: Crystal-Builder `xtal/build/` and sci-form `smiles.rs` (MIT). The
builders produce P1 structures with explicit bonds; they do not infer a space
group.

## SMILES and 3D embedding

`parse_smiles` builds the graph (atoms, bonds, aromatic flags, charges) and
`embed_smiles` places it in 3D. The embedder is a small distance-geometry
relaxation, not ETKDG: the atoms start on a seeded Fibonacci spiral and a
gradient descent drives three terms, bond lengths towards element-pair
values, 1-3 distances towards the ideal bond angle of the shared atom
(180 degrees for a triple bond or cumulated double bonds, 120 degrees for an
aromatic or double-bonded atom, 109.47 degrees otherwise; three- and
four-membered rings keep only their bonds) and a non-bonded repulsion. The
same SMILES and seed give the same coordinates; the seed is an argument,
never a hidden global. Benzene, naphthalene, cyclohexane, neopentane and
propyne are checked for their bond lengths and angles in
`crates/xtalforge-build/tests/embed.rs`. Before the 2026-11-11 audit the angle
term was missing and rings collapsed (benzene had 1.17 Angstrom 1-3
contacts).

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

`build_polymer(smiles, monomers, density)` embeds one repeat unit and rotates
it so that the first SMILES atom (the head) and the last one (the tail) lie on
the `a` axis. The units are joined head to tail, each tail bonded to the next
head at the sum of their covalent radii, and the last unit bonds to the first
through the periodic boundary, so the chain is infinite along `a`. The cell is
orthorhombic: `a` is the chain period and the square cross-section gives the
target density from the real atomic masses. If that cross-section would bring
the chain within 2.5 Angstrom of its periodic image the cell is widened and
the achieved density is reported in the steps. Write the backbone first and
last (`C(c1ccccc1)C` for polystyrene); hydrogens are not added, run
`add_hydrogens` afterwards. Before the 2026-11-11 audit the chain was packed
into a cube smaller than its own length, which overlapped carbons at
0.67 Angstrom.

## Carbon allotropes

![Graphene cell and nanotube roll](/diagrams/carbon-geometry.svg)


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
- `commensurate_twist(n, m, theta, layers)` builds a twisted multilayer in the
  commensurate moire cell `L1 = i a1 + j a2`, `L2 = R(60) L1`, whose twist
  obeys `cos theta = (i^2 + 4ij + j^2) / (2 (i^2 + ij + j^2))` and which holds
  `2 (i^2 + ij + j^2)` atoms per layer. The smallest cell within 0.05 degrees
  of the request is chosen (indices up to 20), otherwise the closest angle;
  `n x m` is the in-plane repeat of the moire cell. Even layers keep the
  reference orientation, odd layers are rotated, 3.35 Angstrom apart.
  `(2, 1)` gives 21.787 degrees and 14 atoms per layer. Before the
  2026-11-11 audit the angle was rounded to 0 or 60 degrees.
- `nanoribbon(width, length, vacuum)` is a honeycomb strip periodic along `a`;
  the `b` vector grows by `vacuum` while the atoms keep their Cartesian
  positions.
- `build_ztc`/`build_schwarzite` sample the gyroid or Schwarz P level set on
  a periodic grid of about 1.2 Angstrom (including the closing plane, and
  nudging samples that sit exactly on the isovalue), weld the marching
  tetrahedra mesh, remesh it towards 2.459 Angstrom edges (the dual of a
  1.42 Angstrom honeycomb), place one carbon per triangle and relax the dual
  graph to 1.42 Angstrom bonds and 120 degree angles. For cells of 12
  Angstrom and above at least 97 % of the carbons are three-coordinate, with
  5-, 6-, 7- and 8-membered rings as in real schwarzites. Before the
  2026-11-11 audit the carbons were the mesh vertices thinned to 1 Angstrom
  (coordination 5-8) and the sampled grid lacked the closing plane, which left
  an open seam on the cell faces.
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
- The polymer builder reports the achieved density; it does not silently
  continue when the target is out of range.
