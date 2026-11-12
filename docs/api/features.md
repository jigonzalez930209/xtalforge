# Core features (Phase 9)

Phase 9 adds the operations a standard crystal library offers beyond the
Crystal-Builder parity of Phases 1-8. Every capability below exists in the
Rust API, the CLI, Python and WASM; the [parity matrix](./parity.md) records
the exact status. The CLI groups the commands under `xtalforge features`,
with `compare`, `group` and `deform` as top-level commands.

## Structure comparison

### `match_structures`

Decides whether two structures are the same up to lattice basis, origin,
atom permutation and the `ltol`, `stol` and `angle_tol` tolerances. The
search works on Niggli-reduced cells and accepts a supercell relationship
when `allow_supercell` is set.

```bash
xtalforge compare a.cif b.cif --stol 0.3 --ltol 0.2
```

```python
result = xtalforge.match_structures(a, b, ltol=0.2, stol=0.3)
result["matched"], result["rms"], result["mapping"]
```

The result carries `matched`, `rms`, `mapping`, `translation`,
`lattice_mapping` and `supercell`. The WASM export is
`match_structures(first, second, format_hint, options_json)`.

### `periodic_rmsd` and `kabsch_alignment`

`periodic_rmsd` measures a matched pair under an explicit atom mapping with
the minimum image of the first cell. `kabsch_alignment` returns the optimal
rigid rotation, translation, RMSD and aligned points for two Cartesian point
sets (Horn quaternion form, no external linear algebra).

### `group_structures`

Groups a list of structures into deterministic equality classes.

```bash
xtalforge group a.cif b.cif c.cif
```

## Structure edits

### `deform`

Applies a strain tensor (six engineering Voigt components), a deformation
gradient (nine components), an isotropic volume change or a target density.

```bash
xtalforge deform a.cif --voigt 0.01,0,0,0,0,0
xtalforge deform a.cif --gradient 1.01,0,0,0,1.01,0,0,0,1.01
xtalforge deform a.cif --density 3.2 --out scaled.cif
```

### `rattle`

Draws Gaussian or bounded displacements with an explicit seed. The report
has `max_displacement` and `rms_displacement`; the same seed always gives the
same structure.

### `edit_structure`

Adds, removes, replaces, moves, sorts and wraps sites, applies a species map
(`{"Si": "Ge"}`), and can replace a fraction of the matching sites with a
seed. The same operations are undoable through the `xtalforge-commands`
crate (`AddSites`, `DeleteSites`, `SetElement`, `MoveSites`,
`ReplaceSpecies`, `SortSites`, `WrapIntoCell`).

```bash
xtalforge features edit a.cif --species Si:Ge --fraction 0.5 --seed 7
```

### `interpolate`

Linear interpolation between two matched structures. The CLI writes every
frame to an `.xyz` path and otherwise emits the first frame; Python returns a
list of structure dictionaries.

### `optimal_supercell`

Searches integer supercell matrices for the one closest to a cube at a target
cell count (ASE `find_optimal_cell_shape`).

## Local environment and chemistry

### `voronoi`

Periodic Voronoi cells built exactly from bisector half-spaces, with the
volume, face area and solid angle of every neighbor.

### `coordination`

Coordination numbers with the `cutoff`, `voronoi` or `crystalnn` strategy,
plus continuous shape measures against linear, trigonal planar, tetrahedral,
square planar, trigonal bipyramid, square pyramid, octahedral, pentagonal
bipyramid, cube, square antiprism and icosahedral ideals.

### `bond_valence_sums` and `oxidation_states`

Bond-valence parameters from O'Keefe and Brese (1991) are embedded with their
citation. Oxidation states come from a maximum-a-posteriori search over the
ICSD bond-valence distributions with a charge-balance penalty; the report
lists the state per site, per element, the total charge and the score.

### `radial_distribution` and `structure_fingerprint`

Total and partial RDF/PDF with optional smoothing, running coordination
numbers, bond-angle distributions and a composition + cell + RDF fingerprint
with a cosine similarity.

## Surfaces

### `unique_miller_indices`

Enumerates symmetry-unique `(h k l)` up to a maximum index under the point
group of the crystal, with multiplicity, family and d-spacing.

### `slab_terminations`

Generates every distinct termination of one `(h k l)`, optionally mirrors the
slab into a symmetric one, and reports the thickness, net dipole and polar
flag of each termination. A minimal slab for a requested thickness is
available through `minimal_slab`.

### `wulff`

Builds the Wulff polyhedron from `(h k l, energy)` pairs, with facet areas,
fractions, volume and surface area. The CLI expands symmetry-equivalent
directions with the detected point group.

### `adsorption_sites`

Enumerates top, bridge and hollow sites of the uppermost slab layer.
`place_molecule` puts a molecule above a chosen site with an offset along the
surface normal.

## Defects, disorder and interfaces

### `vacancies`, `substitutions`, `interstitials`

Symmetry-unique point defects. Vacancies and substitutions use the detected
orbits; interstitial candidates are the Voronoi vertices with a clearance
filter. Each defect is expanded to a supercell that reaches the requested
minimum image distance.

### `orderings`

Enumerates decorations of partial occupancies in a supercell, removes
symmetry-equivalent decorations with the structure matcher, and ranks the
remainder by the Ewald energy of the ordered cell.

### `special_quasirandom`

Generates a binary special quasirandom structure with a fixed seed and reports
the pair correlation of every distance shell and the nearest-neighbor triplet
correlation.

### `lattice_matches` and `stack_slabs`

Zur-McGill-style 2D lattice matching with a strain and area report, and
heterostructure stacking with a gap, in-plane shift and vacuum.

## Mechanical and thermodynamic basics

### `elastic_tensor`

Applies the six Voigt strains with a central difference, reads the analytic
stress of any `Calculator`, symmetrizes the 6x6 stiffness, inverts it to the
compliance and reports Voigt, Reuss and Hill moduli, Young's modulus,
Poisson's ratio and a positive-definiteness stability flag.

```bash
xtalforge features elastic si.cif --strain 0.005
```

### `eos_fit`

Fits Birch-Murnaghan (third order), Vinet or Murnaghan energy-volume data with
Levenberg-Marquardt, reporting `v0`, `e0`, `b0` (GPa), `b0'`, the RMS residual
and parameter errors.

```bash
xtalforge features eos --point 18,1 --point 19,0.4 --point 20,0 --point 21,0.5 --point 22,1.7
```

### `convex_hull`

Formation energies per atom and the distance to the lower convex hull for a
list of phases, with a two-phase simplex that supports any number of
elements. Elemental references are derived from the pure-element entries or
supplied explicitly.

```bash
xtalforge features hull phases.json
```

## Reproducibility

Random algorithms (`rattle`, partial `edit_structure`, `special_quasirandom`)
take an explicit seed. Every search sorts its output deterministically, and
the Ewald, Voronoi and shape-measure kernels use fixed iteration orders, so
the same input and options give the same result on every surface.

## Types and lower-level helpers

The Phase 9 modules expose typed results and a few building blocks that the
higher-level functions share:

- Comparison: `Alignment` (rotation, translation, RMSD and aligned points),
  `structure_fingerprint` returns a `StructureFingerprint` and
  `radial_distribution` returns a `RadialDistribution`; the options are
  `RdfOptions`.
- Edits: `RattleOptions` and `RattleReport`, `SpeciesReplacement`, `SiteSort`
  and `SupercellChoice`; the site helpers are `add_sites`, `remove_sites`,
  `replace_site`, `replace_species`, `move_site`, `sort_sites`, `wrap_sites`,
  `scale_to_volume` and `scale_to_density`, with `voigt_to_tensor` and
  `tensor_to_voigt` for the strain boundary.
- Local environments: `VoronoiReport`, `VoronoiCell` and `VoronoiNeighbor`;
  `CoordinationShell`, `CoordinationNeighbor`, `CoordinationStrategy` and
  `ShapeMeasure` with `shape_measures`; `AngleDistribution` with
  `bond_angle_distribution`.
- Surfaces and geometry: `Polyhedron` and `HalfSpace` (the half-space clipper
  shared by Voronoi and Wulff, with `FULL_SPHERE` as the reference solid
  angle), `SurfaceEnergy`, `WulffShape`, `WulffFacet`, `AdsorptionSite` and
  `AdsorptionKind`.
- Thermodynamics: `EosPoint`, `EosKind` and `EosFit` with `eos_energy`;
  `PhaseEntry`, `PhasePoint` and `ConvexHull` from `convex_hull`.
