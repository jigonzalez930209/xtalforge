# Phase 6 — Builders

**Goal:** everything Crystal-Builder can build from nothing or from parts:
molecules (SMILES, library, coordination geometry), substitution and guest
filling, disordered carbons (ZTC, schwarzites), graphene and stacks with
commensurate twists, nanotubes, fullerenes, amorphous polymers, nets drawn as
structures, and the MOF builder (a Rust port of the vendored PORMAKE plus
Crystal-Builder's extensions: rods, merged nets, layers, supertetrahedral
decoration, pore-space partition, inverse linker design).

**Reference:** `Crystal-Builder/xtal/build/` (molecule, chem 1352 (RDKit),
sketch, library + data/fragments.json, substitute 464, coordination,
clearance, fill 718, topology), `xtal/carbon/` (build, lattice, honeycomb,
mesh 719, surface 434, ribbons), `xtal/graphene/` (build, twist),
`xtal/nanotube/`, `xtal/fullerene/` (build, spiral, geometry, catalogue +
data), `xtal/polymer/` (monomer, chain, pack 755, build, protocol),
`xtal/mof/` (19 modules, about 15 000 lines) and `xtal/mof/pormake/`
(vendored PORMAKE 0.2.3, MIT, see its `PROVENANCE.md`), `xtal/mof/data/`,
`xtal/mof/library/`, `xtal/modules/{build,carbon,graphene,nanotube,fullerene,polymer,mof,net}.py`,
`xtal/commands/connections.py`; sci-form `src/smiles.rs`, `graph.rs`,
`stereo.rs`, `src/materials/{assembly,sbu,topology_data}.rs`.

Crates: `xtalforge-build`, `xtalforge-mof`.

This is the largest phase. Task groups ship independently: 6.1–6.4
(molecules), 6.5–6.8 (carbon family), 6.9 (polymers), 6.10 (nets), 6.11–6.16
(MOF).

---

## Task 6.1 — Molecules from SMILES

- [ ] 6.1.1 SMILES parser and molecular graph
  Port sci-form `smiles.rs` (875 lines), `graph.rs` (627) and the needed parts
  of `stereo.rs` (1121), each split under 500 lines: atoms, bonds,
  aromaticity, charges, isotopes, explicit and implicit hydrogens, ring
  closures, stereo tags, and Crystal-Builder's attachment-point notation
  (`[*:1]`, `[*:2]`) used by the fragment library.
- [ ] 6.1.2 3D embedding
  Crystal-Builder embeds with RDKit (ETKDG) in `build/chem.py`. Use sci-form's
  ETKDG through the optional `sci-form` feature; without the feature, provide
  a reduced distance-geometry embedder good enough for small fragments, or
  return a clear error. Record the decision and the WASM impact.
- [ ] 6.1.3 Coordination geometry for metals
  Port `build/coordination.py`: ETKDG knows no coordination chemistry
  (cisplatin came out with Cl-Pt-Cl angles of 97-120 degrees), so the metal and
  its bonded atoms are placed on an ideal polyhedron (linear, trigonal,
  tetrahedral, square planar, trigonal bipyramidal, square pyramidal,
  octahedral, …) and the rest of the molecule follows.
- [ ] 6.1.4 Built molecule type
  Port `build/molecule.py`: a fragment plus which of its atoms are connection
  points; converts to a `Structure` (molecule in a box) or a clipboard
  `Fragment`.
- [ ] 6.1.5 RDKit behaviour inventory
  Read `build/chem.py` and list everything it asks RDKit for (sanitization,
  hydrogens, embedding, MMFF/UFF cleanup, canonical SMILES, substructure
  matching); map each to a sci-form function or a new subtask, and record
  intentional differences.
- [ ] 6.1.6 Sketch data model
  Port `build/sketch.py` (382 lines) as data only: a molecule as drawn (atoms
  at 2D points, bonds, nothing perceived yet) to graph conversion. The editor
  itself is a GUI and out of scope.

## Task 6.2 — Fragment library

- [ ] 6.2.1 Library data
  Port `build/library.py` and `build/data/fragments.json`: named molecules
  worth not typing twice (DMF `CN(C)C=O`, the phenylene strut
  `[*:1]c1ccc([*:2])cc1`, solvents, linkers), embedded as a data file.

## Task 6.3 — Substitution and clearance

- [ ] 6.3.1 Clearance
  Port `build/clearance.py`: how much room a group of atoms has in a periodic
  cell at a given rotation angle; used to turn a linker about its axis
  (Task 6.13.1) and a substituent about its bond.
- [ ] 6.3.2 Substitute
  Port `build/substitute.py` (464 lines): replace a hydrogen (or any terminal
  atom on one bond, for example the C-F edges of a ZTC) by a group (amino,
  methoxy, nitro, halogens, library fragments), aligned on the bond and turned
  to the angle with most clearance. Exposed as the `SubstituteHydrogens`
  command (Task 3.10.1).

## Task 6.4 — Guest filling

- [ ] 6.4.1 Guest extraction
  Port `guest_molecules` from `build/fill.py` (718 lines; split): lift every
  discrete molecule out of a structure (a solvent drawn in its own box becomes
  a placeable guest).
- [ ] 6.4.2 Placement
  Place copies in the empty space of a host: random positions and
  orientations in the accessible region (distance grid of Task 5.2), clash
  rejection with periodic images, target count or loading, seeded RNG for
  reproducibility. Exposed as `InsertMolecules` (Task 3.10.6).

## Task 6.5 — Honeycomb sheets and graphene

- [x] 6.5.1 Shared honeycomb recipe
  Port `carbon/honeycomb.py` and `carbon/lattice.py`: element validation,
  bond length typed or taken from the elements (Pyykko aromatic radii,
  Task 1.1.3), X-H cap length, and the record turned into a structure; shared
  by the graphene, nanotube, fullerene and carbon builders.
- [x] 6.5.2 Graphene and stacks
  Port `graphene/build.py` (390 lines) and `modules/graphene.py`: sheets of
  carbon, BN or any two elements (two atoms in a hexagonal cell, A at
  (1/3, 2/3), B at (2/3, 1/3), bonds stated, never perceived), one layer or
  stacked AA, AB, ABC, AA', as periodic graphite or a slab in vacuum, in the
  hexagonal or the orthogonal cell.
- [ ] 6.5.3 Commensurate twists
  Port `graphene/twist.py`: twisted bilayers exist only at commensurate angles
  indexed by (m, r) (Lopes dos Santos, Peres and Castro Neto 2007); list the
  angles with their cell sizes and atom counts, and build the supercell with
  both layers' integer cells.

## Task 6.6 — Nanotubes

- [x] 6.6.1 Roll-up geometry
  Port `nanotube/geometry.py`: chiral vector (n, m), diameter, chiral angle,
  translational period, atoms per period, for C, BN or any two elements.
- [x] 6.6.2 Builder
  Port `nanotube/build.py` (472 lines) and `modules/nanotube.py`: periodic
  along c (length rounded to whole periods) or a finite tube cut to length
  with its open ends capped (H) or left; bonds stated; a measured report
  (diameter, length, atom count).

## Task 6.7 — Fullerenes

- [ ] 6.7.1 Face spirals
  Port `fullerene/spiral.py`: a cage as its Fowler-Manolopoulos face spiral
  (the twelve positions of the pentagons; C60-Ih is `1 7 9 11 13 15 18 20 22 24
  26 32`), wound up into the dual graph and then the cage graph.
- [x] 6.7.2 Catalogue and geometry
  Port `fullerene/catalogue.py` and its data (isomer tables from C20 to C100),
  and `fullerene/geometry.py` (3D embedding of the cage graph and
  relaxation).
- [x] 6.7.3 Builder
  Port `fullerene/build.py` and `modules/fullerene.py`: a cage from the table
  (size and isomer number) or from a typed spiral, of one element, alone in a
  box, bonds stated.

## Task 6.8 — Disordered carbons (ZTC, schwarzites)

- [ ] 6.8.1 Convolution surface
  Port `carbon/surface.py` (434 lines): a smooth periodic sheet around a net,
  defined as a level set of a sum of Gaussians over beads placed along the
  net's edges (not a distance field).
- [ ] 6.8.2 Mesh
  Port `carbon/mesh.py` (719 lines; split): the sheet meshed periodically
  (reusing the marching tetrahedra of Task 5.4) and remeshed into a
  honeycomb-like carbon network at the target density.
- [ ] 6.8.3 Ribbons
  Port `carbon/ribbons.py`: cut the closed sheet into ribbons that follow the
  net's edges and still make one connected piece (in the example ZTC four
  carbons in ten are on an edge; each dia edge of 10.3 Angstrom carries about
  eighteen carbons).
- [ ] 6.8.4 Recipe builder
  Port `carbon/build.py` and `modules/carbon.py`: a recipe (net, density,
  fraction of sheet kept, layers, defects, edge termination, seed) gives the
  same carbon every time. Steps as in the source: distorted net, sheet at
  scale, carbon placement, ribbons, termination, report.
- [ ] 6.8.5 Validation
  Atom counts, ring statistics (Task 1.12) and density against
  Crystal-Builder output for fixed seeds (the RNG must be identical or the
  comparison is statistical; decide and document).

## Task 6.9 — Amorphous polymers

- [ ] 6.9.1 Monomers and chains
  Port `polymer/monomer.py` and `polymer/chain.py`: a repeat unit with head
  and tail atoms, torsion choices, tacticity; user monomers read from the
  workspace `monomers/` folder (Task 3.12.1).
- [ ] 6.9.2 Packing
  Port `polymer/pack.py` (755 lines; split): every chain grows at once from
  random seeds, one unit per chain in shuffled rounds, so no chain wanders
  through an empty box while the last ones are squeezed; clash checks with
  periodic images; start density.
- [ ] 6.9.3 Build steps
  Port `polymer/build.py` and `modules/polymer.py`: grow, push apart or
  compress to the target density, state the bonds, report the numbers
  (density, chain lengths, end-to-end distances).
- [ ] 6.9.4 Equilibration protocol as data
  Port `polymer/protocol.py`: how a packed model would be equilibrated, held as
  data, plus the extension point where an equilibrator (an MD engine, not part
  of this project) plugs in. Nothing runs.

## Task 6.10 — Nets drawn as structures

- [ ] 6.10.1 Net builder
  Port `build/topology.py` and `modules/net.py`: an RCSR net as a structure
  (vertices as dummy or chosen atoms, edges as topology bonds) so it can be
  viewed, measured and used as a template. No optional dependency.

## Task 6.11 — MOF data and licensing

- [ ] 6.11.1 License record
  PORMAKE 0.2.3 is MIT ("Copyright (c) 2022 Sangwon"), vendored in
  Crystal-Builder with every modification listed in its `PROVENANCE.md`. Copy
  the notice into `THIRD_PARTY.md`. Verify separately the license of the
  PORMAKE database (2403 topologies, about 867 building blocks; the database
  may come from RCSR and CoRE MOF sources), of `mof/data/` and of
  `mof/library/`.
- [ ] 6.11.2 Compact bundle
  A reproducible script converts the topology `.cgd` files and the building
  block `.xyz` files into one versioned compressed bundle (`data-YYYY.MM.DD`),
  with an index for lazy access. Crystal-Builder measured 2.85 MB for nets and
  blocks; the bundle should not be larger.
- [ ] 6.11.3 Distribution strategy
  Feature `mof-data` embeds the bundle; without it, the bundle is fetched on
  first use and cached with a checksum (the default for WASM, where size
  matters, and an option for Python wheels).

## Task 6.12 — MOF builder: PORMAKE port

- [ ] 6.12.1 Topologies and building blocks
  Port `pormake/topology.py`, `building_block.py`, `local_structure.py`,
  `neighbor_list.py`, `utils.py`: nets with node and edge slots, building
  blocks with connection points (dummy `X` atoms, Task 1.1.1). Replace the
  vendored ASE usage with our own types.
- [ ] 6.12.2 Locating and scaling
  Port `pormake/locator.py` (fit a block onto a slot by RMSD of connection
  directions) and `pormake/scaler.py` (scale the net so every edge length
  matches the blocks), replacing the removed gradient and function call
  described in `PROVENANCE.md` with Rust equivalents.
- [ ] 6.12.3 Builder and framework
  Port `pormake/builder.py` and `framework.py`: place nodes and edges, join
  connection points, remove dummies, produce the framework and its RMSD
  report.
- [ ] 6.12.4 Edge midpoints
  Port `mof/midpoints.py`: PORMAKE's reading of a net (nodes plus edge
  centres) with Crystal-Builder's third method of finding edges.
- [ ] 6.12.5 Compare with sci-form assembly
  sci-form `materials/assembly.rs` (829 lines) also assembles frameworks;
  compare on simple nets (pcu, dia, sql) and reuse what is better tested.

## Task 6.13 — MOF builder: orientation, catalogue, build request

- [ ] 6.13.1 Orientation
  Port `mof/orient.py` (1354 lines; split): a node block is placed by fitting
  its connection directions onto the slot's, and a symmetric node has many
  equally good fits (24 for an octahedral node, RMSD spread 6e-08) that are
  not cosmetic for the rest of the block; choose deterministically and
  explore when asked. `align_edges` turns linkers about their axis using
  clearance (Task 6.3.1).
- [ ] 6.13.2 Catalogue
  Port `mof/catalog.py` (1306 lines; split): what there is to build with, read
  from the bundle files rather than through PORMAKE (import time was 10 to 30
  s in Python); search by topology, connectivity, elements, size; user blocks
  from the workspace `blocks/` folder; the stamp/refresh logic for the
  catalogue cache.
- [ ] 6.13.3 Build request and result
  Port `mof/build.py` (1784 lines; split into resolve, run, read back, net
  overlay): a request is three strings (topology, node assignment, edge
  assignment) so it survives the CLI and saved parameters; the framework is
  read back with its net drawn on it as topology bonds.
- [ ] 6.13.4 Block writing and connection points
  Port `mof/block.py` (connection point placement at `CONNECTION_DISTANCE`,
  block file format) and `commands/connections.py` (`MarkConnectionPoints`,
  `MarkOneConnectionPoint`; one way only, as documented in the source).
- [ ] 6.13.5 Attachment semantics
  Port `mof/attach.py`: what one connection point stands for (one atom, or a
  set of atoms as in MFU-4l's Zn5Cl4 kernel) and which bonds it carries into
  the framework.

## Task 6.14 — MOF builder: rods and layers

- [ ] 6.14.1 Rod SBUs
  Port `mof/rod.py` (1075 lines; split): a rod is one 1-periodic block (MOF-74,
  MIL-53, MOF-303, Ni2Cl2BTDD, Co(BDP)), not a chain of finite nodes (cut into
  nodes, MOF-74's net msf would put carboxylate carbons, not zinc, on the
  vertices).
- [ ] 6.14.2 Rod nets
  Port `mof/rodnet.py` (1178 lines; split): parallel rod axes in one of the
  four plane patterns sql, hxl, hcb, kgm, rod handedness and placement.
- [ ] 6.14.3 Rod builder and donors
  Port `mof/rodbuild.py` (1056 lines; split) and `mof/donors.py`: linkers
  meeting a rod as a core with `X` where each donor joined, recovering which
  `X` was which donor from a whole drawn linker (dobpdc versus pc-dobpdc).
  Validate with the rod test files in `Crystal-Builder/tests/data/rod_*.xyz`.
- [ ] 6.14.4 Layered frameworks
  Port `mof/layers.py`: 2-periodic nets (Ni3(HITP)2, Cu-HHTP, 2D COFs)
  shipped as `.cgd` in a 3D cell and stacked along the free direction.

## Task 6.15 — MOF builder: merged nets, decoration, partition, trimers

- [ ] 6.15.1 Merged nets
  Port `mof/merged.py` (839 lines) and `mof/merged_nets.py` with
  `mof/data/merged_nets.{cgd.gz,json}`: nets with several kinds of edges
  (`sph` is spn plus hxg) need a set of linkers whose reaches put every edge at
  its drawn length at one scale.
- [ ] 6.15.2 Supertetrahedral decoration
  Port `mof/decorate.py` (1708 lines; split): any 4-connected net as an ST-MOF
  net (MIL-100, sod-ZMOFs): four M3O trimers at the corners of a tetrahedron,
  a tritopic linker on each face, every trimer shared by two supertetrahedra.
- [ ] 6.15.3 Trimers
  Port `mof/trimer.py`: an M3O(RCO2)6 trimer with its open metal sites as
  connection points.
- [ ] 6.15.4 Pore-space partition
  Port `mof/partition.py` (624 lines): which partitioner fits which `acs`
  framework (MIL-88, MOF-235; Bu and Feng's pacs frameworks put a C3 pyridyl
  ligand on the channel axis, its nitrogens on the trimers' open sites).

## Task 6.16 — MOF builder: inverse linker design

- [ ] 6.16.1 Linkers by span
  Port `mof/linkers.py` (447 lines) with `mof/data/linker_series.json` and
  `linker_set_nets.json`: a `Requirement` (distance between connection
  points, topicity, geometry) and `suggest` returning molecules that satisfy
  it; used by partition, merged nets and ST-MOF windows.
- [ ] 6.16.2 Invented linkers
  Port `mof/invent.py`: put together new linkers wider than the stored
  isoreticular series on request.
- [ ] 6.16.3 Validation
  Rebuild known MOFs from the corpus (HKUST-1, MOF-5, UiO-66, ZIF-8,
  MIL-100, MOF-74, Ni3(HITP)2) and compare cell, atom count, net, and RMSD
  with Crystal-Builder output and with the corpus CIFs.

## Task 6.17 — Surfaces and docs

- [ ] 6.17.1 CLI
  `build molecule` (SMILES or library name), `build substitute`, `build fill`,
  `build graphene`, `build nanotube`, `build fullerene`, `build carbon`,
  `build polymer`, `build net`, `mof {list,search,build,rod,layer,decorate,partition,linkers}`.
- [ ] 6.17.2 Python and WASM
  Builders return `Structure` plus a report; long builds are cancellable and
  report progress; WASM loads the MOF bundle on demand (Task 6.11.3).
- [ ] 6.17.3 Docs and demos
  Algorithm pages: coordination polyhedra, clearance search, commensurate
  twist angles, nanotube roll-up, face spirals, convolution surfaces, chain
  growth packing, PORMAKE locate and scale, node orientation ambiguity, rods,
  merged nets, supertetrahedral decoration, pore-space partition. Demos:
  nanotube and graphene builders, building a MOF in the browser.
- [ ] 6.17.4 Parity and release
  Parity matrix; changelog; MINOR release at phase completion (task groups
  may ship earlier as PATCH releases).
