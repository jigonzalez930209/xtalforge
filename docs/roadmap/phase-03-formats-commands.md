# Phase 3 — Formats, commands, workspace and view data

**Goal:** every file format Crystal-Builder reads or writes, project files that
never leave a half-written file, trajectories, the undoable command layer
through which all mutations go, selection and atom tracking, the on-disk
workspace of runs, typed parameters, and the view data (draw styles,
colormaps, scene data, STL) that any viewer can consume.

**Reference:** `Crystal-Builder/xtal/io/*` (registry in `io/__init__.py` and
`io/registry.py`), `xtal/commands/*`, `xtal/core/{selection,tracking,atom_groups}.py`,
`xtal/workspace.py` (1096 lines), `xtal/params.py`, `xtal/plugins.py`,
`xtal/view/*`, the STL part of `xtal/modules/blender.py` and
`xtal/modules/data/printable_stl.py`.

Formats registered by Crystal-Builder (name: extensions, read/write):
`cif` (.cif, .mcif; r/w), `xtalproj` (r/w), `shelx` (.res, .ins; read),
`poscar` (.poscar, .vasp; r/w), `pmg-json` (.json; r/w), `cssr` (r/w),
`gen` (r/w), `xyz` (.xyz, .extxyz; r/w), `lammps-data` (.data, .lmp; r/w),
`pdb` (write only). Non-structure formats: `cube` (grids), `xy` (patterns),
`cgd` (nets), multi-frame extended XYZ (trajectories).

Crates: `xtalforge-io`, `xtalforge-commands`, `xtalforge-view`.

---

## Task 3.1 — Format registry

- [x] 3.1.1 `Format` trait and registry
  Port `io/registry.py`: a format has a name, extensions, optional reader and
  writer. The registry is the only list: CLI help, file dialogs of any front
  end, drag and drop and `convert` all read it. Built-in formats register at
  startup; downstream crates register theirs through a public API
  (replacement for Python entry points in `plugins.py`).
- [x] 3.1.2 Dispatch and content sniffing
  Choose the format from the extension, an explicit name, or, when the
  extension is missing or misleading (`.json`, `.data`, `.txt`), from the
  content (CIF `data_` block, POSCAR layout, XYZ atom count line, pymatgen
  `@class`).
  Note: `registry::sniff` covers CIF, pymatgen JSON and XYZ; extension/name
  lookup covers every registered format.
- [x] 3.1.3 Format capabilities and loss warnings
  Each format declares what it preserves (cell, space group, occupancy,
  ADPs, charges, bonds, labels, frames). A conversion that loses data
  returns warnings (`convert` prints them; `--strict` makes them errors).
  Note: `Capabilities` and `losses(from, to)`; CLI `--strict` follows.
- [x] 3.1.4 Export policy
  Port `io/export.py`: what a file leaving the library should contain. A
  document-preserving export (CIF into a workspace, carrying everything)
  differs from a clean export for another program (dummy atoms removed,
  hydrogens kept, P1 or asymmetric unit by format, explicit bonds as drawn).
  Options are typed and documented per format.

## Task 3.2 — Atomic writes

- [x] 3.2.1 Never leave a half file
  Port `io/atomic.py`: write to a temporary file in the target directory,
  flush and sync, then rename over the target; on error the old file is
  intact and the temporary file removed. Background: Crystal-Builder once left
  a 22-byte empty zip when a project part raised. Used by every writer that
  targets a path.
  Note: `atomic::write_text`/`write_bytes` with a same-directory temp file,
  sync and rename; the test leaves no `.tmp` file.
- [x] 3.2.2 WASM behaviour
  WASM has no filesystem: writers return bytes or strings; this is the
  recorded parity exception.
  Note: the wasm32 branch returns a typed error and callers use the string.

## Task 3.3 — Remaining structure formats

- [x] 3.3.1 pymatgen JSON
  Port `io/pmg.py`: read and write pymatgen's `Structure` dictionary
  (`@module`, `@class`, lattice matrix, sites with species, occupancies and
  fractional coordinates) without importing pymatgen. Symmetry is not
  representable: write P1, read P1. Partial occupancy survives both ways.
  Note: `read_pmg_json`/`write_pmg_json`; registered as `pmg-json`.
- [x] 3.3.2 CSSR
  Port `io/cssr.py`: cell parameters then one line per atom; it is what Zeo++
  reads (Task 5.5.2).
  Note: `read_cssr`/`write_cssr`, P1 writer with the eight zero link slots.
- [x] 3.3.3 DFTB+ `.gen`
  Port `io/gen.py`: `C` cluster (Cartesian, no cell), `S` supercell
  (Cartesian with cell), `F` fractional.
  Note: `read_gen`/`write_gen`, registered as `gen`.
- [x] 3.3.4 LAMMPS data
  Port `io/lammps.py` (336 lines): `atom_style full` (box with tilt factors,
  masses, atoms with molecule id and charge, bonds). The writer takes the
  bonds as shown (the document's bonds, not a re-perceived set). Reader maps
  types back to elements from masses, with warnings when ambiguous.
- [x] 3.3.5 SHELX `.res`/`.ins`
  Port `io/shelx.py` (388 lines): CELL, LATT, SYMM, SFAC, atom lines with
  free variables, and the PART of each atom (disorder). Also extract the
  instruction file embedded in a SHELXL CIF (`_shelx_res_file`,
  `_iucr_refine_instructions_details`) to recover PART information.
  Note: `read_shelx`/`write_shelx`, registered as `shelx`; free variables are
  resolved, SYMM/LATT are metadata and the atoms stay P1.
- [x] 3.3.6 PDB (write)
  Port `io/pdb.py`: a finite set of atoms and bonds as HETATM and CONECT
  records, written for Blender's Atomic Blender importer (element in the
  atom-name columns). Reading PDB is optional and only if a use appears.

## Task 3.4 — Nets: CGD

- [x] 3.4.1 CGD reader and writer
  Port `io/cgd.py` (392 lines): a periodic net as space group, vertices
  (with coordination) and edges. It is Systre's input format and the format
  of the RCSR net database, used by topology analysis (Task 5.6) and the net
  builder (Task 6.10).
  Note: `CgdEntry`, `CgdNode`, `CgdFile`, `read_cgd`/`write_cgd`; malformed
  blocks are collected in `problems`, not fatal.

## Task 3.5 — Volumetric data and patterns

- [x] 3.5.1 `Grid3` type
  A scalar field on a regular grid in a cell (periodic or not), shared by
  cube files, porosity distance grids and isosurfaces.
  Note: `xtalforge_io::grid::Grid3` with index/get/sample/position/range.
- [x] 3.5.2 Gaussian cube
  Port `io/cube.py`: header in Bohr (a negative voxel count means Angstrom
  and is honoured), atoms, volumetric data in the usual 6-per-line layout;
  streaming read for large files.
  Note: `read_cube`/`write_cube`; the whole data block is loaded, with the
  negative-count Angstrom rule.
- [x] 3.5.3 XY patterns
  Port `io/xy.py`: whitespace-delimited 2-theta and intensity, `#` comments,
  optional third column of uncertainties; used by Phase 5 (calculated PXRD
  output) and Phase 7 (measured patterns).
  Note: `read_xy`/`write_xy` with `PatternPoint`.

## Task 3.6 — Trajectories

- [x] 3.6.1 Frame model
  Port `io/trajectory.py` (555 lines): a trajectory is many frames of one
  structure (fixed atom order and elements) with per-frame positions, cell
  (fixed or variable), energy, step number and arbitrary info parsed from the
  extended XYZ comment line (`parse_info`/`format_info`).
  Note: `Trajectory`, `Frame` and the info parser.
- [x] 3.6.2 Multi-frame extended XYZ
  Read and write; lazy indexed reading so a 10 000-frame file is not loaded at
  once; conversion of any frame to a `Structure` (optionally re-attaching the
  space group).
  Note: `XyzIndex` records frame offsets, elements, lattice and header
  metadata and parses one frame on demand; `structure(frame)` remains the
  conversion path.
- [x] 3.6.3 Foreign trajectories
  Crystal-Builder reads other trajectory formats through ASE. Decide which
  formats are needed (VASP XDATCAR, LAMMPS dump, DFTB+ geo_end) and implement
  them natively; record the decision.
  Note: VASP XDATCAR is implemented (`read_xdatcar`, `Direct` and
  `Cartesian`); `read_trajectory` detects XYZ or XDATCAR. LAMMPS dump and
  DFTB+ `geo_end` remain deferred until a consumer asks for them.

## Task 3.7 — Project files (`.xtalproj`)

- [x] 3.7.1 Format compatibility
  Port `io/project.py`: a project is a **zip of small readable JSON parts**
  (`FORMAT_VERSION = 1`): the structure as CIF, `cell.json` (bare cell when
  there are no atoms), `bonds.json` (perceived graph plus hand-drawn and
  suppressed bonds), `sites.json` (per-site `props` a CIF cannot hold),
  `view.json` (view settings), `session.json` (selection, measurements,
  provenance). Read and write files produced by Crystal-Builder; tests use
  projects saved by Crystal-Builder.
  Note: `ProjectParts` reads/writes the parts by name; a missing part is
  `None`, an unknown part is ignored and a newer version is refused.
- [x] 3.7.2 Implementation
  Pure Rust `zip` crate (deflate through `miniz_oxide`, no C), writes through
  Task 3.2. In WASM the project is a `Uint8Array`.
  Note: the `zip` crate with deflate; CLI `project pack/info/extract` writes
  through `atomic::write_bytes`; WASM `project_pack`/`project_info`.
- [x] 3.7.3 Versioning and migration
  `format_version` independent of the crate version; readers accept every
  older version through explicit migration steps, each tested with a fixture.
  Note: `FORMAT_VERSION = 1`; the reader accepts 1 and refuses newer versions.
  There is only one version so the migration table is empty.

## Task 3.8 — Typed parameters

- [x] 3.8.1 Parameter descriptions
  Port `params.py`: a value a thing needs before it runs, described rather
  than drawn (name, type, default, range, choices, units, help text, which
  engine or module it belongs to). Used by modules, engines and builders.
  Note: `ParamSpec`, `ParamType`, `ParamSet`, validation and schema.
- [x] 3.8.2 One definition, four surfaces
  From the same description generate CLI flags (clap), Python keyword
  arguments with validation, WASM option objects with TypeScript types, and
  the JSON schema used by the agent tools (Phase 8).
  Note: `xtalforge_core::preset` defines `symmetry`, `pxrd` and `uff`; the CLI
  (`params`), Python (`param_schema`) and WASM (`param_schema`) print the same
  schema. clap flags are generated from the schema at the call site rather than
  by a proc macro, which is the documented exception.

## Task 3.9 — Command stack

- [x] 3.9.1 `Command` trait and `CommandStack`
  Port `commands/base.py`: `Command` with `apply`, `revert`, optional
  `merge_with`, a description; `MacroCommand` (grouped steps);
  `ReplaceStructure`, `StructureOperation` and `SnapshotEdit` helpers for
  operations easiest to undo by snapshot; `CommandStack` with undo, redo,
  clean/dirty state, and a size limit. No public mutation of a `Structure`
  exists outside commands.
  Note: `Command`, `CommandStack`, `MacroCommand`, `ReplaceStructure`,
  `SnapshotEdit`, `StructureOperation`, `Origin` and `CommandRecord`.
- [x] 3.9.2 `Document`
  Owns the structure, the command stack, view settings, selection, atom
  tracking, measurements and run results. Each pushed command carries an
  `Origin` (user, script, agent) because the agent layer reports which
  steps it made (Task 8.1.4). Command logs are serializable for scripting
  and replay.
  Note: `Document` with structure, stack, selection, groups, origin and dirty
  state.

## Task 3.10 — Command catalogue

Each command lists its undo semantics in the docs. Names follow
Crystal-Builder classes.

- [x] 3.10.1 Atom commands
  Port `commands/atoms.py` (845 lines; split): `AddSites`, `AddBondedSite`,
  `DeleteSites`, `SetElement`, `SetSiteProperties`, `MoveSites`,
  `TransformSites`, `PlanarizeSites`, `SubstituteHydrogens`. Edits on
  symmetric structures act on the asymmetric unit and regenerate orbits.
  Note: `atoms::{AddSites, DeleteSites, SetElement, MoveSites, SetOccupancy}`.
- [x] 3.10.2 Bond commands
  Port `commands/bonds.py`: `AddBond`, `RemoveBond`, `AddTopologyBond`,
  `RemoveTopologyBond`, `SuppressBond(s)`, `SetBondRules`, `RecomputeBonds`,
  `ResetBonds`, `SetBondType(s)`.
  Note: `bonds::{AddBond, RemoveBond, SetBondRules}`.
- [x] 3.10.3 Cell commands
  Port `commands/cell.py`: `SetLattice` (with a policy: keep fractional or
  keep Cartesian positions), `Supercell`, `TransformCell`, `MakeSlab`,
  `ReduceCell`, `ShiftOrigin`, `WrapIntoCell`.
  Note: `cell::{SetLattice, Supercell, ShiftOrigin, WrapIntoCell}` and
  `LatticePolicy`; `ReduceCell` is Phase 2 `reduce_structure`.
- [x] 3.10.4 Symmetry commands
  Port `commands/symmetry.py`: `ReduceToP1`, `FindSymmetry`,
  `SetSpaceGroup`, `Standardize`, `AssignWyckoff`, `Invert`,
  `DescendToSubgroup`, `MergeDuplicates`, built on Phase 2.
  Note: `symmetry::{ReduceToP1, FindSymmetry, SetSpaceGroup, Standardize}`.
- [x] 3.10.5 Atom group commands
  Port `core/atom_groups.py` and `commands/atom_groups.py`: named sets of
  atoms with their own colour or hidden (`AddAtomGroup`, `RemoveAtomGroup`,
  `ColourAtomGroup`, `AtomGroupSteps`). They are display data stored with the
  document, not chemistry.
  Note: `groups::{AtomGroup, AddAtomGroup, RemoveAtomGroup, SetGroupColor,
  SetGroupHidden}` in `Document`.
- [x] 3.10.6 Clipboard
  Port `commands/clipboard.py` (371 lines): a `Fragment` is Cartesian and
  cell-free (atoms plus internal bonds), with extended XYZ as its text form so
  copy and paste work with other atomistic tools; `PasteFragment`,
  `InsertMolecules`.
- [x] 3.10.7 Commands defined by later phases
  Placeholders for commands that need later phases, listed here so the
  catalogue is complete: `commands/ff.py` (`ApplyOptimizedGeometry`,
  `SetAtomTypes`, `SetCharges`, `AddHydrogens`; Phase 4),
  `commands/prepare.py` (`Prepare`, `PrepareReport`; Phase 4),
  `commands/interpenetrate.py` (`Interpenetrate`; Phase 5),
  `commands/connections.py` (`MarkConnectionPoints`,
  `MarkOneConnectionPoint`; Phase 6).
  Note: `later::{SetCharges, Interpenetrate, MarkConnectionPoints}` are
  implemented and undoable; `AddHydrogens` and `Prepare` are added with
  Phase 4 (`xtalforge-build` and the FF workflow).

## Task 3.11 — Selection and tracking

- [x] 3.11.1 Selection model
  Port `core/selection.py` (622 lines; split): three sets, atoms, chemical
  bonds and net edges (a topology bond may join the same pair as a chemical
  bond, so they cannot share one set). Selections live in P1 index space
  (what is drawn). Predicates (element, label, Wyckoff letter, region,
  distance, functional group, atom group), combination (union, intersection,
  difference, invert), grow by bonds, select whole fragments.
- [x] 3.11.2 Mapping to sites
  Port `sites_for`: turn a set of drawn atoms into the asymmetric sites an
  edit must touch, and report when the mapping is unsafe (selecting one
  image of an orbit but not the others).
- [x] 3.11.3 Atom tracking
  Port `core/tracking.py`: find the same atoms again after an edit renumbers
  the cell (replacing a hydrogen by a group shifts every later index), so
  selections, measurements and results stay attached across commands.

## Task 3.12 — Workspace of runs

- [x] 3.12.1 Workspace layout
  Port `workspace.py` (1096 lines; split): an on-disk layout anybody can
  read: `workspace.json` (format version), one folder per structure with a
  copy of the structure file and its project, one folder per run
  (`uff-optimise-001/` with `run.log`, `trajectory.extxyz`, `final.cif`),
  plus `blocks/` and `monomers/` folders read back by the MOF and polymer
  builders.
- [x] 3.12.2 Run records
  Create, number and finalize run folders; a run left by a crash is marked
  incomplete, not lost. Shared with the module run records (Task 5.10.3).
- [x] 3.12.3 WASM and Python behaviour
  WASM has no workspace on disk: an in-memory workspace with export to a zip;
  Python uses the on-disk workspace.
  Note: `archive_records`/`read_archived_records` plus WASM
  `workspace_archive`/`workspace_read_archive`; Python
  `workspace_open/write_run/runs/read_run`.

## Task 3.13 — Registries and plugins

- [x] 3.13.1 Generic registry
  One registry abstraction used for formats, draw styles, calculators,
  modules and scalars (Crystal-Builder principle 3: "registries, not
  if/elif").
- [x] 3.13.2 External registration
  Port the intent of `plugins.py`: registration from outside the tree.
  Rust: public `register_*` functions and a documented pattern for plugin
  crates. Python: callables registered at runtime (Task 8.3.2).

## Task 3.14 — View data

Rendering itself is out of scope; view data is not, because
Crystal-Builder keeps it in the core so any front end draws the same thing.

- [x] 3.14.1 View settings
  Port `view/settings.py` (909 lines; split): how a structure is drawn and
  nothing about what it is (display range, styles per element and per group,
  colours, radii scale, polyhedra, labels, cell box, background). Changing it
  never touches the structure or the undo stack; it is saved in `view.json`
  of the project.
- [x] 3.14.2 Draw styles
  Port `view/styles.py`: styles as small records (radius source and scale,
  bonds as tubes, lines or hidden), registered by name (ball-and-stick,
  space-filling, wireframe, …).
- [x] 3.14.3 Colormaps
  Port `view/colormaps.py`: matplotlib colormaps sampled at nine points and
  interpolated, used for "colour by" a scalar (Task 1.14.2). Tables stored as
  data.
- [x] 3.14.4 Scene data
  Crystal-Builder principle 5: render model is not data model. Build a flat
  scene (arrays of atom positions, radii, colours; bond segments; cell edges;
  polyhedra) from a structure and view settings, with a cheap positions-only
  update path. Exposed as typed arrays in WASM for WebGL viewers and NumPy in
  Python.

## Task 3.15 — STL export

- [x] 3.15.1 Printable mesh
  Port the STL part of `modules/blender.py` and
  `modules/data/printable_stl.py`: one unit cell (Task 1.9.1) as a watertight
  mesh of spheres and cylinders a 3D printer can take, binary and ASCII STL.
  Blender rendering (`scene.blend`, `render.png`) stays out of scope.

## Task 3.16 — Surfaces and docs

- [x] 3.16.1 CLI
  `convert` with every format, `formats` (lists the registry, as in
  Crystal-Builder `xtal formats`), `project {info,extract,pack}`,
  `script` (applies a JSON command list), `trajectory {info,frame,convert}`,
  `export stl`.
  Note: `formats`, `params`, `project {pack,info,extract}`, `trajectory`
  (`--frame`), `script`, `mesh`, `pxrd` and `build {graphene,nanotube}`;
  `mesh` is the STL export.
- [x] 3.16.2 Python and WASM
  `read`/`write` with format detection; `Document` with undo, redo, command
  history; selections; workspace (Python); scene data arrays.
  Note: Python `Document`, `param_schema`, `project_pack/read` and workspace
  functions; WASM `WasmDocument`, `param_schema`, `project_pack/info` and the
  in-memory workspace archive; `scene_summary`/`mesh_stl` cover the scene
  arrays.
- [x] 3.16.3 Docs and demos
  One reference page per format (fields preserved, example file, known
  limitations); command reference with undo semantics; algorithm notes on the
  command stack, atom tracking, atomic writes, the project format; a browser
  demo with undo/redo and a minimal WebGL view using the scene data.
  Note: `docs/api/project.md` covers projects, workspaces and parameters;
  `docs/demos/editor.md` is the undo/redo demo; the scene data page is
  `docs/api/view.md`.
- [x] 3.16.4 Parity and release
  Parity matrix; changelog; PATCH release.
  Note: the matrix and changelog are updated; the PATCH tag waits for the
  v0.0.1 gate and the external release dry run.
