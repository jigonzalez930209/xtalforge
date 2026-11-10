# Algorithm coverage

Every ported algorithm has a page with its formulation, references, validation
data and differences from the reference. This page is the index that proves it:
one row per crate module, the page that documents it, and the test that
validates it.

| Module | Page | Validation |
|---|---|---|
| `core::lattice` | [Lattice](./lattice) | crystal-system cases, volume and angles |
| `core::limits` | [Limits](./limits) | `limits_audit.rs` at and beyond every budget |
| `core::p1` | [P1 expansion](./p1) | provenance and operation-orbit golden |
| `core::supercell`, `cell` | [Supercells](./supercell) | duplicate-free expansions, determinism |
| `core::slab` | [Slabs](./slab) | layer counts and vacuum |
| `core::transforms` | [Transforms](./transforms) | Cartesian/frac round-trips, ellipsoids |
| `core::neighbors` | [Neighbor search](./neighbors) | linked-cell scaling and exact lists |
| `core::bonding` | [Bond perception](./bonding) | corpus bonds and CoRE MOF |
| `core::rings` | [Primitive rings](./rings) | diamond six-rings, corpus census |
| `core::groups` | [Functional groups](./groups) | fragment membership |
| `core::measure` | [Measures](./measures) | distances, angles, torsions |
| `core::scalars` | [Measures](./measures) | built-ins and plugin registry |
| `io::cif` | [CIF parsing](./cif) | corpus, CoRE MOF round-trip |
| `io::poscar`, `xyz`, `gen`, `lammps`, `shelx`, `pmg`, `cube`, `cssr` | [POSCAR and XYZ](./formats) | corpus and mutation robustness |
| `io::project` | [Projects](../api/project) | pack/extract round-trip |
| `io::trajectory` | [Formats](./formats) | lazy frames and counts |
| `symmetry::dataset` | [Symmetry detection](./symmetry) | spglib/moyo audit, corpus |
| `symmetry::backend` | [moyo backend](./symmetry-backend) | audit and fallbacks |
| `symmetry::reduction` | [Lattice reduction](./symmetry) | Niggli/Delaunay/Minkowski |
| `symmetry::subgroups` | [Symmetry](./symmetry) | subgroup enumeration and lowering |
| `symmetry::brillouin`, `kpath` | [Symmetry](./symmetry) | Bravais/k-path golden |
| `ff::uff` | [UFF and UFF4MOF](./uff) | sci-form and tmQM |
| `ff::charges` | [Charges and electrostatics](./charges) | EQeq authors' MOF-5 |
| `ff::optimize`, `scan`, `eos` | [Optimization, scans and EOS](./optimization) | analytic minima and scans |
| `ff::vibrations` | [Vibrations and thermodynamics](./vibrations) | analytic Hessian |
| `ff::prepare` | [Preparation](./prepare) | pipeline steps and reports |
| `ff::models` | [ML potentials](../api/ml) | cache fetch/verify/corrupt |
| `ff::external`, `quantum` | [External engines](../api/engines) | recorded DFTB+ output |
| `analysis::pxrd` | [Analysis](./analysis) | rutile and CoRE MOF |
| `analysis::peaks`, `profiles` | [Analysis](./analysis) | synthetic peaks |
| `analysis::grid` | [Analysis](./analysis) | HKUST-1, ZIF-8 bounds |
| `analysis::isosurface` | [Analysis](./analysis) | closed periodic mesh |
| `analysis::nets`, `rcsr`, `merge`, `overlay`, `pareto` | [Nets and RCSR](./nets) | full RCSR index build |
| `analysis::zeo` | [Analysis](./analysis) | recorded Zeo++ output |
| `analysis::indexing` | [Powder indexing](./indexing) | golden M20/F_N |
| `analysis::lebail`, `pawley` | [Le Bail](./lebail), [Pawley](./pawley) | golden fits and esds |
| `analysis::rietveld`, `refine` | [Rietveld](./rietveld) | rutile golden |
| `analysis::extinction`, `automatic` | [Extinction](./extinction), [Powder pipeline](./powder) | class partition and ranking |
| `analysis::cells`, `instrument` | [Powder pipeline](./powder) | free parameters, Caglioti |
| `build::smiles`, `embed`, `fragments`, `substitute` | [Builders](./builder) | SMILES, fragments, substitution |
| `build::hydrogens` | [Builders](./builder) | valence and periodic placement |
| `build::polymer` | [Builders](./builder) | density and determinism |
| `build::carbon`, `ztc`, `ribbon`, `twist` | [Disordered carbons](./carbon) | carbon golden distributions |
| `build::fullerene`, `spiral` | [Builders](./builder) | C20/C60 vertex and edge counts |
| `build::net`, `guest` | [Builders](./builder) | CGD nets and clearance |
| `mof::topology`, `block`, `orient`, `build`, `advanced` | [MOF builder](./mof) | curated subset and determinism |
| `commands` | [Commands](../api/commands) | undo/redo and origin |
| `modules::registry`, `jobs`, `workspace`, `report` | [Modules](./modules) | FIFO, records, archive, report |
| `agent::*` | [Agent session](./agent) | 8.1-8.2 test suite |
| `view::scene`, `mesh` | [View data](../api/view) | scene and STL |
| `models` | [ML potentials](../api/ml) | MACE golden data |

## Known subsets

These are documented subsets, not missing pages:

- The MOF builder ships a curated topology/block set, not the full PORMAKE
  database.
- The RCSR index has 22 unkeyed entries of 2,931.
- The phonon q-mesh and preparation-step list are the subsets recorded in
  `docs/roadmap/index.md`.
- ORB-v3 and MatterSim run through ASE but have no recorded golden comparison;
  MACE does.
- The external engine launchers are tested with simulated executables; the
  in-process xTB methods use `sci-form`.
