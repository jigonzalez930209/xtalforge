# Algorithms

One page per ported algorithm. Template for every page:

1. Problem statement
2. Mathematical formulation
3. References (papers, upstream code)
4. Pseudocode
5. Complexity and numerical caveats
6. Validation data and tolerances

Diagrams live in `docs/public/diagrams/` (29 SVG, one shared style, animated
with SMIL and static under `prefers-reduced-motion`). Every page that describes
a procedure or a geometry embeds one. The inventory and decision pages
([Golden data](./golden-data), [External validation](./external-validation),
[Algorithm coverage](./coverage), [Dependency map](./dependency-map),
[sci-form reuse](./sci-form-reuse), [ML backend](./ml-backend),
[Modules](./modules), [Symmetry backend](./symmetry-backend)) are tables and
records of decisions and have no algorithm to draw.

Pages are added with the feature that introduces them:

- [Lattice conventions](./lattice) (Task 1.2.5)
- [Limits model](./limits) (Task 1.3)
- [P1 expansion with provenance](./p1) (Task 1.6)
- [Supercells and origin shifts](./supercell) (Task 1.7)
- [Slab construction](./slab) (Task 1.8)
- [Cartesian transforms and ellipsoids](./transforms) (Task 1.9)
- [Periodic neighbor search](./neighbors) (Task 1.10)
- [Bond perception](./bonding) (Task 1.11)
- [Primitive rings](./rings) (Task 1.12)
- [Functional groups](./groups) (Task 1.13)
- [Measures](./measures) (Task 1.14)
- [CIF parsing](./cif) (Task 1.15)
- [POSCAR and XYZ](./formats) (Task 1.16)
- [Symmetry detection](./symmetry) (Task 2.1)
- [Agent session and diagnostics](./agent) (Task 8.1)
- [Analysis: patterns, peaks, grids and nets](./analysis) (Task 5.1-5.9)
- [Builders: molecules, polymers and carbon](./builder) (Task 6.1-6.9)
- [MOF builder](./mof) (Task 6.12-6.16)
- [Powder pipeline](./powder) (Task 7.5-7.10)
- [Modules, jobs and reports](./modules) (Task 3.13, 5.10)
- [Algorithm coverage](./coverage) (Task 8.7.2)
- [UFF and UFF4MOF](./uff) (Task 4.2)
- [Charges and electrostatics](./charges) (Task 4.3)
- [Optimization, scans and EOS](./optimization) (Task 4.4-4.5)
- [Vibrations and thermodynamics](./vibrations) (Task 4.8)
- [Powder indexing](./indexing) (Task 7.5)
- [Pawley fit](./pawley) (Task 7.6)
- [Le Bail extraction](./lebail) (Task 7.6)
- [Rietveld refinement](./rietveld) (Task 7.8)
- [External validation](./external-validation) (Task 0.9)
- [Nets and RCSR identification](./nets) (Task 5.6)
- [Extinction screening](./extinction) (Task 7.10)
- [Machine-learning backend](./ml-backend) (Task 4.11)
- [Symmetry backend: moyo](./symmetry-backend) (Task 0.3.5)
- [Validation and tolerances](./validation) (Task 0.9.4)
- [Golden data](./golden-data) (Task 0.9.3)
- [sci-form reuse](./sci-form-reuse) (Task 0.4)
- [Disordered carbons](./carbon) (Task 6.8)
- [Dependency replacement map](./dependency-map) (Task 0.10)

The neighbor-search page records the linked-cell index and the scaling
benchmark in Phase 1.
