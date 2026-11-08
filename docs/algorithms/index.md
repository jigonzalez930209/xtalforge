# Algorithms

One page per ported algorithm. Template for every page:

1. Problem statement
2. Mathematical formulation
3. References (papers, upstream code)
4. Pseudocode
5. Complexity and numerical caveats
6. Validation data and tolerances

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
- [UFF and UFF4MOF](./uff) (Task 4.2)
- [Charges and electrostatics](./charges) (Task 4.3)
- [Optimization, scans and EOS](./optimization) (Task 4.4-4.5)
- [Vibrations and thermodynamics](./vibrations) (Task 4.8)
- [Machine-learning backend](./ml-backend) (Task 4.11)
- [Symmetry backend: moyo](./symmetry-backend) (Task 0.3.5)
- [Validation and tolerances](./validation) (Task 0.9.4)
- [Golden data](./golden-data) (Task 0.9.3)
- [sci-form reuse](./sci-form-reuse) (Task 0.4)
- [Dependency replacement map](./dependency-map) (Task 0.10)

The neighbor-search page records the direct algorithm; the linked-cell
acceleration and the criterion benchmarks are still listed as remaining work
in Phase 1.
