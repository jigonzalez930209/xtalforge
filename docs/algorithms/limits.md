# Limits model

## Problem

A supercell, a porosity grid or a Hessian is expensive before a single atom or
grid point exists. The estimate must be arithmetic and happen before allocation
so the caller can be warned or refused with a concrete alternative.

## Model

![Limits model](/diagrams/limits-model.svg)


Every budget has a soft and a hard limit. Standard is the 8 GB-class table;
Generous doubles both; WarnOnly keeps the soft limit and removes the hard one.
The check returns `Verdict { level, estimate, sentence }`:

- `ok` — estimate at or below the soft limit, empty sentence;
- `warn` — over soft and (if present) at or below hard;
- `refuse` — over hard.

The sentence always names the estimate and what would fit (a smaller cube,
a coarser grid, or fewer atoms).

## Estimates

- Drawn atoms: `n_atoms * product(range_hi - range_lo)`.
- Supercell: `n_atoms * |det P|`.
- Grid: `round(length / spacing)` per axis, minimum 8, product.
- Phonons: the supercell estimate.
- Analytic Hessian: the free atom count, with the memory note
  `13 MiB * n_atoms`.

## Validation

The table and profiles are covered by `crates/xtalforge-core/tests/limits.rs`.
The surfaces are checked by the Python and WASM tests. The carbon, nanotube,
graphene, phonon and Hessian budgets are fixed in the table; a caller decides
which generic check (`check_supercell`, `check_grid`, `check_phonons`,
`check_hessian`) matches its allocation.
