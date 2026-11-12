# Frequently asked questions

## What is xtalforge?

A pure Rust library for building, manipulating, analysing and exporting crystal
structures, with a CLI, Python and WebAssembly bindings. It ports the headless
core of Crystal-Builder; the PySide6/VTK application is out of scope.

## Is it a drop-in replacement for Crystal-Builder?

For the headless core, the behaviour and the numerical results are the target,
not the Python API. Names are Rust-idiomatic and the CLI is command-oriented.
The [migration guide](./migration) maps every module, command and agent verb
and lists the behavioural differences.

## Which parts are pure Rust?

The default features depend only on pure Rust crates. `moyo` supplies symmetry,
`nalgebra` linear algebra, and the parsers, force fields, analysis and builders
are native code. External programs (DFTB+, xTB, ORCA, Zeo++) are optional
probes behind a runtime check; MACE/ORB/MatterSim run through the Python bridge,
not native Rust.

## Why moyo for symmetry?

It is a maintained Rust implementation of the spglib algorithms with the
space-group catalogue, Wyckoff tables and subgroup enumeration. It sits behind
our own types, so the backend can be replaced; the audit and the package
boundary are documented in
[the symmetry backend page](../algorithms/symmetry-backend).

## How accurate are the numbers?

Each algorithm page states its formulation and its validation. The strongest
external checks are:

- EQeq against the authors' program on MOF-5 (maximum difference `1.6e-4 e`);
- UFF against sci-form on small molecules (within `1e-4 kcal/mol`) and over 18
  tmQM metals for typing and finiteness;
- PXRD/Pawley/Rietveld against the Crystal-Builder/rietx rutile data;
- CoRE MOF 2019 and tmQM 2024 subsets for parsing, formula, cell volume,
  coordination and porosity bounds;
- MACE-MP-0 energies, forces and stress within `1e-6 eV`, `1e-5 eV/A` and
  `1e-6 GPa`.

The full matrix is in
[external validation](../algorithms/external-validation).

## Can it replace a DFT code?

No. UFF, QEq/EQeq/EEQ, Ewald and the ML bridge are classical or fitted models.
They do not give electronic structure, band gaps or reaction barriers. The
quantum adapter writes inputs and parses outputs for external codes; it is not
a solver.

## Does it run in the browser?

Yes, through `xtalforge-wasm`. The browser has no filesystem and no
subprocesses: string in, string/JSON out, and the in-process methods only. The
release package is about 4603 KiB (budget 4736 KiB; `scripts/check-wasm-size.sh` is the source of truth). The parity matrix records every
`n/a`.

## Are the results deterministic?

Yes. The same input and options give the same output independent of thread
count and platform. Random algorithms take an explicit seed. The determinism
tests run the same graph, neighbor list and supercell repeatedly and in four
threads.

## What is the largest structure it handles?

The standard limits are 100,000 drawn atoms, 50,000 supercell atoms, a
1.5M-point porosity grid and smaller budgets for carbon, phonons and Hessians.
The limits check refuses before allocating; `Verdict` explains what would fit.
Use `--limits generous` or `warn-only` deliberately.

## How do I know a CIF was read correctly?

Run `info` and `diagnostics`: formula, cell, density, charge, occupancy,
overlaps, coordination and symmetry. The CoRE MOF test checks the archive's own
`_chemical_formula_sum` and `_cell_volume` against the parse for 18 frameworks.

## Does it write project files?

Yes. `.xtalproj` stores the structure, cell, bonds, site properties, view
settings and session data. `format_version` is independent of the crate
version; a change bumps it.

## Can I extend it?

Yes. Rust plugins register formats, calculators, modules, scalars and style
presets; `examples/plugin` is a complete crate. Python callables register as
modules, calculators and formats through `xtalforge.Plugins`. See
[plugins](../api/plugins).

## Why is the agent layer part of the library?

An AI assistant needs a stable verb catalog, typed answers, diagnostics and a
way to review its edits. Putting that in the library (rather than in a script)
makes the CLI, Python, WASM and MCP server share one list and one history. The
session applies edits with `Origin::Agent`, so a user can undo a burst.

## What is out of scope?

The GUI, the VTK window proxy, the full PORMAKE database (a curated subset
ships), a native Rust ML backend, and a full DFT solver. These are recorded as
exceptions or known subsets, not hidden gaps.

## How is it licensed?

`MIT OR Apache-2.0`. Every third-party data set and ported code is recorded in
[`THIRD_PARTY.md`](https://github.com/jigonzalez930209/xtalforge/blob/main/THIRD_PARTY.md)
with its license and citation.
