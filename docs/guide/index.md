# Introduction

`xtalforge` is a pure Rust library for crystal structures, exposed through the
Rust API, a command-line binary, Python bindings and WebAssembly bindings.
Phases 0-15 are complete; `1.8.0` is published, with the public API frozen and
the parity matrix closed. See the [roadmap](/roadmap/) for the current status.

Start here:

- [Quick start](./quickstart) for Rust, the CLI, Python and WebAssembly
- [Performance](./performance): the 1.8.0 comparative benchmarks, methodology
  and regression gate
- [Concepts](./concepts): lattices, sites, fractional and Cartesian
  coordinates, symmetry and bonding
- [File formats](./formats) and the format registry
- [Migrating from ASE](./from-ase), [pymatgen](./from-pymatgen) and
  [spglib](./from-spglib)
- [Automation guide](./automation): MCP server, verbs and skill documents
- [Python plugins](./plugins): entry points and the template repository
- [FAQ](./faq): scope, accuracy, limits, licensing and the questions a user asks first
- [Cargo features](./features) and the WASM budget
- [Data bundles](./data) and the embedded-versus-external policy
- [Contributing](./contributing) and the local check list
- [Release and deployment](./release): tokens, tag flow, CLI packages and
  caches

The [API reference](/api/) documents each capability on all four surfaces, the
[algorithm pages](/algorithms/) record the formulations and validation, and
the [demos](/demos/) run the WebAssembly package in the browser.
