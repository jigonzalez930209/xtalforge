# Introduction

`xtalforge` is a pure Rust library for crystal structures, exposed through the
Rust API, a command-line binary, Python bindings and WebAssembly bindings.
Phases 0-7 are complete locally and Phase 8 (agent and distribution) is the
remaining roadmap work; see the [roadmap](/roadmap/) for the current status.

Start here:

- [Quick start](./quickstart) for Rust, the CLI, Python and WebAssembly
- [Concepts](./concepts): lattices, sites, fractional and Cartesian
  coordinates, symmetry and bonding
- [File formats](./formats) and the format registry
- [Cargo features](./features) and the WASM budget
- [Data bundles](./data) and the embedded-versus-external policy
- [Contributing](./contributing) and the local check list

The [API reference](/api/) documents each capability on all four surfaces, the
[algorithm pages](/algorithms/) record the formulations and validation, and
the [demos](/demos/) run the WebAssembly package in the browser.
