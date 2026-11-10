# Introduction

`xtalforge` is a pure Rust library for crystal structures, exposed through the
Rust API, a command-line binary, Python bindings and WebAssembly bindings.
Phases 0-8 are implemented; `1.0.0-rc.1` is published, with the public API
frozen and the parity matrix closed. See the [roadmap](/roadmap/) for the
current status.

Start here:

- [Quick start](./quickstart) for Rust, the CLI, Python and WebAssembly
- [Concepts](./concepts): lattices, sites, fractional and Cartesian
  coordinates, symmetry and bonding
- [File formats](./formats) and the format registry
- [Automation guide](./automation): MCP server, verbs and skill documents
- [FAQ](./faq): scope, accuracy, limits, licensing and the questions a user asks first
- [Cargo features](./features) and the WASM budget
- [Data bundles](./data) and the embedded-versus-external policy
- [Contributing](./contributing) and the local check list
- [Release and deployment](./release): tokens, tag flow, CLI packages and
  caches

The [API reference](/api/) documents each capability on all four surfaces, the
[algorithm pages](/algorithms/) record the formulations and validation, and
the [demos](/demos/) run the WebAssembly package in the browser.
