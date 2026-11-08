# xtalforge

A pure Rust library for building, manipulating, analysing and exporting crystal
structures, with a CLI, Python bindings and WebAssembly bindings developed at
parity with the Rust API.

> Status: bootstrap (`0.0.0`). The Cargo workspace, lint configuration and
> contributor documentation exist; nothing is released yet. The plan is in
> [`docs/roadmap`](docs/roadmap/index.md); the project rules are in
> [`AGENTS.md`](AGENTS.md).

## Goals

- Everything the headless core of [Crystal-Builder](https://github.com/JulesOpp/Crystal-Builder)
  does: symmetry-aware structure editing, file formats, force fields and
  optimizers, PXRD, porosity, topology, builders (molecules, carbons,
  polymers, MOFs) and powder refinement.
- Symmetry through [moyo](https://github.com/spglib/moyo).
- Proven code reused from [sci-form](https://github.com/jigonzalez930209/sci-form).
- Four surfaces at parity: Rust (`xtalforge` crate), CLI (`xtalforge`
  binary), Python (`xtalforge` on PyPI), WebAssembly (`xtalforge-wasm` on
  npm).
- Documentation site built with VitePress (`docs/`).

## Working on the project

```bash
scripts/check-rules.sh                 # 500-line limit and English-only checks
cargo test --workspace
cargo clippy --workspace --all-targets -- -D warnings
cd docs && npm install && npm run docs:dev
```

See [`docs/guide/contributing.md`](docs/guide/contributing.md) for the full
checklist and the roadmap workflow.

## License

Licensed under either of

- Apache License, Version 2.0 ([`LICENSE-APACHE`](LICENSE-APACHE))
- MIT license ([`LICENSE-MIT`](LICENSE-MIT))

at your option. Third-party code and data: `THIRD_PARTY.md`.
