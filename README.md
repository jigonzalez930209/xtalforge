# xtalforge

A pure Rust library for building, manipulating, analysing and exporting crystal
structures, with a CLI, Python bindings and WebAssembly bindings developed at
parity with the Rust API.

> Status: Phase 0 complete locally (`0.0.0`, nothing published). The workspace,
> rule tooling, hello-world in the four surfaces, corpus and golden data, docs
> site and dependency map are in place. GitHub Actions remain disabled until
> the maintainer enables them. The plan is in
> [`docs/roadmap`](docs/roadmap/index.md); the project rules are in
> [`AGENTS.md`](AGENTS.md).

## Install

Nothing is published yet. After the first release the commands will be:

```bash
cargo add xtalforge@0.0.0
pip install xtalforge==0.0.0
npm install xtalforge-wasm@0.0.0
```

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
