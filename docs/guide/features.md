# Cargo features

The facade crate `xtalforge` is the crate applications depend on. Its features
are additive: enabling one never changes results, only which optional
capabilities are compiled. Every domain crate carries the same switches so a
downstream crate can be built directly with or without them.

| Feature | Default | Domains | Effect |
|---|---|---|---|
| `parallel` | yes | all | Enables optional Rayon support. No Phase 0/1 algorithm uses it yet; it exists so builds can already opt out with `--no-default-features`. |
| `sci-form` | no | force fields | Reserved for quantum-chemistry methods through the `sci-form` crate (Task 4.10). Currently a no-op. |
| `ml` | no | force fields | Reserved for machine-learned potentials (Task 8.4). Currently a no-op. |
| `external` | no | force fields | Reserved for subprocess engines such as DFTB+, xtb and ORCA (Task 4.9). Currently a no-op. |
| `mof-data` | no | MOF builder | Reserved for the embedded MOF data bundle (Task 6.11.3). Currently a no-op. |

## Building with and without defaults

```bash
cargo build -p xtalforge                     # default features (parallel)
cargo build -p xtalforge --no-default-features
cargo build -p xtalforge --features sci-form
cargo build -p xtalforge --all-features
```

## Rules for new features

1. Default features must depend only on pure Rust crates (AGENTS.md rule 3).
   Anything that needs a native library, a subprocess or a large embedded data
   set is opt-in.
2. A feature is documented in the facade crate docs and in this page in the
   same commit that introduces it.
3. The WASM build keeps heavy and external features disabled by default; the
   parity matrix records the exception when a capability cannot exist in the
   browser.
4. Features stay additive. Two features enabled together must behave like the
   union of their capabilities, never like a different code path that changes
   numeric output.
5. The default WASM package has a 512 KiB budget, checked by
   `scripts/check-wasm-size.sh`. The hello-world build was about 41 KiB; the
   Phase 1 readers, builders and analysis brought it to about 424 KiB, still
   inside the budget.
