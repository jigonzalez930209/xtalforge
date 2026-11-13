# Cargo features

The facade crate `xtalforge` is the crate applications depend on. Its features
are additive: enabling one never changes results, only which optional
capabilities are compiled. Every domain crate carries the same switches so a
downstream crate can be built directly with or without them.

| Feature | Default | Domains | Effect |
|---|---|---|---|
| `parallel` | yes | core, io, symmetry | Reserves Rayon for the parallel algorithms. The current algorithms are single-threaded and deterministic; the feature exists so builds can already opt out with `--no-default-features`. |
| `sci-form` | no | force fields | Links the published MIT crate `sci-form` 2.0.0 for in-process GFN0/GFN1/GFN2 (Task 4.10). The methods are molecular only; GFN2 returns forces and the registry lists `gfn0`, `gfn1` and `gfn2`. |
| `ml` | no | force fields | Reserved for machine-learned potentials (Task 8.4). The decision is recorded in `docs/algorithms/ml-backend.md`; currently a no-op. |
| `external` | no | force fields | Reserved; the runtime engine probes and input writers for DFTB+, xTB and ORCA are always compiled and need no feature (Task 4.9). |
| `net` | no | remote data | Adds the blocking `ureq` HTTP client for OPTIMADE, COD and Materials Project fetches (Task 13.3). Parsing and WASM `fetch` do not need it. |
| `mof-data` | no | MOF builder | Reserved for a future full MOF bundle; the curated six-topology/eight-block subset is always embedded (Tasks 6.11-6.16). |

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
5. The default WASM package has a 5120 KiB budget, checked by
   `scripts/check-wasm-size.sh`. The hello-world build was about 41 KiB; the
   Phase 1 readers and builders brought it to about 438 KiB, symmetry to about
   813 KiB, and the Phase 3-7 surfaces (including the polymer/embedding code)
   and the agent session took the release package to 3943 KiB at 1.1.0.
   The Phase 9 core features took it to 4603 KiB at 1.2.0 and the MOF
   platform plus usability work to 4964 KiB at 1.4.0, inside the 5120 KiB
   budget.

## Threaded WebAssembly

The default `xtalforge-wasm` package is single-threaded. A threaded build
needs a shared `WebAssembly.Memory`, `SharedArrayBuffer` and cross-origin
isolation headers from the host page; without them the module cannot start.
The Rust `parallel` feature is not enabled in the published WASM build, and
enabling it never changes a result: every parallel reduction uses a fixed
order, so the single-threaded and threaded paths agree.
