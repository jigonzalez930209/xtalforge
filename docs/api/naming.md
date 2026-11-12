# Option naming across the surfaces

The 1.x API cannot rename options, so this page records the canonical name of
each option and the aliases accepted on each surface. Removals are scheduled
for 2.0 (Task 17.1). The deprecation mechanism is documented in
[`docs/guide/deprecation.md`](../guide/deprecation.md).

| Concept | Rust | CLI | Python | WASM | Alias / note |
|---|---|---|---|---|---|
| Symmetry tolerance (Angstrom) | `SymmetryOptions::symprec` | `--symprec` | `symprec=` | `symprec` | `tolerance` is the spelling in the docs; `symprec` is kept |
| Angle tolerance (degrees) | `SymmetryOptions::angle_tolerance` (radians in Rust) | `--angle-tolerance` (degrees) | `angle_tolerance=` (radians, Python convention) | `angle_tolerance` | CLI converts degrees at the boundary |
| Matcher length tolerance | `MatcherOptions::ltol` | `--ltol` | `ltol=` | `ltol` in options JSON | relative length |
| Matcher site tolerance | `MatcherOptions::stol` | `--stol` | `stol=` | `stol` | Angstrom |
| Matcher angle tolerance | `MatcherOptions::angle_tol` | `--angle-tol` | `angle_tol=` | `angle_tol` | degrees |
| Grid spacing | `PoreOptions::spacing` | `--spacing` | `spacing=` | `spacing` | Angstrom |
| Probe radius | `PoreOptions::probe_radius` | `--probe` | `probe=` | `probe` | Angstrom |
| Number of workers | `ParallelOptions::threads` | `--threads` | `set_threads(n)` | n/a (browser threads) | process-wide in Python |
| Format hint | reader `hint` argument | extension or `-` stdin | `format_hint=` | `hint` | content sniffing when empty |
| Random seed | `seed` | `--seed` | `seed=` | `seed` in options JSON | required for determinism |
| Temperature (K) | `WidomOptions::temperature` | `--temperature` | `temperature=` | `temperature` | adsorption only |
| Pressure (bar) | `GcmcOptions::pressures` | `--pressures` | `pressures=` | `pressures` | adsorption only |

## Deprecated aliases

No alias is removed in 1.x. The following aliases are candidates for removal
in 2.0 and are accepted now:

- Python: `interpolate_structures` (preferred: `interpolate`).
- Python: `mof_read_database` and `mof_screen` (the Rust names are
  `read_core_mof` and `screen_database`; the Python `mof_` prefix is kept).
- WASM: every `*_info` name (the object-returning form is preferred for new
  code; the old names remain).
- CLI: `--angle-tol` and `--angle-tolerance` are synonyms in the matcher.

## Rules for new options

1. One concept, one canonical name per surface; the table above is the source
   of truth.
2. A new surface keeps the Rust name unless it conflicts with the host
   language convention (for example, Python exposes radians while the CLI and
   Rust boundary use radians internally and degrees at the CLI flag).
3. Any alias must be additive, documented here, and must warn once through
   the deprecation mechanism.
