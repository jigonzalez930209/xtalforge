# API stability and SemVer policy

xtalforge follows Semantic Versioning for one version shared by every crate
and all four surfaces. This page restates the rules a 1.0 consumer can rely on;
the release procedure is in `AGENTS.md` section 5.

## What is public API

A change to any of these is a breaking change:

- public items of the facade crate `xtalforge` (names, signatures, types);
- CLI subcommand names, flags, exit codes and the shape of `--json` output;
- Python module, class and function names, keyword arguments and exception
  types;
- JS/TS exported names and TypeScript types of `xtalforge-wasm`;
- error codes (`Error::code()` strings);
- `schema_version` of serialized structures and `format_version` of project
  files.

A change that makes a function return what its documentation always promised
(a wrong geometry, a rejected valid input) is a fix, not a breaking change, even
when the output of the faulty version differed. The 1.0.0-rc.1 stabilization
period used this rule for the graphene, nanotube and nanoribbon builders and
for `structure_info` (see the changelog).

Private modules, test fixtures, benchmark output, log text and the exact
wording of a human-readable message are not public API.

## Version stages

| Stage | Version | Rule |
|---|---|---|
| Pre-stable | `0.y.z` | MINOR for a breaking change on any surface or a finished phase; PATCH for compatible features, fixes and data |
| Stable | `1.0.0` and later | MAJOR for a breaking change; MINOR for a compatible capability; PATCH for a fix |

Independent numbers:

- `schema_version`: bumped on any change to a serialized type;
- `format_version`: bumped on any change to the project archive;
- data bundles: `data-YYYY.MM.DD`.

## Deprecation

1. Announce the replacement in the changelog, the rustdoc and this page.
2. Keep the old item for at least one MINOR release and mark it
   `#[deprecated(since = "...", note = "...")]` (Rust) or the equivalent in
   the bindings.
3. Remove it in a MAJOR release, never in a PATCH.

The repository does not silence deprecation warnings; a deprecated item still
compiles and still has tests.

## Freeze test

`crates/xtalforge/tests/api_freeze.rs` pins representative signatures with
compile-time function-pointer types: the core hello-world, CIF read/write,
scene data, agent diagnostics, the batch runner, the verb catalogue, an answer
and three error codes. Adding an argument, renaming a function or changing a
return type fails that test. It is deliberately a sample, not every item; the
rustdoc and the parity matrix are the full inventory.

## Checklist before a breaking change

- [ ] Decide the version with the table above.
- [ ] Update the four surfaces in the same change.
- [ ] Update `docs/api/parity.md` and the rustdoc.
- [ ] Add a `### Breaking` entry to `CHANGELOG.md`.
- [ ] Bump `schema_version` or `format_version` if a serialized type changed.
- [ ] Run the local gate from `AGENTS.md` section 5.4 on the exact commit.

## Current status

`1.1.0` is published (2026-11-11) on crates.io, PyPI, npm, the GitHub release
and GitHub Pages. It is the first feature release after the stable `1.0.0`
(2026-11-10), which followed the `1.0.0-rc.1` stabilization period; its new
fields and WASM exports are backward compatible. The public API is frozen, the
parity matrix is closed, and later changes follow the table above. `0.0.3` was
the last published pre-1.0 release.
