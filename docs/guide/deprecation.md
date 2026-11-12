# Deprecation policy

`xtalforge` follows the stability policy in `docs/api/stability.md`. During
1.x an old name is never removed; it keeps working and points at its
replacement. The 2.0 migration guide collects the removals.

## Mechanism per surface

- **Rust**: mark the old item with
  `#[deprecated(since = "1.4.0", note = "use `new_name`; see the 2.0 migration guide")]`.
  The attribute is a compile-time warning, never an error, and the old item is
  covered by the freeze test.
- **Python**: emit `DeprecationWarning` once per process and name, with the
  replacement and the migration-guide path. The
  `xtalforge.deprecation_warning(name, replacement)` helper does this for
  extension code; pure-Python wrappers use `warnings.warn(..., DeprecationWarning)`.
- **WASM**: call `deprecation_warning(name, replacement)` once per page; it
  writes a `console.warn` with the same message. The function is exported so
  a wrapper package can use it without duplicating the policy.
- **CLI**: write one `warning[deprecated_name]: ...` line to standard error
  when a deprecated flag is used. Standard output and the JSON payload are
  unchanged, so scripts keep working.

## Removal rule

An alias may be removed only in a MAJOR release, after at least one MINOR
release with the warning, and every removal is listed in
`docs/api/naming.md` and the migration guide. The 2.0 checklist is Task 17.1.

## Checking

`scripts/check-rules.sh` keeps the 500-line and warning-suppression rules;
`scripts/check-parity.sh` keeps the surfaces aligned; the stability test in
`crates/xtalforge/tests` freezes the public names. A deprecated item is still
part of the frozen set until its removal release.
