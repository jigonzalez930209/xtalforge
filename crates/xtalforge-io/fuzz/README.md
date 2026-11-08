# CIF fuzz target

`cargo +nightly fuzz run cif` from `crates/xtalforge-io` feeds arbitrary UTF-8
into the CIF reader and asserts that no input panics. The target is optional:
the stable test suite uses the crafted malformed inputs in
`crates/xtalforge-io/tests/cif.rs`, and CI runs the fuzzer only when the
nightly toolchain is available.
