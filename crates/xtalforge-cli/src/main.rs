//! Command-line interface for `xtalforge`.
//!
//! The command tree is introduced in Task 0.5.3 and grows with each phase.
//! This bootstrap binary only reports the workspace version so the crate,
//! the CI pipeline and the packaging checks have something to build.

fn main() {
    println!("xtalforge {}", env!("CARGO_PKG_VERSION"));
}
