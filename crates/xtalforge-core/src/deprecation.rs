//! Deprecation warnings shared by the surfaces (Task 11.4.2).

use std::collections::BTreeSet;
use std::sync::{Mutex, OnceLock};

static WARNED: OnceLock<Mutex<BTreeSet<String>>> = OnceLock::new();

/// Writes one deprecation warning to standard error, once per process.
///
/// The message names the replacement and the migration guide. It is used by
/// the CLI and by any host that has no native warning channel.
pub fn warn_deprecated(name: &str, replacement: &str) {
    let warned = WARNED.get_or_init(|| Mutex::new(BTreeSet::new()));
    let Ok(mut set) = warned.lock() else {
        return;
    };
    if !set.insert(name.to_owned()) {
        return;
    }
    eprintln!(
        "warning[deprecated_{name}]: use {replacement} instead; \
         see docs/guide/deprecation.md and the 2.0 migration guide"
    );
}
