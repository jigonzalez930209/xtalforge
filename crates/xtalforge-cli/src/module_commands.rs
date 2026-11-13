//! `modules` command: list the module registry (Task 14.5.1).

use xtalforge::modules::builtin_registry;

use crate::{Cli, CliError};

/// Prints the registered modules as text or JSON.
pub(crate) fn run(cli: &Cli) -> Result<(), CliError> {
    let registry = builtin_registry().map_err(CliError::core)?;
    let listing = registry.listing();
    if cli.json {
        let text = registry.listing_json().map_err(CliError::core)?;
        println!("{text}");
        return Ok(());
    }
    if listing.is_empty() {
        println!("no modules registered");
        return Ok(());
    }
    for module in listing {
        let state = if module.available {
            "available"
        } else {
            "unavailable"
        };
        println!("{} ({state}) — {}", module.name, module.description);
    }
    Ok(())
}
