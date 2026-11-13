//! `provenance` command (Task 16.2.2).

use crate::{Cli, CliError, input};

/// Prints the provenance record of a file or standard input.
pub(crate) fn run(cli: &Cli, path: &str) -> Result<(), CliError> {
    let text = input::read_source(path).map_err(CliError::core)?;
    let record = xtalforge::io::provenance(&text)
        .map_err(CliError::core)?
        .ok_or_else(|| CliError::input("the input carries no xtalforge provenance"))?;
    if cli.json {
        println!("{}", record.to_json().map_err(CliError::core)?);
    } else {
        println!(
            "xtalforge {} (schema {})",
            record.xtalforge_version, record.schema_version
        );
        for input in &record.inputs {
            println!(
                "input {}: {}",
                input.name.as_deref().unwrap_or("-"),
                input.sha256
            );
        }
        if !record.options.is_empty() {
            println!(
                "options: {}",
                serde_json::to_string(&record.options)
                    .map_err(|error| CliError::new("typing", error.to_string()))?
            );
        }
        if let Some(seed) = record.seed {
            println!("seed: {seed}");
        }
    }
    Ok(())
}
