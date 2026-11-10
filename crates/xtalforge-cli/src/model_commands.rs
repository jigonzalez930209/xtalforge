//! `models` command: list, fetch, verify and remove ML weights (Task 8.4.2).

use serde_json::json;
use xtalforge::ff::models;

use crate::args::ModelAction;
use crate::{Cli, CliError};

fn error(error: xtalforge::ff::Error) -> CliError {
    CliError::new(error.code(), error.to_string())
}

pub(crate) fn run(cli: &Cli, action: &ModelAction) -> Result<(), CliError> {
    let cache = models::cache_directory();
    match action {
        ModelAction::List => {
            let statuses = models::list(&cache).map_err(error)?;
            if cli.json {
                println!(
                    "{}",
                    serde_json::to_string(&statuses)
                        .map_err(|error| CliError::new("json_error", error.to_string()))?
                );
            } else {
                for status in statuses {
                    let state = match (&status.record, status.verified) {
                        (None, _) => "missing".to_owned(),
                        (Some(_), Some(true)) => "cached, verified".to_owned(),
                        (Some(_), Some(false)) => "cached, CHECKSUM MISMATCH".to_owned(),
                        (Some(_), None) => "cached, artifact missing".to_owned(),
                    };
                    println!(
                        "{:<12} {:<10} {state}  {}",
                        status.spec.name, status.spec.license, status.path
                    );
                }
            }
            Ok(())
        }
        ModelAction::Fetch { name, url } => {
            let record = models::fetch(&cache, name, url).map_err(error)?;
            if cli.json {
                println!(
                    "{}",
                    json!({
                        "name": record.name,
                        "url": record.url,
                        "sha256": record.sha256,
                        "bytes": record.bytes,
                        "license": record.license,
                        "fetched": record.fetched,
                    })
                );
            } else {
                println!(
                    "fetched {} ({} bytes, {}): {}",
                    record.name, record.bytes, record.sha256, record.license
                );
            }
            Ok(())
        }
        ModelAction::Verify { name } => {
            let record = models::verify(&cache, name).map_err(error)?;
            if cli.json {
                println!(
                    "{}",
                    json!({"name": record.name, "sha256": record.sha256, "ok": true})
                );
            } else {
                println!("{} verified: {}", record.name, record.sha256);
            }
            Ok(())
        }
        ModelAction::Remove { name } => {
            let removed = models::remove(&cache, name).map_err(error)?;
            if cli.json {
                println!("{}", json!({"name": name, "removed": removed}));
            } else if removed {
                println!("removed {name}");
            } else {
                println!("{name} was not cached");
            }
            Ok(())
        }
    }
}
