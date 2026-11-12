//! `xtalforge doctor`: installation, structure diagnostics and next commands.

use serde_json::json;

use crate::{Cli, CliError, input};

/// Runs the doctor report.
///
/// # Errors
///
/// Returns [`CliError`] when the optional structure cannot be read or diagnosed.
pub(crate) fn run(cli: &Cli, path: Option<&str>) -> Result<(), CliError> {
    let capabilities = xtalforge::agent::capability_report();
    let diagnostics = match path {
        Some(path) => {
            let structure = input::read_structure_path(path).map_err(CliError::core)?;
            Some(
                xtalforge::agent::diagnose(&structure)
                    .map_err(|error| CliError::new(error.code(), error.to_string()))?,
            )
        }
        None => None,
    };
    let suggestions = suggestions(path, diagnostics.as_ref());
    if cli.json {
        println!(
            "{}",
            json!({
                "version": capabilities.version,
                "features": capabilities.features,
                "engines": capabilities.engines,
                "bundles": capabilities.bundles,
                "formats": capabilities.formats.len(),
                "verbs": capabilities.verbs.len(),
                "diagnostics": diagnostics,
                "suggestions": suggestions,
            })
        );
        return Ok(());
    }
    println!("xtalforge {}", capabilities.version);
    let enabled: Vec<&str> = capabilities
        .features
        .iter()
        .filter(|(_, value)| **value)
        .map(|(name, _)| name.as_str())
        .collect();
    println!("features: {}", enabled.join(", "));
    for engine in &capabilities.engines {
        println!(
            "engine {:<8} {}: {}",
            engine.name,
            if engine.available { "ok" } else { "missing" },
            engine.detail
        );
    }
    if let Some(report) = &diagnostics {
        println!(
            "structure: {} diagnostics ({} errors)",
            report.diagnostics.len(),
            report.counts.get("error").copied().unwrap_or(0)
        );
        for diagnostic in report.diagnostics.iter().take(10) {
            println!(
                "  [{}] {}: {}",
                diagnostic.severity.as_str(),
                diagnostic.code,
                diagnostic.message
            );
        }
    } else {
        println!("structure: not supplied");
    }
    println!("next:");
    for suggestion in suggestions {
        println!("  {suggestion}");
    }
    Ok(())
}

fn suggestions(
    path: Option<&str>,
    diagnostics: Option<&xtalforge::agent::DiagnosticReport>,
) -> Vec<String> {
    let mut suggestions = Vec::new();
    if let Some(path) = path {
        suggestions.push(format!("xtalforge info {path}"));
        suggestions.push(format!("xtalforge symmetry detect {path}"));
        suggestions.push(format!("xtalforge pxrd {path}"));
        suggestions.push(format!("xtalforge ff energy {path}"));
        if diagnostics.is_some_and(|report| {
            report
                .diagnostics
                .iter()
                .any(|diagnostic| diagnostic.code.contains("symmetry"))
        }) {
            suggestions.push(format!(
                "xtalforge convert {path} --format cif --out fixed.cif  # check --spacegroup"
            ));
        }
    } else {
        suggestions.push("xtalforge capabilities".to_owned());
        suggestions.push("xtalforge doctor structure.cif".to_owned());
        suggestions.push("xtalforge --help".to_owned());
    }
    suggestions
}
