//! Agent commands: capabilities, render and scripted sessions (Phase 8).

use serde_json::{Value, json};

use crate::args::SkillAction;
use crate::format::PictureFormat;
use crate::{Cli, CliError, input};

fn agent_error(error: xtalforge::agent::Error) -> CliError {
    CliError::new(error.code(), error.to_string())
}

/// Prints what this installation can do.
pub(crate) fn capabilities(cli: &Cli) -> Result<(), CliError> {
    let report = xtalforge::agent::capability_report();
    if cli.json {
        let payload = serde_json::to_string(&report)
            .map_err(|error| CliError::new("json_error", error.to_string()))?;
        println!("{payload}");
        return Ok(());
    }
    println!("xtalforge {}", report.version);
    let features: Vec<String> = report
        .features
        .iter()
        .filter(|(_, enabled)| **enabled)
        .map(|(name, _)| name.clone())
        .collect();
    println!("features: {}", features.join(", "));
    for engine in &report.engines {
        let state = if engine.available {
            "available"
        } else {
            "missing"
        };
        println!("engine {:<8} {state}: {}", engine.name, engine.detail);
    }
    for bundle in &report.bundles {
        let state = if bundle.available {
            "available"
        } else {
            "missing"
        };
        println!("bundle {:<8} {state}: {}", bundle.name, bundle.detail);
    }
    println!("formats: {}", report.formats.len());
    println!("verbs: {}", report.verbs.len());
    Ok(())
}

/// Draws a structure to SVG or PNG.
pub(crate) fn render(
    cli: &Cli,
    path: &str,
    out: Option<&str>,
    format: PictureFormat,
    width: u32,
    height: u32,
) -> Result<(), CliError> {
    let structure = input::read_structure_path(path).map_err(CliError::core)?;
    let options = xtalforge::agent::RenderOptions {
        width,
        height,
        ..Default::default()
    };
    let (name, bytes) = match format {
        PictureFormat::Svg => (
            "svg",
            xtalforge::agent::render_svg(&structure, &options)
                .map_err(agent_error)?
                .into_bytes(),
        ),
        PictureFormat::Png => (
            "png",
            xtalforge::agent::render_png(&structure, &options).map_err(agent_error)?,
        ),
    };
    if let Some(out) = out {
        std::fs::write(out, &bytes)
            .map_err(|error| CliError::new("io", format!("{out}: {error}")))?;
    } else {
        use std::io::Write;
        std::io::stdout()
            .write_all(&bytes)
            .map_err(|error| CliError::new("io", error.to_string()))?;
    }
    if cli.json {
        println!(
            "{}",
            json!({
                "format": name,
                "width": width,
                "height": height,
                "bytes": bytes.len(),
                "path": out,
            })
        );
    }
    Ok(())
}

/// Serves the Model Context Protocol over stdio.
pub(crate) fn mcp() -> Result<(), CliError> {
    xtalforge::agent::serve_stdio().map_err(agent_error)
}

/// Prints, lists or installs the skill documents.
pub(crate) fn skill(cli: &Cli, action: &SkillAction) -> Result<(), CliError> {
    match action {
        SkillAction::Path => {
            let path = xtalforge::agent::default_directory();
            if cli.json {
                println!("{}", json!({"path": path.to_string_lossy()}));
            } else {
                println!("{}", path.display());
            }
            Ok(())
        }
        SkillAction::List => {
            let entries: Vec<Value> = xtalforge::agent::skill::DOCUMENTS
                .iter()
                .map(|document| {
                    json!({
                        "name": document.name,
                        "title": document.title,
                    })
                })
                .collect();
            if cli.json {
                println!("{}", json!(entries));
            } else {
                for document in xtalforge::agent::skill::DOCUMENTS {
                    println!("{:<32} {}", document.name, document.title);
                }
            }
            Ok(())
        }
        SkillAction::Show { name } => {
            let document = xtalforge::agent::skill::find(name).ok_or_else(|| {
                CliError::input(format!(
                    "unknown skill document {name:?}; try `xtalforge skill list`"
                ))
            })?;
            if cli.json {
                println!(
                    "{}",
                    json!({"name": document.name, "title": document.title, "body": document.body})
                );
            } else {
                print!("{}", document.body);
            }
            Ok(())
        }
        SkillAction::Install { dir } => {
            let directory = dir
                .as_ref()
                .map(std::path::PathBuf::from)
                .unwrap_or_else(xtalforge::agent::default_directory);
            let written = xtalforge::agent::skill::install(&directory).map_err(agent_error)?;
            if cli.json {
                println!(
                    "{}",
                    json!({
                        "directory": directory.to_string_lossy(),
                        "files": written.iter().map(|path| path.to_string_lossy()).collect::<Vec<_>>(),
                    })
                );
            } else {
                println!(
                    "installed {} files in {}",
                    written.len(),
                    directory.display()
                );
            }
            Ok(())
        }
    }
}

/// Runs a JSON list of verbs against one session.
pub(crate) fn session(cli: &Cli, path: &str) -> Result<(), CliError> {
    let text = input::read_source(path).map_err(CliError::core)?;
    let document: Value = serde_json::from_str(&text)
        .map_err(|error| CliError::new("json_error", format!("{path}: {error}")))?;
    let calls = document
        .get("verbs")
        .and_then(Value::as_array)
        .or_else(|| document.as_array())
        .ok_or_else(|| {
            CliError::input("an agent script is an array of {verb, arguments} objects")
        })?;
    let mut session = xtalforge::agent::Session::new();
    let mut answers = Vec::with_capacity(calls.len());
    for call in calls {
        let verb = call
            .get("verb")
            .and_then(Value::as_str)
            .ok_or_else(|| CliError::input("every script entry needs a string `verb`"))?;
        let arguments = call.get("arguments").cloned().unwrap_or_else(|| json!({}));
        let answer = session.run(verb, &arguments).map_err(agent_error)?;
        if !cli.json {
            println!(
                "{verb}: {}",
                serde_json::to_string(&answer)
                    .map_err(|error| CliError::new("json_error", error.to_string()))?
            );
        }
        answers.push(answer);
    }
    if cli.json {
        println!(
            "{}",
            serde_json::to_string(&answers)
                .map_err(|error| CliError::new("json_error", error.to_string()))?
        );
    }
    Ok(())
}
