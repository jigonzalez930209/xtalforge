//! Remote structure fetchers (Task 13.3).

use clap::Subcommand;
use serde_json::json;

use crate::format::StructureFormat;
use crate::{Cli, CliError, output};

/// Provider selected by `xtalforge fetch`.
#[derive(Debug, Subcommand)]
pub(crate) enum FetchAction {
    /// Query an OPTIMADE provider.
    Optimade {
        /// Provider root, without `/structures`.
        #[arg(long, default_value = "https://optimade.materialsproject.org")]
        base: String,
        /// OPTIMADE filter string.
        #[arg(long, default_value = "")]
        filter: String,
        /// Maximum number of structures.
        #[arg(long, default_value_t = 1)]
        limit: usize,
        /// Zero-based page offset.
        #[arg(long, default_value_t = 0)]
        offset: usize,
        /// Output path; a multi-structure result is written as JSON.
        #[arg(long, value_name = "PATH")]
        out: Option<String>,
        /// Output format for a single structure.
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
    /// Fetch a COD entry by numeric identifier.
    Cod {
        /// COD identifier, for example `1000041`.
        id: String,
        /// Output path; defaults to standard output.
        #[arg(long, value_name = "PATH")]
        out: Option<String>,
        /// Output format; defaults to the output extension or CIF.
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
    /// Fetch a Materials Project structure by identifier.
    Mp {
        /// Materials Project identifier, for example `mp-149`.
        id: String,
        /// Materials Project API key.
        #[arg(long, env = "MP_API_KEY")]
        api_key: String,
        /// Output path; defaults to standard output.
        #[arg(long, value_name = "PATH")]
        out: Option<String>,
        /// Output format; defaults to the output extension or CIF.
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
}

pub(crate) fn run(cli: &Cli, action: &FetchAction) -> Result<(), CliError> {
    match action {
        FetchAction::Optimade {
            base,
            filter,
            limit,
            offset,
            out,
            format,
        } => {
            let limit = (*limit).max(1);
            let page = xtalforge::remote::optimade_query(base, filter, limit, *offset)
                .map_err(CliError::core)?;
            if page.structures.len() == 1
                && out.is_some()
                && let Some(item) = page.structures.first()
            {
                output::emit_structure(cli, &item.structure, out.as_deref(), *format)?;
                return Ok(());
            }
            let documents = page
                .structures
                .iter()
                .map(|item| {
                    json!({
                        "id": item.id,
                        "source": item.source,
                        "license": item.license,
                        "structure": item.structure.to_document(),
                    })
                })
                .collect::<Vec<_>>();
            let payload = json!({
                "provider": base,
                "count": documents.len(),
                "next": page.next,
                "structures": documents,
            });
            let text = serde_json::to_string_pretty(&payload)
                .map_err(|error| CliError::new("typing", error.to_string()))?;
            if let Some(path) = out.as_deref().filter(|path| *path != "-") {
                output::write_destination(path, &text)?;
                if !cli.quiet {
                    println!("wrote {} structures to {path}", documents.len());
                }
            } else {
                println!("{text}");
            }
            Ok(())
        }
        FetchAction::Cod { id, out, format } => {
            let item = xtalforge::remote::cod_fetch(id).map_err(CliError::core)?;
            output::emit_structure(cli, &item.structure, out.as_deref(), *format)
        }
        FetchAction::Mp {
            id,
            api_key,
            out,
            format,
        } => {
            let item =
                xtalforge::remote::materials_project_fetch(id, api_key).map_err(CliError::core)?;
            output::emit_structure(cli, &item.structure, out.as_deref(), *format)
        }
    }
}
