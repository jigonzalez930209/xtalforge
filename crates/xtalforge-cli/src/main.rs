//! Command-line interface for `xtalforge`.
//!
//! Phase 1 provides structure input and conversion, supercells, slabs,
//! connectivity reports and geometry queries. Every command accepts the global
//! flags, emits one JSON value in `--json` mode and reports runtime failures
//! with exit code 3.

mod agent_commands;
mod analysis_commands;
mod analyze;
mod args;
mod basic;
mod build;
mod build_structures;
mod commands_extra;
mod convert;
mod dispatch;
mod engine_commands;
mod ff_commands;
mod format;
mod input;
mod limits;
mod measure;
mod mesh;
mod mof_commands;
mod output;
mod powder_commands;
mod pxrd;
mod spacegroup;
mod symmetry;
mod symmetry_analysis;

use clap::Parser;
use std::process::ExitCode;

use args::Command;
use limits::Limits;

/// Exit code used for every runtime failure.
const ERROR_EXIT_CODE: u8 = 3;

#[derive(Debug, Parser)]
#[command(
    name = "xtalforge",
    version,
    about = "Build, manipulate, analyse and export crystal structures",
    disable_help_subcommand = true
)]
struct Cli {
    /// Emit machine-readable JSON instead of human-readable text.
    #[arg(long, global = true)]
    json: bool,
    /// Suppress headers and non-essential output.
    #[arg(long, global = true)]
    quiet: bool,
    /// Turn loss warnings into errors.
    #[arg(long, global = true)]
    strict: bool,
    /// Worker threads (reserved for the parallel algorithms).
    #[arg(long, global = true, value_name = "N")]
    threads: Option<usize>,
    /// Limit policy used by the size checks.
    #[arg(long, global = true, value_enum, default_value_t = Limits::Standard)]
    limits: Limits,
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug)]
pub(crate) struct CliError {
    pub(crate) code: &'static str,
    pub(crate) message: String,
}

impl CliError {
    pub(crate) fn new(code: &'static str, message: impl Into<String>) -> Self {
        Self {
            code,
            message: message.into(),
        }
    }

    pub(crate) fn input(message: impl Into<String>) -> Self {
        Self::new("invalid_input", message)
    }

    pub(crate) fn core(error: xtalforge::core::Error) -> Self {
        Self::new(error.code(), error.message())
    }

    pub(crate) fn ff(error: xtalforge::ff::Error) -> Self {
        Self::new(error.code(), error.message())
    }

    pub(crate) fn analysis(error: xtalforge::analysis::Error) -> Self {
        Self::new(error.code(), error.message())
    }

    pub(crate) fn mof(error: xtalforge::mof::Error) -> Self {
        Self::new(error.code(), error.message())
    }
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    match dispatch::run(&cli) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            dispatch::report_error(&cli, &error);
            ExitCode::from(ERROR_EXIT_CODE)
        }
    }
}
