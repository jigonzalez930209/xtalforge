//! Command-line interface for `xtalforge`.
//!
//! Phase 1 provides structure input and conversion, supercells, slabs,
//! connectivity reports and geometry queries. Every command accepts the global
//! flags, emits one JSON value in `--json` mode and reports runtime failures
//! with exit code 3.

mod analyze;
mod basic;
mod build;
mod convert;
mod format;
mod input;
mod measure;
mod output;

use clap::{Parser, Subcommand, ValueEnum};
use std::process::ExitCode;

use format::StructureFormat;

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
    /// Worker threads (reserved for the parallel algorithms).
    #[arg(long, global = true, value_name = "N")]
    threads: Option<usize>,
    /// Limit policy used by the size checks.
    #[arg(long, global = true, value_enum, default_value_t = Limits::Standard)]
    limits: Limits,
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Print the xtalforge version.
    Version,
    /// Print reference data for an element.
    Element {
        /// Element symbol or atomic number.
        input: String,
    },
    /// Print a summary of a structure or structure JSON document.
    Info {
        /// Input path, or `-` for standard input.
        path: String,
    },
    /// Build a cell from its six parameters.
    Cell {
        /// a, b, c in Angstrom and alpha, beta, gamma in degrees.
        #[arg(long, num_args = 6, value_names = ["a", "b", "c", "alpha", "beta", "gamma"])]
        params: Vec<f64>,
    },
    /// Compute the volume of a cell from its six parameters.
    CellVolume {
        /// Cell length a in Angstrom.
        a: f64,
        /// Cell length b in Angstrom.
        b: f64,
        /// Cell length c in Angstrom.
        c: f64,
        /// Angle alpha in degrees.
        alpha: f64,
        /// Angle beta in degrees.
        beta: f64,
        /// Angle gamma in degrees.
        gamma: f64,
    },
    /// Convert a structure between CIF, POSCAR and XYZ.
    Convert {
        /// Input path, or `-` for standard input.
        path: String,
        /// Output format.
        #[arg(long, value_enum)]
        format: StructureFormat,
        /// Output path, or `-` for standard output.
        #[arg(long, value_name = "PATH")]
        out: Option<String>,
    },
    /// Build an na x nb x nc supercell.
    Supercell {
        /// Input path, or `-` for standard input.
        path: String,
        /// Repeats along a.
        na: u32,
        /// Repeats along b.
        nb: u32,
        /// Repeats along c.
        nc: u32,
        /// Output path; defaults to standard output.
        #[arg(long, value_name = "PATH")]
        out: Option<String>,
        /// Output format; defaults to the output extension or CIF.
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
    /// Cut a slab along the (h k l) plane.
    Slab {
        /// Input path, or `-` for standard input.
        path: String,
        /// First Miller index.
        h: i32,
        /// Second Miller index.
        k: i32,
        /// Third Miller index.
        l: i32,
        /// Number of layers.
        #[arg(long, value_name = "N")]
        layers: u32,
        /// Vacuum above the slab in Angstrom.
        #[arg(long, value_name = "A")]
        vacuum: f64,
        /// Output path; defaults to standard output.
        #[arg(long, value_name = "PATH")]
        out: Option<String>,
        /// Output format; defaults to the output extension or CIF.
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
    /// List the perceived bonds of the P1 cell.
    Bonds {
        /// Input path, or `-` for standard input.
        path: String,
    },
    /// List the chordless rings of the P1 cell.
    Rings {
        /// Input path, or `-` for standard input.
        path: String,
        /// Largest ring size in atoms.
        #[arg(long, value_name = "N", default_value_t = 6)]
        max_size: usize,
    },
    /// List the matched functional groups.
    Groups {
        /// Input path, or `-` for standard input.
        path: String,
    },
    /// List the periodic neighbors within a cutoff.
    Neighbors {
        /// Input path, or `-` for standard input.
        path: String,
        /// Cutoff radius in Angstrom.
        #[arg(long, value_name = "A")]
        cutoff: f64,
    },
    /// Minimum-image distance between two P1 atoms.
    Distance {
        /// Input path, or `-` for standard input.
        path: String,
        /// First atom index.
        i: usize,
        /// Second atom index.
        j: usize,
    },
    /// Angle i-j-k in degrees.
    Angle {
        /// Input path, or `-` for standard input.
        path: String,
        /// First atom index.
        i: usize,
        /// Pivot atom index.
        j: usize,
        /// Third atom index.
        k: usize,
    },
    /// Torsion i-j-k-l in degrees.
    Torsion {
        /// Input path, or `-` for standard input.
        path: String,
        /// First atom index.
        i: usize,
        /// Second atom index.
        j: usize,
        /// Third atom index.
        k: usize,
        /// Fourth atom index.
        l: usize,
    },
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum Limits {
    Standard,
    Generous,
    WarnOnly,
}

impl Limits {
    fn profile(self) -> xtalforge::core::Profile {
        match self {
            Self::Standard => xtalforge::core::Profile::Standard,
            Self::Generous => xtalforge::core::Profile::Generous,
            Self::WarnOnly => xtalforge::core::Profile::WarnOnly,
        }
    }
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
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    match run(&cli) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            report_error(&cli, &error);
            ExitCode::from(ERROR_EXIT_CODE)
        }
    }
}

fn run(cli: &Cli) -> Result<(), CliError> {
    if let Some(threads) = cli.threads
        && threads == 0
    {
        return Err(CliError::input("--threads must be a positive number"));
    }
    match &cli.command {
        Command::Version => {
            let version = xtalforge::core::version();
            if cli.json {
                println!("{}", serde_json::json!({ "version": version }));
            } else if cli.quiet {
                println!("{version}");
            } else {
                println!("xtalforge {version}");
            }
            Ok(())
        }
        Command::Element { input } => {
            let element = xtalforge::core::element_info(input).map_err(CliError::core)?;
            if cli.json {
                println!(
                    "{}",
                    serde_json::json!({
                        "z": element.z,
                        "symbol": element.symbol,
                        "name": element.name,
                        "mass": element.mass,
                        "covalent_radius": element.covalent_radius,
                        "vdw_radius": element.vdw_radius,
                        "color": element.color,
                        "is_metal": element.is_metal,
                        "valence": element.valence,
                        "electronegativity": element.electronegativity,
                        "pyykko_single": element.pyykko_single,
                        "pyykko_double": element.pyykko_double,
                        "pyykko_aromatic": element.pyykko_aromatic,
                    })
                );
            } else if cli.quiet {
                println!("{}", element.symbol);
            } else {
                let [red, green, blue] = element.color;
                println!(
                    "{} {} (Z={}) mass={:.4} covalent={:.3} vdw={:.3} color=#{red:02X}{green:02X}{blue:02X}",
                    element.symbol,
                    element.name,
                    element.z,
                    element.mass,
                    element.covalent_radius,
                    element.vdw_radius
                );
            }
            Ok(())
        }
        Command::Info { path } => basic::info(cli, path),
        Command::Cell { params } => basic::cell(cli, params),
        Command::CellVolume {
            a,
            b,
            c,
            alpha,
            beta,
            gamma,
        } => basic::cell_volume(cli, *a, *b, *c, *alpha, *beta, *gamma),
        Command::Convert { path, format, out } => convert::run(cli, path, *format, out.as_deref()),
        Command::Supercell {
            path,
            na,
            nb,
            nc,
            out,
            format,
        } => build::supercell(cli, path, [*na, *nb, *nc], out.as_deref(), *format),
        Command::Slab {
            path,
            h,
            k,
            l,
            layers,
            vacuum,
            out,
            format,
        } => build::slab(
            cli,
            path,
            [*h, *k, *l],
            *layers,
            *vacuum,
            out.as_deref(),
            *format,
        ),
        Command::Bonds { path } => analyze::bonds(cli, path),
        Command::Rings { path, max_size } => analyze::rings(cli, path, *max_size),
        Command::Groups { path } => analyze::groups(cli, path),
        Command::Neighbors { path, cutoff } => analyze::neighbors(cli, path, *cutoff),
        Command::Distance { path, i, j } => measure::distance(cli, path, *i, *j),
        Command::Angle { path, i, j, k } => measure::angle(cli, path, *i, *j, *k),
        Command::Torsion { path, i, j, k, l } => measure::torsion(cli, path, *i, *j, *k, *l),
    }
}

fn report_error(cli: &Cli, error: &CliError) {
    if cli.json {
        let payload = serde_json::json!({
            "error": {
                "code": error.code,
                "message": error.message,
            }
        });
        eprintln!("{payload}");
    } else {
        eprintln!("error: {}", error.message);
    }
}
