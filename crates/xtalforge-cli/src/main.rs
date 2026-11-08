//! Command-line interface for `xtalforge`.
//!
//! The command tree grows with each phase. Phase 0 provides the global flags
//! and the hello-world commands `version` and `cell-volume`; the rest of the
//! surface is added as the domain crates land.

use clap::{Parser, Subcommand, ValueEnum};
use std::process::ExitCode;

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
    /// Limit policy (reserved for the limits feature).
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
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum Limits {
    Standard,
    Generous,
    WarnOnly,
}

#[derive(Debug)]
struct CliError {
    code: &'static str,
    message: String,
    exit_code: u8,
}

impl CliError {
    fn input(message: impl Into<String>) -> Self {
        Self {
            code: "invalid_input",
            message: message.into(),
            exit_code: 3,
        }
    }
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    match run(&cli) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            report_error(&cli, &error);
            ExitCode::from(error.exit_code)
        }
    }
}

fn run(cli: &Cli) -> Result<(), CliError> {
    if let Some(threads) = cli.threads
        && threads == 0
    {
        return Err(CliError::input("--threads must be a positive number"));
    }
    // Both flags are part of the stable CLI contract and are consumed by the
    // parallel and limits work in later phases; reading them here keeps the
    // skeleton warning-free.
    let _ = cli.limits;

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
            let element = xtalforge::core::element_info(input).map_err(|error| CliError {
                code: error.code(),
                message: error.message(),
                exit_code: 3,
            })?;
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
        Command::Cell { params } => {
            let [a, b, c, alpha, beta, gamma] = match params.as_slice() {
                [a, b, c, alpha, beta, gamma] => [*a, *b, *c, *alpha, *beta, *gamma],
                _ => {
                    return Err(CliError::input(
                        "cell needs six parameters: a b c alpha beta gamma",
                    ));
                }
            };
            let lattice =
                xtalforge::core::lattice_info(a, b, c, alpha, beta, gamma).map_err(|error| {
                    CliError {
                        code: error.code(),
                        message: error.message(),
                        exit_code: 3,
                    }
                })?;
            if cli.json {
                println!(
                    "{}",
                    serde_json::json!({
                        "parameters": lattice.parameters(),
                        "volume": lattice.volume(),
                        "matrix": lattice.matrix(),
                        "right_handed": lattice.is_right_handed(),
                    })
                );
            } else {
                let [a, b, c, alpha, beta, gamma] = lattice.parameters();
                println!(
                    "a={a:.6} b={b:.6} c={c:.6} alpha={alpha:.6} beta={beta:.6} gamma={gamma:.6} volume={:.6}",
                    lattice.volume()
                );
                for row in lattice.matrix() {
                    let [x, y, z] = row;
                    println!("{x:12.6} {y:12.6} {z:12.6}");
                }
            }
            Ok(())
        }
        Command::CellVolume {
            a,
            b,
            c,
            alpha,
            beta,
            gamma,
        } => {
            let volume = xtalforge::core::cell_volume(*a, *b, *c, *alpha, *beta, *gamma).map_err(
                |error| CliError {
                    code: error.code(),
                    message: error.message(),
                    exit_code: 3,
                },
            )?;
            if cli.json {
                println!("{}", serde_json::json!({ "volume": volume }));
            } else {
                println!("{volume}");
            }
            Ok(())
        }
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
