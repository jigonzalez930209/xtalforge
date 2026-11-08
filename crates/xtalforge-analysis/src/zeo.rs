//! Zeo++ bridge: result records, output parsers and an optional runner
//! (Task 5.5).
//!
//! Zeo++ is an external program. Parsing its captured output needs no binary,
//! so CI tests the parsers against fabricated files; [`run_zeo`] probes for
//! `zeo++` on `PATH` and is compiled out of WASM targets. The parsers accept
//! both the labelled `Channels` form (`Pore volume`, `Accessible volume`,
//! `Surface area`, the three sphere/diameter lines) and the older
//! `POAV_*`/`AV_*` summary spelling.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Structure};

/// Options for one Zeo++ run.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ZeoOptions {
    /// Probe radius in Angstrom (Zeo++ `-probe_radius`).
    pub probe_radius: f64,
    /// Monte Carlo samples per atom (Zeo++ `-samples`).
    pub samples: u32,
    /// Whether to request the channel output (Zeo++ `-chan`).
    pub channels: bool,
}

impl Default for ZeoOptions {
    fn default() -> Self {
        Self {
            probe_radius: 1.2,
            samples: 1000,
            channels: false,
        }
    }
}

/// Zeo++ pore metrics, all per unit cell unless stated otherwise.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ZeoResult {
    /// Pore volume the probe centre can occupy, in Angstrom^3 per cell.
    pub pore_volume: f64,
    /// Pore volume reachable from outside the crystal, in Angstrom^3 per cell.
    pub accessible_volume: f64,
    /// Accessible surface area, in Angstrom^2 per cell.
    pub surface_area: f64,
    /// Largest included sphere diameter `D_i`, in Angstrom: the biggest
    /// cavity anywhere in the crystal, reachable or not.
    pub largest_included_sphere: f64,
    /// Largest free sphere diameter `D_f`, in Angstrom: the biggest sphere
    /// that can pass all the way through the framework.
    pub largest_free_sphere: f64,
    /// Pore limiting diameter, in Angstrom: the narrowest window the free
    /// sphere must pass through.
    pub pore_limiting_diameter: f64,
}

/// Parses Zeo++ `Channels` output into a [`ZeoResult`].
///
/// Labels are matched case-insensitively and tolerate units in parentheses.
/// The older `POAV_A^3`/`AV_A^3` spellings supply the pore and accessible
/// volumes, and `ASA_A^2` supplies the surface area.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] listing every required quantity whose line
/// is missing.
pub fn parse_zeo_channels(text: &str) -> Result<ZeoResult, Error> {
    let values = fields(text);
    let mut missing = Vec::new();
    let mut take = |key: &str, label: &'static str| -> f64 {
        match values.get(key) {
            Some(value) => *value,
            None => {
                missing.push(label);
                0.0
            }
        }
    };
    let pore_volume = take("pore_volume", "pore volume");
    let accessible_volume = take("accessible_volume", "accessible volume");
    let surface_area = take("surface_area", "surface area");
    let largest_included_sphere = take("largest_included_sphere", "largest included sphere");
    let largest_free_sphere = take("largest_free_sphere", "largest free sphere");
    let pore_limiting_diameter = take("pore_limiting_diameter", "pore limiting diameter");
    if !missing.is_empty() {
        return Err(Error::InvalidInput(format!(
            "Zeo++ channels output is missing: {}",
            missing.join(", ")
        )));
    }
    Ok(ZeoResult {
        pore_volume,
        accessible_volume,
        surface_area,
        largest_included_sphere,
        largest_free_sphere,
        pore_limiting_diameter,
    })
}

/// Reads the accessible volume from a Zeo++ `.vol` file, in Angstrom^3.
///
/// Both the `-vol` spelling (`AV_A^3`) and the occupyable `-volpo` spelling
/// (`POAV_A^3`) are accepted; when both appear, the first one read wins.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] when neither quantity is present.
pub fn parse_zeo_volume(text: &str) -> Result<f64, Error> {
    fields(text)
        .get("accessible_volume")
        .copied()
        .ok_or_else(|| {
            Error::InvalidInput(
                "Zeo++ .vol output has no accessible volume (AV_A^3 or POAV_A^3)".to_owned(),
            )
        })
}

/// Writes the CSSR input and the Zeo++ argument string for one structure.
///
/// The first string is a P1 CSSR file (all symmetry images expanded) for
/// `structure.cssr`. The second is a space-separated argument list for the
/// `network` program: `-probe_radius`, `-samples`, `-res result.res`,
/// `-sa result.sa`, `-volpo result.vol`, and `-chan result.chan` when
/// [`ZeoOptions::channels`] is set. The caller appends the CSSR path as the
/// positional argument.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an invalid option, an empty structure,
/// or a structure that cannot be expanded to P1.
pub fn write_zeo_inputs(
    structure: &Structure,
    options: &ZeoOptions,
) -> Result<(String, String), Error> {
    validate_options(options)?;
    let cssr = xtalforge_io::write_cssr(structure)?;
    let mut arguments = vec![
        format!("-probe_radius {}", options.probe_radius),
        format!("-samples {}", options.samples),
        "-res result.res".to_owned(),
        "-sa result.sa".to_owned(),
        "-volpo result.vol".to_owned(),
    ];
    if options.channels {
        arguments.push("-chan result.chan".to_owned());
    }
    Ok((cssr, arguments.join(" ")))
}

/// Runs Zeo++ on a structure and parses its channel output.
///
/// The binary is looked up as `zeo++` on `PATH` and runs in a fresh temporary
/// directory with a 120 second timeout. Its standard output is parsed with
/// [`parse_zeo_channels`]; if that fails, `result.chan` is tried.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] when the options are invalid, the structure
/// cannot be written, `zeo++` is not available, the process fails or times
/// out, or its output cannot be parsed.
#[cfg(not(target_arch = "wasm32"))]
pub fn run_zeo(structure: &Structure, options: &ZeoOptions) -> Result<ZeoResult, Error> {
    let (cssr, arguments) = write_zeo_inputs(structure, options)?;
    probe_zeo()?;
    let directory = TempDir::create()?;
    let input = directory.path.join("structure.cssr");
    std::fs::write(&input, cssr).map_err(|error| {
        Error::InvalidInput(format!("could not write {}: {error}", input.display()))
    })?;
    let mut command = std::process::Command::new("zeo++");
    command.current_dir(&directory.path);
    for argument in arguments.split_whitespace() {
        command.arg(argument);
    }
    command.arg("structure.cssr");
    let captured = run_with_timeout(&mut command, ZEO_TIMEOUT)?;
    if !captured.status.success() {
        let reason = captured.stderr.trim();
        return Err(Error::InvalidInput(format!(
            "zeo++ exited with {}: {reason}",
            captured.status
        )));
    }
    let parsed = parse_zeo_channels(&captured.stdout);
    if parsed.is_ok() {
        return parsed;
    }
    let channel_file = directory.path.join("result.chan");
    if let Ok(text) = std::fs::read_to_string(&channel_file)
        && let Ok(result) = parse_zeo_channels(&text)
    {
        return Ok(result);
    }
    parsed
}

fn validate_options(options: &ZeoOptions) -> Result<(), Error> {
    if !options.probe_radius.is_finite() || options.probe_radius < 0.0 {
        return Err(Error::InvalidInput(
            "the Zeo++ probe radius must be finite and non-negative".to_owned(),
        ));
    }
    if options.samples == 0 {
        return Err(Error::InvalidInput(
            "the Zeo++ sample count must be positive".to_owned(),
        ));
    }
    Ok(())
}

/// Canonical keys and the label spellings that produce them.
const ALIASES: &[(&str, &[&str])] = &[
    ("pore_volume", &["pore_volume", "poav_a_3"]),
    (
        "accessible_volume",
        &["accessible_volume", "av_a_3", "poav_a_3"],
    ),
    ("surface_area", &["surface_area", "asa_a_2"]),
    (
        "largest_included_sphere",
        &["largest_included_sphere", "included_sphere"],
    ),
    (
        "largest_free_sphere",
        &["largest_free_sphere", "free_sphere"],
    ),
    (
        "pore_limiting_diameter",
        &[
            "pore_limiting_diameter",
            "pore_limiting",
            "included_sphere_along_free_path",
        ],
    ),
];

fn fields(text: &str) -> BTreeMap<String, f64> {
    let mut values: BTreeMap<String, f64> = BTreeMap::new();
    for raw in text.lines() {
        let normalized = raw.replace('=', ":");
        let line = normalized.trim();
        if line.is_empty() {
            continue;
        }
        let mut offset = 0;
        while let Some(relative) = line.get(offset..).and_then(|rest| rest.find(':')) {
            let colon = offset + relative;
            let prefix = canonical(line.get(..colon).unwrap_or(""));
            let after = line.get(colon + 1..).unwrap_or("");
            let end = after.find(':').unwrap_or(after.len());
            if let Some(value) = first_number(after.get(..end).unwrap_or("")) {
                for key in alias_keys(&prefix) {
                    values.entry(key.to_owned()).or_insert(value);
                }
            }
            offset = colon + 1 + end;
            if end == 0 {
                offset = colon + 1;
            }
        }
    }
    values
}

fn alias_keys(prefix: &str) -> Vec<&'static str> {
    ALIASES
        .iter()
        .filter(|(_, aliases)| aliases.iter().any(|alias| prefix.ends_with(alias)))
        .map(|(key, _)| *key)
        .collect()
}

fn canonical(text: &str) -> String {
    let mut out = String::new();
    let mut depth = 0u32;
    for character in text.chars() {
        match character {
            '(' | '[' => depth += 1,
            ')' | ']' => depth = depth.saturating_sub(1),
            _ if depth > 0 => {}
            _ if character.is_ascii_alphanumeric() => {
                out.push(character.to_ascii_lowercase());
            }
            _ => {
                if !out.ends_with('_') {
                    out.push('_');
                }
            }
        }
    }
    out.trim_matches('_').to_owned()
}

fn first_number(text: &str) -> Option<f64> {
    text.split_whitespace().find_map(|token| {
        let cleaned = token.trim_matches(|character: char| {
            !(character.is_ascii_digit() || matches!(character, '.' | '-' | '+' | 'e' | 'E'))
        });
        cleaned.parse::<f64>().ok()
    })
}

#[cfg(not(target_arch = "wasm32"))]
const ZEO_TIMEOUT: std::time::Duration = std::time::Duration::from_secs(120);

#[cfg(not(target_arch = "wasm32"))]
fn probe_zeo() -> Result<(), Error> {
    use std::process::{Command, Stdio};

    match Command::new("zeo++")
        .arg("-h")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
    {
        Ok(_) => Ok(()),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => Err(Error::InvalidInput(
            "zeo++ is not available on PATH; install Zeo++ or parse captured output with parse_zeo_channels".to_owned(),
        )),
        Err(error) => Err(Error::InvalidInput(format!(
            "could not probe zeo++ on PATH: {error}"
        ))),
    }
}

#[cfg(not(target_arch = "wasm32"))]
struct Captured {
    status: std::process::ExitStatus,
    stdout: String,
    stderr: String,
}

#[cfg(not(target_arch = "wasm32"))]
fn run_with_timeout(
    command: &mut std::process::Command,
    timeout: std::time::Duration,
) -> Result<Captured, Error> {
    use std::process::Stdio;

    command.stdout(Stdio::piped()).stderr(Stdio::piped());
    let mut child = command
        .spawn()
        .map_err(|error| Error::InvalidInput(format!("could not run zeo++: {error}")))?;
    let stdout = child.stdout.take().map(read_stream);
    let stderr = child.stderr.take().map(read_stream);
    let started = std::time::Instant::now();
    loop {
        match child.try_wait() {
            Ok(Some(status)) => {
                return Ok(Captured {
                    status,
                    stdout: join_reader(stdout),
                    stderr: join_reader(stderr),
                });
            }
            Ok(None) if started.elapsed() >= timeout => {
                let _ = child.kill();
                let _ = child.wait();
                let _ = join_reader(stdout);
                let _ = join_reader(stderr);
                return Err(Error::InvalidInput(format!(
                    "zeo++ did not finish within {} seconds",
                    timeout.as_secs()
                )));
            }
            Ok(None) => std::thread::sleep(std::time::Duration::from_millis(20)),
            Err(error) => {
                let _ = child.kill();
                let _ = child.wait();
                let _ = join_reader(stdout);
                let _ = join_reader(stderr);
                return Err(Error::InvalidInput(format!(
                    "could not wait for zeo++: {error}"
                )));
            }
        }
    }
}

#[cfg(not(target_arch = "wasm32"))]
fn read_stream(mut stream: impl std::io::Read + Send + 'static) -> std::thread::JoinHandle<String> {
    std::thread::spawn(move || {
        let mut text = String::new();
        let _ = std::io::Read::read_to_string(&mut stream, &mut text);
        text
    })
}

#[cfg(not(target_arch = "wasm32"))]
fn join_reader(handle: Option<std::thread::JoinHandle<String>>) -> String {
    handle
        .map(|inner| inner.join().unwrap_or_default())
        .unwrap_or_default()
}

#[cfg(not(target_arch = "wasm32"))]
struct TempDir {
    path: std::path::PathBuf,
}

#[cfg(not(target_arch = "wasm32"))]
impl TempDir {
    fn create() -> Result<Self, Error> {
        use std::sync::atomic::{AtomicU64, Ordering};

        static COUNTER: AtomicU64 = AtomicU64::new(0);
        let base = std::env::temp_dir();
        for _ in 0..64 {
            let counter = COUNTER.fetch_add(1, Ordering::Relaxed);
            let nanos = std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map(|value| value.as_nanos())
                .unwrap_or(0);
            let path = base.join(format!(
                "xtalforge-zeo-{}-{nanos}-{counter}",
                std::process::id()
            ));
            match std::fs::create_dir(&path) {
                Ok(()) => return Ok(Self { path }),
                Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
                Err(error) => {
                    return Err(Error::InvalidInput(format!(
                        "could not create temporary directory {}: {error}",
                        path.display()
                    )));
                }
            }
        }
        Err(Error::InvalidInput(
            "could not create a unique temporary directory for zeo++".to_owned(),
        ))
    }
}

#[cfg(not(target_arch = "wasm32"))]
impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.path);
    }
}
