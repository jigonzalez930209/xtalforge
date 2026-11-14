//! Runs every recipe of `docs/guide/recipes.md` against the CLI.
//!
//! Task 17.1.5: the recipe table is executable. The test parses the CLI
//! column, substitutes the documented fixture names with the repository
//! corpus in a temporary directory, and runs each command with `--json`.
//! Commands must exit with the documented success code (0), emit only JSON on
//! standard output and never print a panic report.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::thread;
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

const RECIPES: &str = include_str!("../../../docs/guide/recipes.md");
const NACL: &str = include_str!("../../../tests/corpus/simple/NaCl.cif");
const SI: &str = include_str!("../../../tests/corpus/simple/Si.cif");
const CU: &str = include_str!("../../../tests/corpus/simple/Cu.cif");

/// Wall-clock budget for one recipe. The slowest documented command
/// (`group`) takes about three seconds in a debug build.
const TIMEOUT: Duration = Duration::from_secs(60);

/// Recipes that need an external engine or the network.
///
/// Every documented recipe runs natively: UFF, EEQ, MOF checking, native
/// porosity and Widom insertion are all pure Rust, and no recipe fetches a
/// structure or invokes a compiled program. The list is kept as the single
/// allowed exception (Task 17.1.5): an entry here is exercised with
/// `--help` instead, so a broken subcommand still fails the suite.
const EXTERNAL: &[&str] = &[];

/// Captured result of one CLI invocation.
struct Captured {
    status: Option<std::process::ExitStatus>,
    stdout: Vec<u8>,
    stderr: Vec<u8>,
    timed_out: bool,
}

/// Temporary directory removed when the test ends.
struct Scratch(PathBuf);

impl Scratch {
    fn new() -> Self {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("the clock must be after the epoch")
            .as_nanos();
        let path =
            std::env::temp_dir().join(format!("xtalforge-recipes-{}-{nanos}", std::process::id()));
        fs::create_dir_all(&path).expect("the recipe directory must be creatable");
        Self(path)
    }

    fn path(&self) -> &Path {
        &self.0
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

/// Returns the `(number, command)` pairs of the recipe table.
fn table_recipes() -> Vec<(u32, String)> {
    let mut recipes = Vec::new();
    for line in RECIPES.lines() {
        let cells: Vec<&str> = line
            .trim()
            .trim_matches('|')
            .split('|')
            .map(str::trim)
            .collect();
        if cells.len() != 6 {
            continue;
        }
        let Ok(number) = cells[0].parse::<u32>() else {
            continue;
        };
        let cli = cells[2].trim_matches('`');
        assert!(
            cli.starts_with("xtalforge "),
            "recipe {number} has no CLI column: {cli}"
        );
        recipes.push((number, cli.to_owned()));
    }
    recipes
}

/// Splits a CLI column into arguments, honouring single and double quotes.
fn tokenize(command: &str) -> Vec<String> {
    let mut tokens = Vec::new();
    let mut current = String::new();
    let mut quote = None;
    for character in command.chars() {
        match (quote, character) {
            (Some(open), value) if value == open => quote = None,
            (Some(_), value) => current.push(value),
            (None, '\'' | '"') => quote = Some(character),
            (None, value) if value.is_whitespace() => {
                if !current.is_empty() {
                    tokens.push(std::mem::take(&mut current));
                }
            }
            (None, value) => current.push(value),
        }
    }
    assert!(quote.is_none(), "unclosed quote in recipe: {command}");
    if !current.is_empty() {
        tokens.push(current);
    }
    tokens
}

/// Applies the two documented argument adaptations.
///
/// The recipe page itself is not modified by this task, so every deviation is
/// recorded here:
///
/// - Recipe 6: the CLI requires `--vacuum` and the page omits it; 15 Angstrom
///   is the value the Python column of the same recipe uses.
/// - Recipe 17: the documented fixtures are NaCl (`a.cif`) and Si (`b.cif`),
///   which are not isomorphous, so interpolation cannot map them. The test
///   builds two seeded frames of `a.cif` and interpolates those instead.
fn adapt(number: u32, args: &mut Vec<String>, scratch: &Path) {
    match number {
        6 => {
            args.push("--vacuum".to_owned());
            args.push("15".to_owned());
        }
        17 => {
            let first = scratch.join("interpolate_a.cif");
            let second = scratch.join("interpolate_b.cif");
            for argument in args.iter_mut() {
                match argument.as_str() {
                    "a.cif" => *argument = first.display().to_string(),
                    "b.cif" => *argument = second.display().to_string(),
                    _ => {}
                }
            }
        }
        _ => {}
    }
}

/// Replaces the documented fixture names with absolute temporary paths.
fn resolve_fixtures(args: &mut [String], scratch: &Path) {
    for argument in args {
        let fixture = match argument.as_str() {
            "a.cif" => Some("a.cif"),
            "b.cif" => Some("b.cif"),
            "c.cif" => Some("c.cif"),
            _ => None,
        };
        if let Some(name) = fixture {
            *argument = scratch.join(name).display().to_string();
        }
    }
}

/// Runs the binary with a wall-clock limit, draining both pipes concurrently.
fn run_cli(args: &[String], scratch: &Path) -> Captured {
    let mut child = Command::new(env!("CARGO_BIN_EXE_xtalforge"))
        .args(args)
        .current_dir(scratch)
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("the CLI binary must run");
    let mut stdout = child.stdout.take().expect("stdout must be piped");
    let mut stderr = child.stderr.take().expect("stderr must be piped");
    let stdout_reader = thread::spawn(move || {
        let mut buffer = Vec::new();
        let _ = stdout.read_to_end(&mut buffer);
        buffer
    });
    let stderr_reader = thread::spawn(move || {
        let mut buffer = Vec::new();
        let _ = stderr.read_to_end(&mut buffer);
        buffer
    });
    let start = Instant::now();
    let mut timed_out = false;
    loop {
        match child.try_wait().expect("the child status must be readable") {
            Some(_) => break,
            None if start.elapsed() > TIMEOUT => {
                let _ = child.kill();
                timed_out = true;
                break;
            }
            None => thread::sleep(Duration::from_millis(10)),
        }
    }
    let status = child.wait().expect("the child must be reaped");
    Captured {
        status: Some(status),
        stdout: stdout_reader.join().expect("the stdout reader must finish"),
        stderr: stderr_reader.join().expect("the stderr reader must finish"),
        timed_out,
    }
}

/// Asserts that a command succeeded, printed JSON and did not panic.
fn assert_success(number: u32, command: &str, captured: &Captured) {
    let stderr = String::from_utf8_lossy(&captured.stderr);
    let stdout = String::from_utf8_lossy(&captured.stdout);
    assert!(
        !captured.timed_out,
        "recipe {number} timed out after {TIMEOUT:?}: {command}"
    );
    assert_eq!(
        captured
            .status
            .as_ref()
            .and_then(std::process::ExitStatus::code),
        Some(0),
        "recipe {number} failed: {command}\nstdout:\n{stdout}\nstderr:\n{stderr}"
    );
    assert!(
        !stderr.contains("panicked at") && !stdout.contains("panicked at"),
        "recipe {number} panicked: {command}\n{stderr}"
    );
    for line in stdout.lines().filter(|line| !line.trim().is_empty()) {
        serde_json::from_str::<serde_json::Value>(line).unwrap_or_else(|error| {
            panic!("recipe {number} printed non-JSON output: {error}\nline: {line}")
        });
    }
}

#[test]
fn the_recipe_table_is_complete() {
    let recipes = table_recipes();
    assert!(
        recipes.len() >= 30,
        "the recipe page must keep at least thirty rows, found {}",
        recipes.len()
    );
    let numbers: Vec<u32> = recipes.iter().map(|(number, _)| *number).collect();
    assert_eq!(numbers, (1..=recipes.len() as u32).collect::<Vec<u32>>());
}

#[test]
fn external_recipes_answer_help() {
    for entry in EXTERNAL {
        let mut arguments = tokenize(entry);
        arguments.push("--help".to_owned());
        let captured = run_cli(&arguments, Path::new("."));
        assert!(
            captured
                .status
                .as_ref()
                .is_some_and(std::process::ExitStatus::success),
            "the external recipe `{entry}` has no --help: {}",
            String::from_utf8_lossy(&captured.stderr)
        );
    }
}

#[test]
fn every_recipe_runs() {
    let recipes = table_recipes();
    let scratch = Scratch::new();
    // Documented fixture names: `a.cif` is the NaCl corpus structure and
    // `b.cif` is Si. Recipe 14 also uses `c.cif`, which has no documented
    // source; Cu keeps it a third distinct cubic structure.
    fs::write(scratch.path().join("a.cif"), NACL).expect("NaCl fixture");
    fs::write(scratch.path().join("b.cif"), SI).expect("Si fixture");
    fs::write(scratch.path().join("c.cif"), CU).expect("Cu fixture");

    // Recipe 17 needs two frames of one structure; build shifted frames with
    // the documented seeded rattle command.
    for (seed, name) in [("7", "interpolate_a.cif"), ("8", "interpolate_b.cif")] {
        let setup = run_cli(
            &[
                "features".to_owned(),
                "rattle".to_owned(),
                "a.cif".to_owned(),
                "--seed".to_owned(),
                seed.to_owned(),
                "--out".to_owned(),
                name.to_owned(),
            ],
            scratch.path(),
        );
        assert!(
            setup
                .status
                .as_ref()
                .is_some_and(std::process::ExitStatus::success),
            "the interpolation fixture {name} must be built: {}",
            String::from_utf8_lossy(&setup.stderr)
        );
    }

    for (number, command) in &recipes {
        if EXTERNAL.iter().any(|entry| command.contains(entry)) {
            continue;
        }
        let mut args = tokenize(command);
        assert_eq!(args.first().map(String::as_str), Some("xtalforge"));
        args.remove(0);
        adapt(*number, &mut args, scratch.path());
        resolve_fixtures(&mut args, scratch.path());
        let mut full = vec!["--json".to_owned()];
        full.extend(args);
        let captured = run_cli(&full, scratch.path());
        assert_success(*number, command, &captured);
    }
}
