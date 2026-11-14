//! Integration tests for workflow compute steps, checkpoints and schedulers
//! (Task 18.4).

#![cfg(unix)]
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::fs;
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

#[cfg(feature = "external")]
use xtalforge_engines::Error;
use xtalforge_engines::input_set::InputFile;
use xtalforge_engines::job::{
    ComputeStep, EngineCommand, InputFiles, JobCheckpoint, JobOptions, JobRequest, JobState,
    ParserRegistry, RunOptions, Scheduler, SchedulerConfig, checkpoint_path, corrective_action,
    detect_failure, failure_catalogue, poll, probe_scheduler, process_execution_available, run_job,
    scheduler_script, submit,
};
#[cfg(feature = "external")]
use xtalforge_engines::job::{JOB_CHECKPOINT_VERSION, collect};
#[cfg(feature = "external")]
use xtalforge_engines::result::{CalculationResult, Convergence};

static COUNTER: AtomicU64 = AtomicU64::new(0);

fn temp_dir(tag: &str) -> PathBuf {
    let id = COUNTER.fetch_add(1, Ordering::Relaxed);
    let dir =
        std::env::temp_dir().join(format!("xtalforge-jobs-{}-{tag}-{id}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    dir
}

fn script(dir: &Path, name: &str, body: &str) -> PathBuf {
    let path = dir.join(name);
    fs::write(&path, format!("#!/bin/sh\nset -eu\n{body}\n")).unwrap();
    let mut permissions = fs::metadata(&path).unwrap().permissions();
    permissions.set_mode(0o755);
    fs::set_permissions(&path, permissions).unwrap();
    path
}

fn options(scheduler: Scheduler) -> JobOptions {
    let mut options = JobOptions::new(scheduler, "mpirun vasp_std");
    options.queue = "gpu".to_owned();
    options.nodes = 2;
    options.walltime = "02:30:00".to_owned();
    options.job_name = "mof".to_owned();
    options.output = "run.out".to_owned();
    options.error = "run.err".to_owned();
    options
}

fn fake_config(scheduler: Scheduler, submit: &Path, status: &Path) -> SchedulerConfig {
    SchedulerConfig {
        options: options(scheduler),
        submit_command: Some(submit.to_string_lossy().into_owned()),
        status_command: Some(status.to_string_lossy().into_owned()),
    }
}

fn compute_request(engine: &Path, parser: &str) -> JobRequest {
    JobRequest {
        name: "calc".to_owned(),
        inputs: InputFiles::one("input.txt", "payload"),
        command: EngineCommand::new(engine.to_string_lossy().into_owned(), Vec::new()),
        parser: parser.to_owned(),
        output: "result.txt".to_owned(),
    }
}

#[cfg(feature = "external")]
fn fake_parser(text: &str) -> Result<CalculationResult, Error> {
    if !text.contains("payload") {
        return Err(Error::Parse("missing payload".to_owned()));
    }
    Ok(CalculationResult {
        energies: vec![1.5],
        convergence: Convergence {
            scf: true,
            ionic: true,
            ..Convergence::default()
        },
        ..CalculationResult::default()
    })
}

#[test]
fn slurm_script_has_every_directive() {
    let script = scheduler_script(&options(Scheduler::Slurm)).unwrap();
    let expected = r#"#!/bin/bash
#SBATCH --job-name=mof
#SBATCH --partition=gpu
#SBATCH --nodes=2
#SBATCH --time=02:30:00
#SBATCH --output=run.out
#SBATCH --error=run.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
mpirun vasp_std
"#;
    assert_eq!(script, expected);
}

#[test]
fn pbs_script_has_every_directive() {
    let script = scheduler_script(&options(Scheduler::Pbs)).unwrap();
    let expected = r#"#!/bin/bash
#PBS -N mof
#PBS -q gpu
#PBS -l nodes=2
#PBS -l walltime=02:30:00
#PBS -o run.out
#PBS -e run.err
set -euo pipefail
cd "${PBS_O_WORKDIR:-$PWD}"
mpirun vasp_std
"#;
    assert_eq!(script, expected);
}

#[test]
fn invalid_scheduler_options_are_refused() {
    let unknown = scheduler_script(&options(Scheduler::Unknown)).unwrap_err();
    assert_eq!(unknown.code(), "engine_job");

    let mut empty = options(Scheduler::Slurm);
    empty.command = "   ".to_owned();
    assert!(scheduler_script(&empty).is_err());

    let mut injected = options(Scheduler::Slurm);
    injected.walltime = "01:00:00\n#SBATCH --nodes=99".to_owned();
    assert!(scheduler_script(&injected).is_err());

    let mut escaping = options(Scheduler::Slurm);
    escaping.output = "../escape.out".to_owned();
    assert!(scheduler_script(&escaping).is_err());

    let mut no_nodes = options(Scheduler::Slurm);
    no_nodes.nodes = 0;
    assert!(scheduler_script(&no_nodes).is_err());

    let mut path_job = options(Scheduler::Slurm);
    path_job.job_name = "../escape".to_owned();
    assert!(scheduler_script(&path_job).is_err());
}

#[test]
fn dry_run_writes_only_the_script() {
    let dir = temp_dir("dry-submit");
    let marker = dir.join("submitted");
    let submitter = script(&dir, "fake-sbatch", &format!("touch {}", marker.display()));
    let status = script(&dir, "fake-squeue", "echo PENDING");
    let mut config = fake_config(Scheduler::Slurm, &submitter, &status);
    config.options.dry_run = true;
    let submitted = submit(&config, &dir).unwrap();
    assert!(submitted.dry_run);
    assert_eq!(submitted.job.id, "dry-run");
    assert_eq!(submitted.job.state, JobState::Prepared);
    assert!(submitted.script.ends_with("mof.slurm"));
    assert_eq!(
        fs::read_to_string(&submitted.script).unwrap(),
        scheduler_script(&config.options).unwrap()
    );
    assert!(!marker.exists());
}

#[test]
fn probes_report_their_commands() {
    for scheduler in [Scheduler::Slurm, Scheduler::Pbs] {
        let probe = probe_scheduler(scheduler);
        assert_eq!(probe.scheduler, scheduler);
        assert!(!probe.submit_command.is_empty());
        assert!(!probe.status_command.is_empty());
        if probe.available {
            assert_eq!(probe.reason, "scheduler commands are available");
        } else {
            assert!(!probe.reason.is_empty());
        }
    }
    let unknown = probe_scheduler(Scheduler::Unknown);
    assert!(!unknown.available);
    assert!(unknown.reason.contains("unknown scheduler"));
}

#[test]
fn submit_without_a_scheduler_command_is_unavailable() {
    let dir = temp_dir("missing-submit");
    let mut config = SchedulerConfig::new(options(Scheduler::Slurm));
    config.submit_command = Some("xtalforge-no-such-sbatch".to_owned());
    let error = submit(&config, &dir).unwrap_err();
    assert_eq!(error.code(), "engine_unavailable");
    assert!(dir.join("mof.slurm").exists());
}

#[test]
fn every_documented_failure_pattern_is_detected() {
    let cases = [
        (
            "scf_not_converged",
            "Error: SCF not converged after 60 iterations",
        ),
        ("zbrent", "ZBRENT: fatal error in bracketing"),
        (
            "walltime_exceeded",
            "slurmstepd: JOB 42 CANCELLED AT 2026-01-01 DUE TO TIME LIMIT",
        ),
        ("out_of_memory", "slurmstepd: Detected 1 oom-kill event(s)"),
        (
            "missing_pseudopotential",
            "ERROR: the following species are missing in the POTCAR: Hf",
        ),
    ];
    for (code, text) in cases {
        let pattern = detect_failure(text).unwrap();
        assert_eq!(pattern.code, code);
        assert!(!pattern.pattern.is_empty());
        assert!(!pattern.action.is_empty());
        assert_eq!(corrective_action(code), Some(pattern.action.as_str()));
    }
    assert!(detect_failure("normal completion, 1 SCF step").is_none());
    let catalogue = failure_catalogue();
    for expected in [
        "scf_not_converged",
        "zbrent",
        "walltime_exceeded",
        "out_of_memory",
        "missing_pseudopotential",
    ] {
        assert!(catalogue.iter().any(|pattern| pattern.code == expected));
    }
}

#[test]
fn poll_rejects_dry_runs_and_bad_ids() {
    let mut dry_options = options(Scheduler::Slurm);
    dry_options.dry_run = true;
    let config = SchedulerConfig::new(dry_options);
    assert_eq!(poll(&config, "42").unwrap_err().code(), "engine_job");
    let config = SchedulerConfig::new(options(Scheduler::Slurm));
    assert_eq!(poll(&config, "  ").unwrap_err().code(), "engine_job");
}

#[test]
fn input_files_reject_unsafe_names() {
    let dir = temp_dir("input-paths");
    assert_eq!(
        InputFiles::one("../escape.txt", "x")
            .write_to(&dir)
            .unwrap_err()
            .code(),
        "engine_job"
    );
    let duplicate = InputFiles::new(vec![
        InputFile {
            name: "a.txt".to_owned(),
            text: "1".to_owned(),
        },
        InputFile {
            name: "a.txt".to_owned(),
            text: "2".to_owned(),
        },
    ]);
    assert!(duplicate.validate().is_err());
}

#[test]
fn compute_dry_run_writes_inputs_and_never_runs() {
    let dir = temp_dir("compute-dry");
    let marker = dir.join("ran");
    let engine = script(&dir, "engine", &format!("touch {}", marker.display()));
    let request = compute_request(&engine, "fake");
    let mut options = RunOptions::new(&dir);
    options.dry_run = true;
    let outcome = run_job(&request, &ParserRegistry::new(), &options).unwrap();
    assert_eq!(outcome.state, JobState::Prepared);
    assert!(outcome.dry_run);
    assert_eq!(outcome.attempts, 0);
    assert_eq!(
        fs::read_to_string(dir.join("calc").join("input.txt")).unwrap(),
        "payload"
    );
    assert!(!marker.exists());
    assert!(!checkpoint_path(&dir, "calc").unwrap().exists());
}

#[cfg(not(feature = "external"))]
#[test]
fn compute_without_external_writes_inputs_then_reports_unavailable() {
    let dir = temp_dir("compute-no-external");
    let request = compute_request(Path::new("true"), "fake");
    let options = RunOptions::new(&dir);
    let error = run_job(&request, &ParserRegistry::new(), &options).unwrap_err();
    assert_eq!(error.code(), "engine_unavailable");
    assert!(dir.join("calc").join("input.txt").exists());
    let checkpoint = JobCheckpoint::read(checkpoint_path(&dir, "calc").unwrap()).unwrap();
    assert_eq!(checkpoint.state, JobState::Prepared);
    assert_eq!(checkpoint.attempts, 0);
}

#[test]
fn process_gating_matches_the_feature() {
    #[cfg(feature = "external")]
    assert!(process_execution_available());
    #[cfg(not(feature = "external"))]
    assert!(!process_execution_available());
}

#[test]
fn compute_step_uses_the_input_set_registry() {
    let lattice = xtalforge_core::Lattice::cubic(5.0).unwrap();
    let structure = xtalforge_core::StructureBuilder::new("test", lattice)
        .build()
        .unwrap();
    let step = ComputeStep::new(
        "calc",
        "not-a-preset",
        EngineCommand::new("true", Vec::new()),
        "vasp_outcar",
        "OUTCAR",
    );
    let error = step.request(&structure).unwrap_err();
    assert_eq!(error.code(), "unknown_preset");
    let request = step.request_with_inputs(vec![InputFile {
        name: "INCAR".to_owned(),
        text: "x".to_owned(),
    }]);
    assert_eq!(request.inputs.names(), vec!["INCAR"]);
}

#[cfg(feature = "external")]
#[test]
fn slurm_fake_scheduler_submits_polls_and_collects() {
    let dir = temp_dir("slurm-run");
    let submitter = script(&dir, "sbatch", "echo 'Submitted batch job 4242'");
    let state = dir.join("state");
    fs::write(&state, "RUNNING\n").unwrap();
    let status = script(&dir, "squeue", &format!("cat {}", state.display()));
    let config = fake_config(Scheduler::Slurm, &submitter, &status);

    let submitted = submit(&config, &dir).unwrap();
    assert!(!submitted.dry_run);
    assert_eq!(submitted.job.id, "4242");
    assert_eq!(submitted.job.state, JobState::Submitted);

    let running = poll(&config, "4242").unwrap();
    assert_eq!(running.state, JobState::Running);
    fs::write(&state, "COMPLETED\n").unwrap();
    let completed = poll(&config, "4242").unwrap();
    assert_eq!(completed.state, JobState::Completed);

    fs::write(dir.join("run.out"), "total energy -10.0\n").unwrap();
    fs::write(dir.join("run.err"), "warning\n").unwrap();
    let collected = collect(&config, &completed, &dir).unwrap();
    assert_eq!(collected.stdout, "total energy -10.0\n");
    assert_eq!(collected.stderr, "warning\n");

    fs::write(&state, "FAILED\n").unwrap();
    assert_eq!(poll(&config, "4242").unwrap().state, JobState::Failed);
}

#[cfg(feature = "external")]
#[test]
fn pbs_fake_scheduler_maps_states() {
    let dir = temp_dir("pbs-run");
    let submitter = script(&dir, "qsub", "echo 77.server");
    let state = dir.join("state");
    fs::write(&state, "job_state = Q\n").unwrap();
    let status = script(&dir, "qstat", &format!("cat {}", state.display()));
    let config = fake_config(Scheduler::Pbs, &submitter, &status);

    let submitted = submit(&config, &dir).unwrap();
    assert_eq!(submitted.job.id, "77.server");
    assert!(submitted.script.ends_with("mof.pbs"));
    assert_eq!(poll(&config, "77.server").unwrap().state, JobState::Pending);
    for (text, expected) in [
        ("job_state = R\n", JobState::Running),
        ("job_state = C\n", JobState::Completed),
        ("job_state = F\n", JobState::Failed),
        ("job_state = H\n", JobState::Pending),
    ] {
        fs::write(&state, text).unwrap();
        assert_eq!(poll(&config, "77.server").unwrap().state, expected);
    }
}

#[cfg(feature = "external")]
#[test]
fn compute_runs_parses_and_resumes_without_relaunch() {
    let dir = temp_dir("compute-run");
    let calls = dir.join("calls");
    let engine = script(
        &dir,
        "engine",
        &format!(
            "printf 'call\\n' >> {}\ncat input.txt > result.txt",
            calls.display()
        ),
    );
    let request = compute_request(&engine, "fake");
    let mut registry = ParserRegistry::new();
    registry.register("fake", fake_parser);
    let options = RunOptions::new(&dir);

    let outcome = run_job(&request, &registry, &options).unwrap();
    assert_eq!(outcome.state, JobState::Completed);
    assert!(!outcome.resumed);
    assert_eq!(outcome.result.unwrap().energies, vec![1.5]);
    assert_eq!(fs::read_to_string(&calls).unwrap().lines().count(), 1);

    let resumed = run_job(&request, &registry, &options).unwrap();
    assert!(resumed.resumed);
    assert_eq!(resumed.state, JobState::Completed);
    assert_eq!(resumed.result.unwrap().energies, vec![1.5]);
    assert_eq!(fs::read_to_string(&calls).unwrap().lines().count(), 1);

    let checkpoint = JobCheckpoint::read(checkpoint_path(&dir, "calc").unwrap()).unwrap();
    assert_eq!(checkpoint.state, JobState::Completed);
    assert_eq!(checkpoint.attempts, 1);
    assert_eq!(checkpoint.checkpoint_version, JOB_CHECKPOINT_VERSION);
}

#[cfg(feature = "external")]
#[test]
fn failed_run_is_checkpointed_and_restarts() {
    let dir = temp_dir("compute-fail");
    let engine = script(&dir, "engine", "echo 'SCF not converged' >&2\nexit 1");
    let request = compute_request(&engine, "fake");
    let mut registry = ParserRegistry::new();
    registry.register("fake", fake_parser);
    let options = RunOptions::new(&dir);

    let failed = run_job(&request, &registry, &options).unwrap();
    assert_eq!(failed.state, JobState::Failed);
    assert_eq!(failed.failure.unwrap().code, "scf_not_converged");
    let checkpoint = JobCheckpoint::read(checkpoint_path(&dir, "calc").unwrap()).unwrap();
    assert_eq!(checkpoint.state, JobState::Failed);
    assert_eq!(checkpoint.attempts, 1);

    fs::write(&engine, "#!/bin/sh\nset -eu\ncat input.txt > result.txt\n").unwrap();
    let restarted = run_job(&request, &registry, &options).unwrap();
    assert_eq!(restarted.state, JobState::Completed);
    assert_eq!(restarted.attempts, 2);
    let checkpoint = JobCheckpoint::read(checkpoint_path(&dir, "calc").unwrap()).unwrap();
    assert_eq!(checkpoint.state, JobState::Completed);
    assert_eq!(checkpoint.attempts, 2);
}

#[cfg(feature = "external")]
#[test]
fn resume_refuses_a_different_request() {
    let dir = temp_dir("compute-mismatch");
    let engine = script(&dir, "engine", "cat input.txt > result.txt");
    let request = compute_request(&engine, "fake");
    let mut registry = ParserRegistry::new();
    registry.register("fake", fake_parser);
    let options = RunOptions::new(&dir);
    run_job(&request, &registry, &options).unwrap();

    let mut changed = request.clone();
    changed.inputs = InputFiles::one("input.txt", "payload changed");
    let error = run_job(&changed, &registry, &options).unwrap_err();
    assert_eq!(error.code(), "engine_job");
}

#[cfg(feature = "external")]
#[test]
fn unknown_parser_and_parse_failure_are_checkpointed() {
    let dir = temp_dir("compute-parse");
    let engine = script(&dir, "engine", "echo nonsense > result.txt");
    let request = compute_request(&engine, "fake");
    let options = RunOptions::new(&dir);

    let error = run_job(&request, &ParserRegistry::new(), &options).unwrap_err();
    assert_eq!(error.code(), "engine_missing");
    let checkpoint = JobCheckpoint::read(checkpoint_path(&dir, "calc").unwrap()).unwrap();
    assert_eq!(checkpoint.state, JobState::Failed);
    assert!(checkpoint.error.is_some());

    let mut registry = ParserRegistry::new();
    registry.register("fake", fake_parser);
    let error = run_job(&request, &registry, &options).unwrap_err();
    assert_eq!(error.code(), "engine_parse");
    let checkpoint = JobCheckpoint::read(checkpoint_path(&dir, "calc").unwrap()).unwrap();
    assert_eq!(checkpoint.state, JobState::Failed);
    assert_eq!(checkpoint.attempts, 2);
}
