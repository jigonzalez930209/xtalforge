//! Integration tests for opportunistic external-engine runs (Task 14.2.6).
//!
//! No test needs a real DFTB+ or xTB installation.  The probe path is
//! exercised with a missing binary, the runner path with a small shell script
//! on unix, and a real engine runs only when the probe finds one.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::sync::Arc;
use std::sync::atomic::AtomicBool;

use xtalforge_core::Structure;
use xtalforge_ff::{
    EngineKind, EngineRunOptions, EngineRunStatus, EngineStatus, PROBE_ENGINES, available_engines,
    probe, try_run_engine,
};

fn water() -> Structure {
    let points: Vec<(String, [f64; 3])> = [
        ("O", [0.0, 0.0, 0.0]),
        ("H", [0.96, 0.0, 0.0]),
        ("H", [-0.24, 0.93, 0.0]),
    ]
    .iter()
    .map(|(element, point)| ((*element).to_owned(), *point))
    .collect();
    Structure::molecule("water", &points, 8.0).expect("water")
}

fn cancel_flag() -> Arc<AtomicBool> {
    Arc::new(AtomicBool::new(false))
}

#[test]
fn probes_cover_every_known_engine() {
    assert_eq!(PROBE_ENGINES.len(), 3);
    let probes = available_engines();
    assert_eq!(probes.len(), PROBE_ENGINES.len());
    for (kind, status) in probes {
        match status {
            EngineStatus::Available { path, version } => {
                assert!(!path.is_empty(), "{kind:?} has an empty path");
                assert!(!version.is_empty(), "{kind:?} has an empty version");
            }
            EngineStatus::Missing { reason } => {
                assert!(!reason.is_empty(), "{kind:?} has an empty reason");
            }
        }
    }
}

#[test]
fn a_missing_binary_is_skipped_not_failed() {
    let options = EngineRunOptions {
        program: Some("xtalforge-no-such-engine-binary".to_owned()),
        ..EngineRunOptions::default()
    };
    let outcome = try_run_engine(EngineKind::Xtb, &water(), &options, cancel_flag()).unwrap();
    assert!(outcome.is_skipped());
    assert!(!outcome.is_completed() && !outcome.is_failed());
    assert!(
        outcome
            .reason()
            .is_some_and(|reason| reason.contains("not found"))
    );
    assert_eq!(
        EngineRunStatus::Skipped {
            reason: "missing".to_owned()
        }
        .reason(),
        Some("missing")
    );
}

#[test]
fn orca_has_no_runner() {
    let error = try_run_engine(
        EngineKind::Orca,
        &water(),
        &EngineRunOptions::default(),
        cancel_flag(),
    )
    .unwrap_err();
    assert_eq!(error.code(), "invalid_input");
    assert!(error.message().contains("ORCA"));
}

#[test]
fn invalid_engine_options_are_refused() {
    let options = EngineRunOptions {
        timeout_secs: 0.0,
        ..EngineRunOptions::default()
    };
    assert!(try_run_engine(EngineKind::Xtb, &water(), &options, cancel_flag()).is_err());
}

#[cfg(unix)]
mod scripts {
    use std::fs;
    use std::os::unix::fs::PermissionsExt;
    use std::path::PathBuf;

    use super::{cancel_flag, try_run_engine, water};
    use xtalforge_ff::{EngineKind, EngineRunOptions, HARTREE_TO_EV};

    fn script(name: &str, body: &str) -> PathBuf {
        let path = std::env::temp_dir().join(format!(
            "xtalforge-available-{}-{name}.sh",
            std::process::id()
        ));
        fs::write(&path, format!("#!/bin/sh\n{body}\n")).unwrap();
        let mut permissions = fs::metadata(&path).unwrap().permissions();
        permissions.set_mode(0o755);
        fs::set_permissions(&path, permissions).unwrap();
        path
    }

    #[test]
    fn a_present_binary_completes() {
        let program = script(
            "ok",
            "test -f coord.xyz || exit 3\n\
             cat <<'OUT'\n\
             TOTAL ENERGY -5.0 Eh\n\
             TOTAL CHARGE 0.0 e\n\
             OUT\n",
        );
        let options = EngineRunOptions {
            program: Some(program.to_string_lossy().into_owned()),
            ..EngineRunOptions::default()
        };
        let outcome = try_run_engine(EngineKind::Xtb, &water(), &options, cancel_flag()).unwrap();
        let result = outcome.result().expect("completed run");
        assert!(outcome.is_completed());
        assert!((result.energy - (-5.0 * HARTREE_TO_EV)).abs() < 1e-9);
        fs::remove_file(program).unwrap();
    }

    #[test]
    fn a_broken_binary_is_reported_as_failed() {
        let program = script("broken", "echo broken >&2\nexit 7");
        let options = EngineRunOptions {
            program: Some(program.to_string_lossy().into_owned()),
            timeout_secs: 10.0,
            ..EngineRunOptions::default()
        };
        let outcome =
            try_run_engine(EngineKind::DftbPlus, &water(), &options, cancel_flag()).unwrap();
        assert!(outcome.is_failed());
        assert!(outcome.result().is_none());
        assert!(outcome.reason().is_some_and(|reason| !reason.is_empty()));
        fs::remove_file(program).unwrap();
    }
}

#[test]
fn a_real_engine_runs_only_when_it_is_installed() {
    if let EngineStatus::Missing { reason } = probe(EngineKind::Xtb) {
        eprintln!("skipping the real xTB run: {reason}");
        return;
    }
    let options = EngineRunOptions::default();
    let outcome = try_run_engine(EngineKind::Xtb, &water(), &options, cancel_flag()).unwrap();
    match outcome {
        EngineRunStatus::Completed(result) => {
            assert!(result.energy.is_finite());
        }
        other => panic!("xTB was probed as available but the run gave {other:?}"),
    }
}
