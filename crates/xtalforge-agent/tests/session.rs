//! Session, documents and verbs (Tasks 8.1.2 to 8.1.5).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use serde_json::json;
use xtalforge_agent::answers::Answer;
use xtalforge_agent::{Error, Session};

const NACL: &str = include_str!("../../../tests/corpus/simple/NaCl.cif");
const DIAMOND: &str = include_str!("../../../tests/corpus/simple/diamond.cif");

fn loaded() -> (Session, String) {
    let mut session = Session::new();
    let answer = session
        .run(
            "load",
            &json!({"text": NACL, "format": "cif", "name": "nacl"}),
        )
        .expect("load");
    let Answer::Loaded { document } = answer else {
        panic!("expected a loaded answer");
    };
    (session, document.name)
}

#[test]
fn capabilities_report_the_catalogue() {
    let mut session = Session::new();
    let answer = session
        .run("capabilities", &json!({}))
        .expect("capabilities");
    let Answer::Capabilities { capabilities } = answer else {
        panic!("expected capabilities");
    };
    assert_eq!(
        capabilities.verbs.len(),
        xtalforge_agent::verbs::VERBS.len()
    );
    assert!(capabilities.features.contains_key("sci-form"));
    assert!(
        capabilities
            .formats
            .iter()
            .any(|format| format.name == "cif")
    );
}

#[test]
fn load_info_and_diagnostics() {
    let (mut session, _) = loaded();
    let answer = session.run("info", &json!({})).expect("info");
    let Answer::Info {
        document,
        inspection,
        diagnostics,
    } = answer
    else {
        panic!("expected info");
    };
    assert_eq!(document, "nacl");
    assert_eq!(inspection.summary.formula, "Na4 Cl4");
    assert!(inspection.n_p1_atoms >= 2);
    assert_eq!(inspection.sites.len(), inspection.summary.n_sites);
    assert_eq!(
        diagnostics.counts.values().sum::<usize>(),
        diagnostics.diagnostics.len()
    );
    assert!(!diagnostics.has_errors(), "NaCl must not have errors");
}

#[test]
fn unknown_verb_and_document_are_typed_errors() {
    let mut session = Session::new();
    let error = session.run("teleport", &json!({})).unwrap_err();
    assert_eq!(error.code(), "unknown_verb");
    let error = session.run("info", &json!({})).unwrap_err();
    assert_eq!(error.code(), "no_document");
    let (mut session, _) = loaded();
    let error = session
        .run("info", &json!({"document": "missing"}))
        .unwrap_err();
    assert_eq!(error.code(), "unknown_document");
}

#[test]
fn supercell_undo_and_change_bursts() {
    let (mut session, _) = loaded();
    let before = session
        .active()
        .expect("document")
        .structure()
        .sites()
        .len();
    let answer = session
        .run("supercell", &json!({"repeats": [2, 2, 2]}))
        .expect("supercell");
    let Answer::Structure { summary, .. } = answer else {
        panic!("expected a structure");
    };
    assert_eq!(summary.n_sites, before * 8);
    let answer = session.run("changes", &json!({})).expect("changes");
    let Answer::Changes { bursts, .. } = answer else {
        panic!("expected changes");
    };
    assert_eq!(bursts.len(), 1);
    assert_eq!(bursts[0].len, 1);
    assert_eq!(bursts[0].undo_steps, Some(1));
    let answer = session.run("undo", &json!({})).expect("undo");
    let Answer::History { applied, .. } = answer else {
        panic!("expected history");
    };
    assert!(applied);
    assert_eq!(
        session
            .active()
            .expect("document")
            .structure()
            .sites()
            .len(),
        before
    );
    let answer = session.run("redo", &json!({})).expect("redo");
    let Answer::History { applied, .. } = answer else {
        panic!("expected history");
    };
    assert!(applied);
}

#[test]
fn history_records_the_agent_origin() {
    let (mut session, _) = loaded();
    session.run("wrap", &json!({})).expect("wrap");
    let answer = session.run("history", &json!({})).expect("history");
    let Answer::History {
        entries,
        can_undo,
        can_redo,
        ..
    } = answer
    else {
        panic!("expected history");
    };
    assert!(can_undo);
    assert!(!can_redo);
    assert_eq!(entries.last().expect("entry").origin, "agent");
}

#[test]
fn save_round_trips_through_a_file() {
    let (mut session, _) = loaded();
    let mut path = std::env::temp_dir();
    path.push(format!("xtalforge-agent-{}.cif", std::process::id()));
    let answer = session
        .run("save", &json!({"path": path.to_string_lossy()}))
        .expect("save");
    let Answer::Saved { format, bytes, .. } = answer else {
        panic!("expected saved");
    };
    assert_eq!(format, "cif");
    assert!(bytes > 0);
    let text = std::fs::read_to_string(&path).expect("written file");
    assert!(text.contains("_cell_length_a"));
    let _ = std::fs::remove_file(&path);
}

#[test]
fn symmetry_bonds_and_rings() {
    let (mut session, _) = loaded();
    let answer = session.run("symmetry", &json!({})).expect("symmetry");
    let Answer::Symmetry { symmetry, .. } = answer else {
        panic!("expected symmetry");
    };
    assert_eq!(symmetry.number, 225, "NaCl is Fm-3m");
    assert_eq!(symmetry.hm_symbol, "F m -3 m");
    let answer = session.run("bonds", &json!({})).expect("bonds");
    let Answer::Bonds { bonds, .. } = answer else {
        panic!("expected bonds");
    };
    assert!(!bonds.is_empty());
    assert!(bonds.iter().all(|bond| bond.distance > 0.0));

    let mut session = Session::new();
    session
        .run(
            "load",
            &json!({"text": DIAMOND, "format": "cif", "name": "diamond"}),
        )
        .expect("load diamond");
    let answer = session
        .run("rings", &json!({"max_size": 8}))
        .expect("rings");
    let Answer::Rings { total, counts, .. } = answer else {
        panic!("expected rings");
    };
    assert!(total > 0, "diamond has six-membered rings");
    assert!(counts.iter().any(|count| count.size == 6));
}

#[test]
fn pxrd_and_porosity_are_finite() {
    let (mut session, _) = loaded();
    let answer = session.run("pxrd", &json!({})).expect("pxrd");
    let Answer::Pxrd {
        n_reflections,
        reflections,
        points,
        ..
    } = answer
    else {
        panic!("expected pxrd");
    };
    assert!(n_reflections > 0);
    assert!(!reflections.is_empty());
    assert!(!points.is_empty());
    assert!(points.iter().all(|point| point.intensity.is_finite()));
    let answer = session
        .run("porosity", &json!({"spacing": 0.6}))
        .expect("porosity");
    let Answer::Porosity { porosity, .. } = answer else {
        panic!("expected porosity");
    };
    assert!(porosity.cell_volume > 0.0);
    assert!(porosity.accessible_volume >= 0.0);
    assert!(porosity.accessible_fraction <= 1.0 + 1e-9);
    assert!(porosity.largest_included_diameter >= 0.0);
}

#[test]
fn uff_energy_is_finite_and_optimize_changes_the_geometry() {
    let (mut session, _) = loaded();
    let answer = session.run("uff_energy", &json!({})).expect("energy");
    let Answer::Energy {
        energy_ev,
        max_force_ev_per_angstrom,
        ..
    } = answer
    else {
        panic!("expected energy");
    };
    assert!(energy_ev.is_finite());
    assert!(max_force_ev_per_angstrom.is_finite());

    // A distorted water molecule relaxes to a lower energy.
    let water = "3\nwater\nO 0.0 0.0 0.0\nH 1.10 0.0 0.0\nH -0.30 0.95 0.0\n";
    let mut session = Session::new();
    session
        .run(
            "load",
            &json!({"text": water, "format": "xyz", "name": "water"}),
        )
        .expect("load water");
    let before = session.run("uff_energy", &json!({})).expect("energy");
    let Answer::Energy {
        energy_ev: first, ..
    } = before
    else {
        panic!("expected energy");
    };
    let answer = session
        .run("optimize", &json!({"max_steps": 60}))
        .expect("optimize");
    let Answer::Optimization { optimization, .. } = answer else {
        panic!("expected optimization");
    };
    assert!(optimization.steps > 0);
    assert!(
        optimization.energy_ev <= first + 1e-6,
        "optimization raised the energy: {} -> {}",
        first,
        optimization.energy_ev
    );
}

#[test]
fn size_limits_are_checked_before_allocating() {
    let (mut session, _) = loaded();
    let error = session
        .run("supercell", &json!({"repeats": [1000, 1000, 1000]}))
        .unwrap_err();
    assert_eq!(error.code(), "limit_exceeded");
    assert_eq!(
        session
            .active()
            .expect("document")
            .structure()
            .sites()
            .len(),
        2
    );
}

#[test]
fn a_script_runs_against_one_session() {
    let script = serde_json::json!({
        "verbs": [
            {"verb": "load", "arguments": {"text": NACL, "format": "cif"}},
            {"verb": "info", "arguments": {}},
            {"verb": "symmetry", "arguments": {}}
        ]
    });
    let answers = xtalforge_agent::run_verbs(&script.to_string()).expect("script");
    assert_eq!(answers.len(), 3);
    assert_eq!(answers[0].kind(), "loaded");
    assert_eq!(answers[1].kind(), "info");
    assert_eq!(answers[2].kind(), "symmetry");
    let error = xtalforge_agent::run_verbs("not json").unwrap_err();
    assert_eq!(error.code(), "json");
}

#[test]
fn close_removes_the_document() {
    let (mut session, _) = loaded();
    let answer = session.run("close", &json!({})).expect("close");
    let Answer::Documents { documents, active } = answer else {
        panic!("expected documents");
    };
    assert!(documents.is_empty());
    assert!(active.is_none());
    assert!(matches!(
        session.run("info", &json!({})),
        Err(Error::NoDocument)
    ));
}
