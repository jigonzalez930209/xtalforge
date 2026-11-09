//! Integration tests for the external-engine helpers (Task 4.9).
//!
//! No test launches a real engine.  The writers are checked for their required
//! blocks and atom counts, the parsers read fabricated outputs in the units the
//! engines print, and the probe test uses a binary name that cannot exist.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::{Lattice, Site, Structure, StructureBuilder};
use xtalforge_ff::{
    DftbOptions, EngineKind, EngineStatus, OrcaOptions, XtbOptions, orca_catalogue,
    parse_dftb_output, parse_xtb_output, probe, probe_binary, write_dftb_hsd, write_orca_input,
    write_xtb_input,
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

fn periodic_nacl() -> Structure {
    let lattice = Lattice::cubic(5.0).unwrap();
    let mut builder = StructureBuilder::new("periodic-nacl", lattice);
    builder = builder.site(Site::new("Na1", "Na", [0.98, 0.5, 0.5]).unwrap());
    builder = builder.site(Site::new("Cl1", "Cl", [0.10, 0.5, 0.5]).unwrap());
    builder.build().unwrap()
}

#[test]
fn dftb_writer_contains_the_required_blocks() {
    let text = write_dftb_hsd(&water(), &DftbOptions::default()).unwrap();
    assert!(text.contains("Geometry = GenFormat"));
    assert!(text.contains("Hamiltonian = DFTB"));
    assert!(text.contains("SlaterKosterFiles = Type2FileNames"));
    assert!(text.contains("Prefix = \"slakos/\""));
    assert!(text.contains("KPointsAndWeights = SupercellFolding"));
    assert!(text.contains("SCC = Yes"));
    assert!(text.contains("2  S F"), "{text}");
    let species_line = text
        .lines()
        .find(|line| line.trim_start().starts_with('3'))
        .expect("atom count line");
    assert_eq!(species_line.split_whitespace().count(), 4, "{text}");
    assert_eq!(
        text.lines()
            .filter(
                |line| line.trim_start().starts_with(|c: char| c.is_ascii_digit())
                    && line.split_whitespace().count() == 5
            )
            .count(),
        3,
        "{text}"
    );
}

#[test]
fn dftb_parser_reads_energy_forces_and_stress() {
    let text = "\
Total energy:                       -4.0000000000 H
Total Forces
    1     0.1000000000     0.2000000000     0.0000000000
    2    -0.1000000000     0.0000000000     0.3000000000
Total stress tensor (GPa):
    1.0 0.0 0.0
    0.0 2.0 0.0
    0.0 0.0 3.0
";
    let result = parse_dftb_output(text).unwrap();
    assert!((result.energy - (-4.0 * xtalforge_ff::HARTREE_TO_EV)).abs() < 1e-9);
    assert_eq!(result.forces.len(), 2);
    assert!((result.forces[0][0] - 0.1 * 51.4220674763).abs() < 1e-6);
    let stress = result.stress.expect("stress");
    assert!((stress[0][0] - 1.0).abs() < 1e-12);
    assert!((stress[2][2] - 3.0).abs() < 1e-12);
}

#[test]
fn dftb_parser_allows_missing_forces_but_not_energy() {
    let result = parse_dftb_output("Total energy: -1.0 H\n").unwrap();
    assert!(result.forces.is_empty());
    assert!(result.stress.is_none());
    assert!(parse_dftb_output("no energy here\n").is_err());
}

#[test]
fn xtb_writer_contains_xyz_and_the_blocks() {
    let text = write_xtb_input(&water(), &XtbOptions::default()).unwrap();
    assert!(text.starts_with("3\n"), "{text}");
    assert!(text.contains("$coord"));
    assert!(text.contains("$end"));
    assert!(text.contains("$chrg 0.000000"));
    assert!(text.contains("$uhf 0"));
    assert!(text.contains("method=GFN2-xTB"));
    let coord_start = text.find("$coord").unwrap();
    let coord_end = text.find("$end").unwrap();
    let rows = text[coord_start..coord_end]
        .lines()
        .skip(1)
        .filter(|line| !line.trim().is_empty())
        .count();
    assert_eq!(rows, 3);
}

#[test]
fn xtb_parser_reads_energy_charge_and_forces() {
    let text = "\
TOTAL ENERGY        -5.000000000000 Eh
TOTAL CHARGE        -0.000000000000 e
:: forces (Eh/a0) ::
    1     0.1000000000     0.0000000000     0.0000000000
    2     0.0000000000    -0.2000000000     0.0000000000
";
    let result = parse_xtb_output(text).unwrap();
    assert!((result.energy - (-5.0 * xtalforge_ff::HARTREE_TO_EV)).abs() < 1e-9);
    assert_eq!(result.forces.len(), 2);
    assert!((result.forces[1][1] + 0.2 * 51.4220674763).abs() < 1e-6);
    assert!(result.stress.is_none());
    assert!(parse_xtb_output("TOTAL ENERGY nonsense Eh\n").is_err());
}

#[test]
fn orca_writer_contains_the_simple_input_and_atom_count() {
    let text = write_orca_input(&water(), &OrcaOptions::default()).unwrap();
    assert!(text.starts_with("! B3LYP def2-SVP\n"), "{text}");
    assert!(text.contains("* xyz 0 1"));
    let block = text
        .lines()
        .skip_while(|line| !line.starts_with("* xyz"))
        .skip(1)
        .take_while(|line| *line != "*");
    assert_eq!(block.count(), 3);
}

#[test]
fn orca_writer_makes_a_molecule_whole_across_the_cell_face() {
    let structure = periodic_nacl();
    let text = write_orca_input(&structure, &OrcaOptions::default()).unwrap();
    let coordinates: Vec<[f64; 3]> = text
        .lines()
        .skip_while(|line| !line.starts_with("* xyz"))
        .skip(1)
        .take_while(|line| *line != "*")
        .map(|line| {
            let values: Vec<f64> = line
                .split_whitespace()
                .skip(1)
                .map(|token| token.parse::<f64>().unwrap())
                .collect();
            [values[0], values[1], values[2]]
        })
        .collect();
    assert_eq!(coordinates.len(), 2);
    let delta = [
        coordinates[1][0] - coordinates[0][0],
        coordinates[1][1] - coordinates[0][1],
        coordinates[1][2] - coordinates[0][2],
    ];
    let distance = (delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2]).sqrt();
    assert!(
        distance < 1.0,
        "bonded Na-Cl images should stay together, got {distance}"
    );
}

#[test]
fn orca_catalogue_offers_functionals_and_bases() {
    let catalogue = orca_catalogue();
    assert!(catalogue.len() >= 8);
    assert!(catalogue.contains(&("B3LYP", "functional")));
    assert!(catalogue.contains(&("def2-SVP", "basis")));
}

#[test]
fn engine_kinds_serialize_with_their_documented_names() {
    assert_eq!(
        serde_json::to_string(&EngineKind::DftbPlus).unwrap(),
        "\"dftb+\""
    );
    assert_eq!(serde_json::to_string(&EngineKind::Xtb).unwrap(), "\"xtb\"");
    assert_eq!(
        serde_json::to_string(&EngineKind::Orca).unwrap(),
        "\"orca\""
    );
}

#[test]
fn missing_binary_probes_report_missing_and_kinds_do_not_panic() {
    let status = probe_binary("xtalforge-no-such-engine-4-9");
    assert!(matches!(status, EngineStatus::Missing { .. }));
    for kind in [EngineKind::DftbPlus, EngineKind::Xtb, EngineKind::Orca] {
        let status = probe(kind);
        match status {
            EngineStatus::Available { path, version } => {
                assert!(!path.is_empty() && !version.is_empty());
            }
            EngineStatus::Missing { reason } => assert!(!reason.is_empty()),
        }
    }
}

#[test]
fn writers_reject_bad_inputs() {
    let bad_dftb = DftbOptions {
        kpoints: [1, 0, 1],
        ..DftbOptions::default()
    };
    assert!(write_dftb_hsd(&water(), &bad_dftb).is_err());
    let bad_xtb = XtbOptions {
        charge: f64::NAN,
        ..XtbOptions::default()
    };
    assert!(write_xtb_input(&water(), &bad_xtb).is_err());
    let bad_orca = OrcaOptions {
        multiplicity: 0,
        ..OrcaOptions::default()
    };
    assert!(write_orca_input(&water(), &bad_orca).is_err());
    let injected = OrcaOptions {
        method: "B3LYP\n* xyz".to_owned(),
        ..OrcaOptions::default()
    };
    assert!(write_orca_input(&water(), &injected).is_err());
}

#[test]
fn the_recorded_dftb_output_parses() {
    let text = include_str!("fixtures/dftb/si_scc_detailed.out");
    let result = parse_dftb_output(text).unwrap();
    assert!(
        (result.energy - (-69.2989)).abs() < 1e-3,
        "{}",
        result.energy
    );
    assert_eq!(result.forces.len(), 2);
    let sum = result.forces.iter().fold([0.0; 3], |total, force| {
        [
            total[0] + force[0],
            total[1] + force[1],
            total[2] + force[2],
        ]
    });
    assert!(sum.iter().all(|value| value.abs() < 1e-6));
    let stress = result.stress.expect("a stress tensor");
    assert!(stress.iter().flatten().all(|value| value.is_finite()));
}
