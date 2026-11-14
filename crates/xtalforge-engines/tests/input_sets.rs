//! Integration tests for Task 18.1 input sets.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::indexing_slicing,
    clippy::panic
)]

use std::collections::BTreeSet;

use serde_json::{Value, json};
use xtalforge_core::{Lattice, Site, Structure, StructureBuilder};
use xtalforge_engines::Error;
use xtalforge_engines::input_set::{
    InputSetRequest, InputSetSummary, list_presets, preset, preset_description, preset_provenance,
    write_input_set,
};
use xtalforge_ff::OrcaOptions;
use xtalforge_ff::external::write_orca_input;

const RASPA_SIMULATION: &str = include_str!("fixtures/input-sets/raspa-simulation.input");
const RASPA_FRAMEWORK: &str = include_str!("fixtures/input-sets/raspa-framework.cif");
const RASPA_FORCE_FIELD: &str = include_str!("fixtures/input-sets/raspa-force-field");
const RASPA_MIXING: &str = include_str!("fixtures/input-sets/raspa-mixing-rules");

fn nacl() -> Structure {
    StructureBuilder::new("NaCl", Lattice::cubic(5.64).unwrap())
        .site(Site::new("Na1", "Na", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("Cl1", "Cl", [0.5, 0.5, 0.5]).unwrap())
        .build()
        .unwrap()
}

fn charged_nacl() -> Structure {
    StructureBuilder::new("NaCl", Lattice::cubic(5.64).unwrap())
        .site(
            Site::new("Na1", "Na", [0.0, 0.0, 0.0])
                .unwrap()
                .with_charge(1.0)
                .unwrap(),
        )
        .site(
            Site::new("Cl1", "Cl", [0.5, 0.5, 0.5])
                .unwrap()
                .with_charge(-1.0)
                .unwrap(),
        )
        .build()
        .unwrap()
}

fn feo() -> Structure {
    let iron = Site::new("Fe1", "Fe", [0.0, 0.0, 0.0])
        .unwrap()
        .with_prop("magmom", json!(4.0))
        .with_prop("ldau_u", json!(4.0))
        .with_prop("ldau_orbital", json!("d"));
    let oxygen = Site::new("O1", "O", [0.5, 0.5, 0.5])
        .unwrap()
        .with_prop("magmom", json!(0.0))
        .with_prop("ldau_u", json!(0.0))
        .with_prop("ldau_orbital", json!("p"));
    StructureBuilder::new("FeO", Lattice::cubic(4.33).unwrap())
        .site(iron)
        .site(oxygen)
        .build()
        .unwrap()
}

fn request(preset: &str, overrides: Value) -> InputSetRequest {
    InputSetRequest {
        preset: preset.to_owned(),
        overrides: match overrides {
            Value::Object(map) => map,
            _ => panic!("overrides must be an object"),
        },
    }
}

fn build(structure: &Structure, preset: &str, overrides: Value) -> InputSetSummary {
    write_input_set(structure, &request(preset, overrides)).unwrap()
}

fn text<'a>(summary: &'a InputSetSummary, name: &str) -> &'a str {
    summary
        .files
        .iter()
        .find(|file| file.name == name)
        .map(|file| file.text.as_str())
        .expect("expected input file")
}

#[test]
fn registry_is_real_and_every_preset_builds_byte_identically() {
    let names = list_presets();
    let unique: BTreeSet<&str> = names.iter().copied().collect();
    assert_eq!(names.len(), unique.len());
    assert_eq!(names.len(), 20);
    let structure = nacl();
    for name in names {
        let writer = preset(name).expect("preset must be registered");
        let first = writer.build(&structure, &Default::default()).unwrap();
        let second = writer.build(&structure, &Default::default()).unwrap();
        assert_eq!(first, second, "{name} is not deterministic");
        assert_eq!(first.preset, name);
        assert!(!first.files.is_empty());
        assert!(preset_description(name).is_some());
        let provenance = preset_provenance(name).expect("provenance");
        assert_eq!(provenance.version, 1);
        assert!(!provenance.source.is_empty());
        assert!(!provenance.reference.is_empty());
        assert!(!provenance.license.is_empty());
    }
    assert!(preset("vasp/unknown").is_none());
}

#[test]
fn vasp_presets_write_incar_kpoints_poscar_and_potcar_spec() {
    let structure = nacl();
    for name in [
        "vasp/relax",
        "vasp/static",
        "vasp/band",
        "vasp/dos",
        "vasp/mof",
    ] {
        let summary = build(&structure, name, json!({}));
        let files: Vec<&str> = summary
            .files
            .iter()
            .map(|file| file.name.as_str())
            .collect();
        assert_eq!(files, ["INCAR", "KPOINTS", "POSCAR", "POTCAR.spec"]);
        assert!(text(&summary, "POSCAR").contains("Direct"));
        assert!(text(&summary, "POTCAR.spec").contains("PAW_PBE Na_pv"));
        assert!(text(&summary, "POTCAR.spec").contains("PAW_PBE Cl"));
        assert!(!text(&summary, "POTCAR.spec").contains("PAW_PBE PAW"));
    }
    let relax = build(&structure, "vasp/relax", json!({}));
    assert!(text(&relax, "INCAR").contains("ISIF = 3"));
    assert!(text(&relax, "INCAR").contains("IVDW = 12"));
    assert!(text(&relax, "POSCAR").contains("Selective dynamics"));
    let static_run = build(&structure, "vasp/static", json!({}));
    assert!(text(&static_run, "INCAR").contains("NSW = 0"));
    assert!(text(&static_run, "INCAR").contains("LCHARG = .TRUE."));
    let band = build(&structure, "vasp/band", json!({}));
    assert!(text(&band, "KPOINTS").contains("Line-mode"));
    assert!(text(&band, "KPOINTS").contains("! G"));
    let dos = build(&structure, "vasp/dos", json!({}));
    assert!(text(&dos, "INCAR").contains("ISMEAR = -5"));
}

#[test]
fn vasp_kpoint_density_follows_the_target_spacing() {
    let structure = nacl();
    let coarse = build(
        &structure,
        "vasp/relax",
        json!({"kpoints": {"spacing": 0.2}}),
    );
    let dense = build(
        &structure,
        "vasp/relax",
        json!({"kpoints": {"spacing": 0.1}}),
    );
    let lines_coarse: Vec<&str> = text(&coarse, "KPOINTS").lines().collect();
    let lines_dense: Vec<&str> = text(&dense, "KPOINTS").lines().collect();
    assert_eq!(lines_coarse[3], "6 6 6");
    assert_eq!(lines_dense[3], "12 12 12");
    let explicit = build(
        &structure,
        "vasp/relax",
        json!({"kpoints": {"mode": "mesh", "mesh": [3, 4, 5]}}),
    );
    assert!(text(&explicit, "KPOINTS").contains("3 4 5"));
}

#[test]
fn vasp_magmom_and_ldau_come_from_site_properties() {
    let summary = build(&feo(), "vasp/mof", json!({}));
    let incar = text(&summary, "INCAR");
    assert!(incar.contains("ISPIN = 2"), "{incar}");
    assert!(incar.contains("MAGMOM = 4.0000 0.0000"), "{incar}");
    assert!(incar.contains("LDAU = .TRUE."));
    assert!(incar.contains("LDAUL = 2 1"), "{incar}");
    assert!(incar.contains("LDAUU = 4.0000 0.0000"), "{incar}");
    assert!(incar.contains("LDAUJ = 0.0000 0.0000"));
}

#[test]
fn vasp_metadata_is_the_documented_magmom_fallback() {
    let structure = StructureBuilder::new("FeO", Lattice::cubic(4.33).unwrap())
        .site(Site::new("Fe1", "Fe", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("O1", "O", [0.5, 0.5, 0.5]).unwrap())
        .metadata("magmom", json!([5.0, 0.0]))
        .metadata("ldau_u", json!([5.0, 0.0]))
        .build()
        .unwrap();
    let summary = build(&structure, "vasp/relax", json!({}));
    let incar = text(&summary, "INCAR");
    assert!(incar.contains("MAGMOM = 5.0000 0.0000"), "{incar}");
    assert!(incar.contains("LDAUU = 5.0000 0.0000"), "{incar}");
}

#[test]
fn overrides_merge_deeply_and_reject_unknown_keys() {
    let structure = nacl();
    let summary = build(
        &structure,
        "vasp/relax",
        json!({"encut": 600.0, "incar": {"NCORE": 4}}),
    );
    let incar = text(&summary, "INCAR");
    assert!(incar.contains("ENCUT = 600"));
    assert!(incar.contains("NCORE = 4"));

    let error = write_input_set(
        &structure,
        &request("vasp/relax", json!({"kpoints": {"spacings": 0.2}})),
    )
    .unwrap_err();
    assert!(matches!(error, Error::UnknownOverride(ref key) if key == "kpoints.spacings"));
    assert_eq!(error.code(), "unknown_override");

    let error =
        write_input_set(&structure, &request("vasp/relax", json!({"not_a_tag": 1}))).unwrap_err();
    assert!(matches!(error, Error::UnknownOverride(ref key) if key == "not_a_tag"));

    let allowed = build(
        &structure,
        "vasp/relax",
        json!({"allow_unknown": true, "not_a_tag": {"nested": 1}}),
    );
    let baseline = build(&structure, "vasp/relax", json!({}));
    assert_eq!(allowed, baseline);
}

#[test]
fn quantum_espresso_presets_and_pseudopotentials() {
    let structure = feo();
    let scf = build(&structure, "qe/scf", json!({}));
    assert!(text(&scf, "pw.in").contains("calculation = 'scf'"));
    assert!(text(&scf, "pw.in").contains("K_POINTS automatic"));
    assert!(text(&scf, "pw.in").contains("Fe.pbe-spn-kjpaw_psl.1.0.0.UPF"));
    let relax = build(&structure, "qe/relax", json!({}));
    assert!(text(&relax, "pw.in").contains("calculation = 'relax'"));
    assert!(text(&relax, "pw.in").contains("&IONS"));
    let vc = build(&structure, "qe/vc-relax", json!({}));
    assert!(text(&vc, "pw.in").contains("calculation = 'vc-relax'"));
    assert!(text(&vc, "pw.in").contains("&CELL"));
    let nscf = build(&structure, "qe/nscf", json!({}));
    assert!(text(&nscf, "pw.in").contains("calculation = 'nscf'"));
    let bands = build(&structure, "qe/bands", json!({}));
    assert!(text(&bands, "pw.in").contains("K_POINTS crystal_b"));
    assert!(!text(&bands, "pw.in").contains("&IONS"));
    let mapped = build(
        &structure,
        "qe/scf",
        json!({"pseudopotentials": {"Fe": "Fe.custom.UPF"}}),
    );
    assert!(text(&mapped, "pw.in").contains("Fe.custom.UPF"));
}

#[test]
fn cp2k_presets_write_gpw_input_with_molopt_and_d3() {
    let structure = feo();
    let geo = build(&structure, "cp2k/geo-opt", json!({}));
    let input = text(&geo, "cp2k.inp");
    assert!(input.contains("RUN_TYPE GEO_OPT"));
    assert!(input.contains("METHOD Quickstep"));
    assert!(input.contains("TYPE DFTD3"));
    assert!(input.contains("BASIS_SET DZVP-MOLOPT-SR-GTH"));
    assert!(input.contains("POTENTIAL GTH-PBE-q16"));
    assert!(input.contains("POTENTIAL GTH-PBE-q6"));
    assert!(input.contains("SCALED"));
    let single = build(&structure, "cp2k/single-point", json!({}));
    assert!(text(&single, "cp2k.inp").contains("RUN_TYPE ENERGY"));
    assert!(!text(&single, "cp2k.inp").contains("&GEO_OPT"));
}

#[test]
fn lammps_writes_data_and_uff_script() {
    let structure = nacl();
    for name in ["lammps/uff", "lammps/uff4mof"] {
        let summary = build(&structure, name, json!({}));
        let files: Vec<&str> = summary
            .files
            .iter()
            .map(|file| file.name.as_str())
            .collect();
        assert_eq!(files, ["system.data", "input.lammps"]);
        assert!(text(&summary, "system.data").contains("LAMMPS data file"));
        let script = text(&summary, "input.lammps");
        assert!(script.contains("read_data system.data"));
        assert!(script.contains("pair_style lj/cut 12.0000"));
        assert!(script.contains("pair_coeff 1 1"));
        assert!(script.contains("pair_coeff 1 2"));
        assert!(script.contains("minimize"));
    }
    let single = build(
        &structure,
        "lammps/uff",
        json!({"mode": "run", "steps": 10}),
    );
    assert!(text(&single, "input.lammps").contains("run 10"));
    let uff = build(&structure, "lammps/uff", json!({}));
    let uff4mof = build(&structure, "lammps/uff4mof", json!({}));
    assert_ne!(text(&uff, "input.lammps"), text(&uff4mof, "input.lammps"));
}

#[test]
fn gulp_writes_lattice_relaxation_input() {
    let summary = build(&charged_nacl(), "gulp/lattice-relax", json!({}));
    let input = text(&summary, "gulp.gin");
    assert!(input.contains("opti conp"));
    assert!(input.contains("vectors"));
    assert!(input.contains("fractional"));
    assert!(input.contains("species"));
    assert!(input.contains("Na core 1.0000"));
    assert!(input.contains("Cl core -1.0000"));
    assert!(input.contains("uff"));
    let relax = build(&charged_nacl(), "gulp/relax", json!({}));
    assert!(text(&relax, "gulp.gin").contains("opti\n"));
}

#[test]
fn raspa_matches_the_phase_10_export_byte_for_byte() {
    let summary = build(&nacl(), "raspa/widom", json!({}));
    assert_eq!(text(&summary, "simulation.input"), RASPA_SIMULATION);
    assert_eq!(text(&summary, "framework.cif"), RASPA_FRAMEWORK);
    assert_eq!(text(&summary, "force_field"), RASPA_FORCE_FIELD);
    assert_eq!(text(&summary, "mixing_rules"), RASPA_MIXING);
}

#[test]
fn orca_matches_the_existing_writer_byte_for_byte() {
    let structure = nacl();
    let summary = build(&structure, "orca/single-point", json!({}));
    let expected = write_orca_input(&structure, &OrcaOptions::default()).unwrap();
    assert_eq!(text(&summary, "orca.inp"), expected);
    let pbe0 = build(&structure, "orca/pbe0-def2-tzvp", json!({}));
    assert!(text(&pbe0, "orca.inp").starts_with("! PBE0 def2-TZVP"));
}

#[test]
fn invalid_settings_return_errors_instead_of_panicking() {
    let structure = nacl();
    let error = write_input_set(
        &structure,
        &request("vasp/relax", json!({"encut": "large"})),
    )
    .unwrap_err();
    assert_eq!(error.code(), "engine_job");
    let error = write_input_set(
        &structure,
        &request("vasp/relax", json!({"kpoints": {"spacing": 0.0}})),
    )
    .unwrap_err();
    assert_eq!(error.code(), "engine_job");
    let error = write_input_set(
        &structure,
        &request("vasp/relax", json!({"allow_unknown": "yes"})),
    )
    .unwrap_err();
    assert_eq!(error.code(), "engine_job");
    let error = write_input_set(
        &structure,
        &request("qe/scf", json!({"calculation": "unknown"})),
    )
    .unwrap_err();
    assert_eq!(error.code(), "engine_job");
    let error = write_input_set(
        &structure,
        &request("cp2k/geo-opt", json!({"run_type": "UNKNOWN"})),
    )
    .unwrap_err();
    assert_eq!(error.code(), "engine_job");
    let error = write_input_set(
        &structure,
        &request("raspa/widom", json!({"guest": "unobtainium"})),
    )
    .unwrap_err();
    assert_eq!(error.code(), "engine_job");
    let error = write_input_set(
        &structure,
        &request("lammps/uff", json!({"parameter_set": "unknown"})),
    )
    .unwrap_err();
    assert_eq!(error.code(), "engine_job");
    let error = write_input_set(
        &structure,
        &request("orca/single-point", json!({"multiplicity": 0})),
    )
    .unwrap_err();
    assert_eq!(error.code(), "engine_job");
    let error = write_input_set(&structure, &request("vasp/missing", json!({}))).unwrap_err();
    assert_eq!(error.code(), "unknown_preset");
}

#[test]
fn empty_structures_are_rejected_by_the_dft_writers() {
    let empty = StructureBuilder::new("empty", Lattice::cubic(5.0).unwrap())
        .build()
        .unwrap();
    for name in [
        "vasp/relax",
        "qe/scf",
        "cp2k/geo-opt",
        "lammps/uff",
        "gulp/relax",
    ] {
        let error = write_input_set(&empty, &request(name, json!({}))).unwrap_err();
        assert_eq!(error.code(), "engine_job", "{name}");
    }
}
