//! Integration tests for the MOF platform (Phase 10).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::collections::BTreeMap;

use xtalforge_core::{Lattice, Site, Structure, StructureBuilder};
use xtalforge_mof::{
    CapSpecies, DefectOptions, DisorderOptions, FunctionalGroup, FunctionalizeOptions, GcmcOptions,
    MixtureOptions, OpenMetalOptions, PoreOptions, ScreeningOptions, SolventOptions, WidomOptions,
    charge_balance, check_mof, deconstruct, deduplicate, energy_grid, framework_descriptors,
    functionalize, gcmc_isotherm, interpenetrate, mixed_framework, mofid, open_metal_sites,
    pore_descriptors, raspa_input, read_core_mof, remove_solvent, resolve_disorder,
    screen_database, stack_layers, widom_insertion,
};

fn framework() -> Structure {
    // A coarse pcu-like Zn-O-C rod: bonds are perceived on the 1.8-2.1 A
    // contacts; nodes and linkers are small but topologically connected.
    StructureBuilder::new("test-mof", Lattice::cubic(12.0).unwrap())
        .site(Site::new("Zn1", "Zn", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("O1", "O", [0.125, 0.0, 0.0]).unwrap())
        .site(Site::new("C1", "C", [0.208, 0.0, 0.0]).unwrap())
        .site(Site::new("C2", "C", [0.292, 0.0, 0.0]).unwrap())
        .site(Site::new("O2", "O", [0.375, 0.0, 0.0]).unwrap())
        .site(Site::new("Zn2", "Zn", [0.5, 0.0, 0.0]).unwrap())
        .build()
        .unwrap()
}

fn porous() -> Structure {
    StructureBuilder::new("porous", Lattice::cubic(20.0).unwrap())
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .build()
        .unwrap()
}

#[test]
fn checker_flags_overlaps_and_reports_clean_frameworks() {
    let clean = check_mof(&framework()).unwrap();
    assert!(clean.clean || clean.counts.get("warning").copied().unwrap_or(0) > 0);
    let overlapping = StructureBuilder::new("bad", Lattice::cubic(10.0).unwrap())
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("C2", "C", [0.001, 0.0, 0.0]).unwrap())
        .build()
        .unwrap();
    let report = check_mof(&overlapping).unwrap();
    assert!(
        report
            .findings
            .iter()
            .any(|finding| finding.code == "overlapping_atoms")
    );
}

#[test]
fn solvent_removal_drops_free_water() {
    let mut builder = StructureBuilder::new("wet", Lattice::cubic(12.0).unwrap());
    for site in framework().sites() {
        builder = builder.site(site.clone());
    }
    let with_water = builder
        .site(Site::new("O2", "O", [0.7, 0.7, 0.7]).unwrap())
        .site(Site::new("H1", "H", [0.78, 0.7, 0.7]).unwrap())
        .site(Site::new("H2", "H", [0.7, 0.78, 0.7]).unwrap())
        .build()
        .unwrap();
    let report = remove_solvent(&with_water, &SolventOptions::default()).unwrap();
    assert_eq!(report.removed.len(), 1);
    assert_eq!(report.structure.sites().len(), 6);
}

#[test]
fn deconstruction_and_mofid_are_non_empty() {
    let deconstruction = deconstruct(&framework()).unwrap();
    assert!(!deconstruction.units.is_empty());
    let id = mofid(&framework()).unwrap();
    assert!(!id.mofid.is_empty());
    assert_eq!(id.mofkey.len(), 64);
    let descriptors = framework_descriptors(&framework()).unwrap();
    assert_eq!(descriptors.n_atoms, 6);
    assert!(!descriptors.formula.is_empty());
}

#[test]
fn porosity_and_open_metal_sites_run() {
    let pores = pore_descriptors(
        &porous(),
        &PoreOptions {
            spacing: 1.0,
            ..PoreOptions::default()
        },
    )
    .unwrap();
    assert!(pores.descriptors.lcd > 0.0);
    assert!(pores.descriptors.accessible_volume > 0.0);
    assert!(pores.descriptors.channel_dimension <= 3);
    let sites = open_metal_sites(&framework(), &OpenMetalOptions::default()).unwrap();
    assert!(
        sites
            .iter()
            .all(|site| site.expected_coordination >= site.coordination)
    );
}

#[test]
fn widom_gcmc_grid_and_raspa_are_finite() {
    let widom = widom_insertion(
        &porous(),
        "CO2",
        &WidomOptions {
            temperature: 298.0,
            samples: 100,
            seed: 1,
        },
    )
    .unwrap();
    assert!(widom.henry.is_finite());
    let isotherm = gcmc_isotherm(
        &porous(),
        "CO2",
        &GcmcOptions {
            temperature: 298.0,
            pressures: vec![1.0, 5.0],
            steps: 20,
            seed: 1,
        },
    )
    .unwrap();
    assert_eq!(isotherm.points.len(), 2);
    let grid = energy_grid(&porous(), "CO2", 2.0).unwrap();
    assert!(grid.maximum >= grid.minimum);
    let raspa = raspa_input(&porous(), "CO2", 298.0, 1.0).unwrap();
    assert!(raspa.simulation_input.contains("MonteCarlo"));
    assert!(raspa.framework_cif.contains("data_"));
}

#[test]
fn modification_and_variants_compile_and_run() {
    let functionalized = functionalize(
        &StructureBuilder::new("benzene", Lattice::cubic(10.0).unwrap())
            .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
            .site(Site::new("H1", "H", [0.1, 0.0, 0.0]).unwrap())
            .build()
            .unwrap(),
        &FunctionalizeOptions {
            group: FunctionalGroup::Amino,
            ..FunctionalizeOptions::default()
        },
    )
    .unwrap();
    assert!(
        functionalized
            .sites()
            .iter()
            .any(|site| site.element() == "N")
    );
    let interpenetrated = interpenetrate(&framework(), 2).unwrap();
    assert!(interpenetrated.sites().len() >= framework().sites().len());
    let layered = stack_layers(&framework(), 2, [0.5, 0.0]).unwrap();
    assert_eq!(layered.sites().len(), 12);
    let mixed = mixed_framework(
        &framework(),
        &MixtureOptions {
            map: BTreeMap::from([("Zn".to_owned(), "Co".to_owned())]),
            ordered: true,
            fraction: 1.0,
            seed: 2,
        },
    )
    .unwrap();
    assert!(mixed.sites().iter().any(|site| site.element() == "Co"));
    let defects = xtalforge_mof::missing_defects(
        &framework(),
        &DefectOptions {
            linker_fraction: 0.0,
            cluster_fraction: 0.0,
            cap: CapSpecies::Water,
            seed: 0,
        },
    )
    .unwrap();
    assert_eq!(defects.structure.sites().len(), framework().sites().len());
}

#[test]
fn database_reader_and_screening_are_deterministic() {
    let root =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/external/core_mof");
    let (entries, info) = read_core_mof(&root).unwrap();
    assert!(!entries.is_empty());
    assert_eq!(info.entries, entries.len());
    let records = screen_database(
        &entries[..1],
        &ScreeningOptions {
            spacing: 2.0,
            samples: 4,
            limit: 1,
            ..ScreeningOptions::default()
        },
    )
    .unwrap();
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].id, entries[0].id);
    let groups =
        deduplicate(&[entries[0].structure.clone(), entries[0].structure.clone()]).unwrap();
    assert_eq!(groups.len(), 1);
    let _ = charge_balance(&framework()).unwrap();
    let _ = resolve_disorder(&framework(), &DisorderOptions::default()).unwrap();
}
