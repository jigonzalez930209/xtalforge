//! Integration tests for the format registry, atomic writes and the new
//! structure formats (Tasks 3.1-3.3).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::{Lattice, Site, StructureBuilder};
use xtalforge_io::registry::{find, losses, sniff};
use xtalforge_io::{read_cssr, read_gen, read_pmg_json, write_cssr, write_gen, write_pmg_json};

fn nacl() -> xtalforge_core::Structure {
    StructureBuilder::new("NaCl", Lattice::cubic(5.64).unwrap())
        .site(Site::new("Na1", "Na", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("Cl1", "Cl", [0.5, 0.5, 0.5]).unwrap())
        .build()
        .unwrap()
}

#[test]
fn registry_finds_names_and_extensions() {
    assert_eq!(find("cif").unwrap().name, "cif");
    assert_eq!(find(".vasp").unwrap().name, "poscar");
    assert_eq!(find("pmg-json").unwrap().name, "pmg-json");
    assert!(find("unknown").is_none());
    assert_eq!(sniff("data_x\n_cell_length_a 4").unwrap().name, "cif");
    let lost = losses(
        find("cif").unwrap().capabilities,
        find("xyz").unwrap().capabilities,
    );
    assert!(lost.contains(&"space-group symmetry"));
    assert!(lost.contains(&"partial occupancy"));
}

#[test]
fn pymatgen_json_round_trip() {
    let text = write_pmg_json(&nacl()).unwrap();
    assert!(text.contains("@class"));
    let back = read_pmg_json(&text).unwrap();
    assert_eq!(back.formula(), "Na Cl");
    assert_eq!(back.sites().len(), 2);
    assert!((back.lattice().volume() - 179.406_144).abs() < 1e-5);
}

#[test]
fn cssr_round_trip_is_p1() {
    let text = write_cssr(&nacl()).unwrap();
    assert!(text.contains("SPGR =  1 P 1"));
    let back = read_cssr(&text).unwrap();
    assert_eq!(back.sites().len(), 2);
    assert!((back.lattice().lengths()[0] - 5.64).abs() < 1e-3);
}

#[test]
fn gen_round_trip_keeps_the_cell() {
    let text = write_gen(&nacl()).unwrap();
    assert!(text.starts_with("2 F"));
    let back = read_gen(&text).unwrap();
    assert_eq!(back.sites().len(), 2);
    assert!((back.lattice().volume() - 179.406_144).abs() < 1e-5);
    assert_eq!(back.formula(), "Na Cl");
}

#[test]
fn invalid_format_text_returns_typed_errors() {
    assert!(read_pmg_json("{}").is_err());
    assert!(read_cssr("short").is_err());
    assert!(read_gen("1 C\nC\n1 1 0 0").is_err());
}

#[test]
fn atomic_write_replaces_and_leaves_no_temporary_file() {
    let directory = std::env::temp_dir().join(format!("xtalforge-atomic-{}", std::process::id()));
    std::fs::create_dir_all(&directory).unwrap();
    let path = directory.join("structure.cif");
    xtalforge_io::atomic::write_text(&path, "first").unwrap();
    xtalforge_io::atomic::write_text(&path, "second").unwrap();
    assert_eq!(std::fs::read_to_string(&path).unwrap(), "second");
    let leftovers: Vec<_> = std::fs::read_dir(&directory)
        .unwrap()
        .filter_map(|entry| entry.ok())
        .filter(|entry| entry.file_name().to_string_lossy().ends_with(".tmp"))
        .collect();
    assert!(leftovers.is_empty());
    std::fs::remove_dir_all(&directory).unwrap();
}

#[test]
fn export_policy_applies_representation_and_filters() {
    let structure = nacl();
    let p1 = xtalforge_io::registry::ExportPolicy {
        expand_to_p1: true,
        ..Default::default()
    }
    .apply(&structure)
    .unwrap();
    assert_eq!(p1.sites().len(), 2);
    let without_h = xtalforge_io::registry::ExportPolicy {
        keep_hydrogens: false,
        ..Default::default()
    };
    let filtered = without_h.apply(&p1).unwrap();
    assert_eq!(filtered.sites().len(), 2);
}
