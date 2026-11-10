//! External validation against the CoRE MOF 2019 subset (Task 0.9).
//!
//! The files under `tests/external/core_mof/` are unmodified CoRE MOF 2019
//! structures (CC BY 4.0; citation in that directory's README). They exercise
//! the CIF reader, the P1 expansion, cell volume, formula and bond perception
//! on real frameworks, including charged and mixed-metal ones.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_io::{read_cif, write_cif};

macro_rules! cif {
    ($name:literal) => {
        (
            $name,
            include_str!(concat!("../../../tests/external/core_mof/", $name)),
        )
    };
}

const FILES: &[(&str, &str)] = &[
    cif!("ABAVIJ_clean_pacman.cif"),
    cif!("ABEXEM_clean_pacman.cif"),
    cif!("ACOGAB_clean_pacman.cif"),
    cif!("ADATOQ_clean_pacman.cif"),
    cif!("AGUWUV_clean_pacman.cif"),
    cif!("AHOKAJ_clean_pacman.cif"),
    cif!("ALIPEQ_clean_pacman.cif"),
    cif!("AVIPAX_clean_pacman.cif"),
    cif!("BEKYOG_clean_pacman.cif"),
    cif!("CELZEA_clean_pacman.cif"),
    cif!("ELOFIU_charged_pacman.cif"),
    cif!("FALZOI_clean_pacman.cif"),
    cif!("FIQCEN_clean_pacman.cif"),
    cif!("HAPBEI_clean_pacman.cif"),
    cif!("HUCBAJ_charged_pacman.cif"),
    cif!("IRISAD01_charged_pacman.cif"),
    cif!("OFERUN_clean_pacman.cif"),
    cif!("RUBTAK03_auto_pacman.cif"),
];

const METALS: &[(&str, &str)] = &[
    ("ABAVIJ_clean_pacman.cif", "Co"),
    ("ABEXEM_clean_pacman.cif", "La"),
    ("ACOGAB_clean_pacman.cif", "Cd"),
    ("ADATOQ_clean_pacman.cif", "Ag"),
    ("AGUWUV_clean_pacman.cif", "Ni"),
    ("AHOKAJ_clean_pacman.cif", "Fe"),
    ("ALIPEQ_clean_pacman.cif", "Zr"),
    ("AVIPAX_clean_pacman.cif", "Mg"),
    ("BEKYOG_clean_pacman.cif", "Cr"),
    ("CELZEA_clean_pacman.cif", "Mo"),
    ("ELOFIU_charged_pacman.cif", "Ca"),
    ("FALZOI_clean_pacman.cif", "Na"),
    ("FIQCEN_clean_pacman.cif", "Cu"),
    ("HAPBEI_clean_pacman.cif", "Ti"),
    ("HUCBAJ_charged_pacman.cif", "Li"),
    ("IRISAD01_charged_pacman.cif", "Bi"),
    ("OFERUN_clean_pacman.cif", "Zn"),
    ("RUBTAK03_auto_pacman.cif", "Zr"),
];

#[test]
fn every_core_mof_file_parses_and_expands() {
    assert_eq!(FILES.len(), 18);
    for (name, text) in FILES {
        let structure = read_cif(text).unwrap_or_else(|error| panic!("{name}: {error}"));
        let summary = structure
            .summary()
            .unwrap_or_else(|error| panic!("{name}: {error}"));
        assert!(summary.n_sites > 0, "{name} has no sites");
        assert!(!summary.formula.trim().is_empty(), "{name} has no formula");
        assert!(
            summary.volume > 0.0 && summary.volume.is_finite(),
            "{name} volume"
        );
        assert!(
            summary.density > 0.0 && summary.density < 30.0,
            "{name} density {}",
            summary.density
        );
        let p1 = structure
            .p1()
            .unwrap_or_else(|error| panic!("{name}: {error}"));
        assert!(p1.atoms.len() >= summary.n_sites, "{name} P1 expansion");
        let graph = structure
            .bond_graph()
            .unwrap_or_else(|error| panic!("{name}: {error}"));
        assert!(!graph.bonds.is_empty(), "{name} has no perceived bonds");
        let expected = METALS
            .iter()
            .find(|(file, _)| file == name)
            .map(|(_, metal)| *metal)
            .expect("a metal for every file");
        assert!(
            structure
                .sites()
                .iter()
                .any(|site| site.element() == expected),
            "{name} has no {expected}"
        );
    }
}

#[test]
fn cell_volume_matches_the_declared_value() {
    for (name, text) in FILES {
        let structure = read_cif(text).unwrap();
        let summary = structure.summary().unwrap();
        let declared = text.lines().find_map(|line| {
            let line = line.trim();
            line.strip_prefix("_cell_volume")
                .and_then(|rest| rest.split_whitespace().next())
                .and_then(|value| value.parse::<f64>().ok())
        });
        if let Some(declared) = declared {
            assert!(
                (summary.volume - declared).abs() < 0.05,
                "{name}: computed {} vs declared {declared}",
                summary.volume
            );
        }
    }
}

#[test]
fn formulas_match_the_cif_formula_sum() {
    let expected = [
        ("FIQCEN_clean_pacman.cif", "Cu12 H24 C72 O48"),
        ("OFERUN_clean_pacman.cif", "Zn6 H60 C48 N24"),
        ("RUBTAK03_auto_pacman.cif", "Zr24 H112 C192 O128"),
        ("IRISAD01_charged_pacman.cif", "Bi2 H8 C8 N2 O16"),
    ];
    for (name, formula) in expected {
        let text = FILES
            .iter()
            .find(|(file, _)| *file == name)
            .map(|(_, text)| *text)
            .expect("file");
        let structure = read_cif(text).unwrap();
        assert_eq!(structure.formula(), formula, "{name}");
    }
}

#[test]
fn cif_round_trip_preserves_sites_and_volume() {
    for (name, text) in FILES {
        let structure = read_cif(text).unwrap();
        let written = write_cif(&structure).unwrap();
        let again = read_cif(&written).unwrap_or_else(|error| panic!("{name}: {error}"));
        assert_eq!(
            again.sites().len(),
            structure.sites().len(),
            "{name} site count"
        );
        let first = structure.summary().unwrap().volume;
        let second = again.summary().unwrap().volume;
        assert!(
            (first - second).abs() < 1e-6 * first.max(1.0),
            "{name}: round-trip volume {first} -> {second}"
        );
        let elements: Vec<&str> = structure
            .sites()
            .iter()
            .map(|site| site.element())
            .collect();
        let again_elements: Vec<&str> = again.sites().iter().map(|site| site.element()).collect();
        assert_eq!(elements, again_elements, "{name} elements");
    }
}
