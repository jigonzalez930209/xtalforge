//! External validation against the QMOF subset (Task 0.9).
//!
//! Fourteen DFT-relaxed MOF structures from the QMOF Database (CC BY 4.0; see
//! `tests/external/qmof/README.md`). The tabulated QMOF atom count, volume and
//! density are external checks of the CIF reader; four structures also check
//! that symmetry detection recovers the QMOF space group from the P1 file.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::collections::BTreeMap;

use xtalforge_io::{read_cif, write_cif};

macro_rules! cif {
    ($name:literal) => {
        (
            $name,
            include_str!(concat!("../../../tests/external/qmof/", $name)),
        )
    };
}

/// File, QMOF P1 atom count, volume (A^3) and density (g/cm3).
const FILES: &[(&str, usize, f64, f64)] = &[
    ("qmof-0331d6b.cif", 52, 634.4377163, 2.183926517),
    ("qmof-36ade61.cif", 26, 330.7816054, 3.065290240),
    ("qmof-40f119f.cif", 17, 210.8922878, 2.045622412),
    ("qmof-421f50d.cif", 26, 283.3960503, 2.141538839),
    ("qmof-5d3a22f.cif", 20, 193.3513918, 3.090675263),
    ("qmof-7a9d846.cif", 29, 270.8557715, 2.182629188),
    ("qmof-7aa58f3.cif", 52, 573.2955044, 2.044036911),
    ("qmof-832a1ad.cif", 17, 220.9817634, 1.887541028),
    ("qmof-89d2f92.cif", 23, 326.1823619, 2.213976994),
    ("qmof-9247cfe.cif", 58, 635.3724487, 2.322573422),
    ("qmof-b6b9d01.cif", 20, 390.4659742, 2.003314030),
    ("qmof-c3fa563.cif", 24, 236.5909663, 1.853758960),
    ("qmof-e1456f8.cif", 28, 338.0063212, 2.181555382),
    ("qmof-ea5d153.cif", 26, 294.2946500, 2.031370271),
];

const TEXT: &[(&str, &str)] = &[
    cif!("qmof-0331d6b.cif"),
    cif!("qmof-36ade61.cif"),
    cif!("qmof-40f119f.cif"),
    cif!("qmof-421f50d.cif"),
    cif!("qmof-5d3a22f.cif"),
    cif!("qmof-7a9d846.cif"),
    cif!("qmof-7aa58f3.cif"),
    cif!("qmof-832a1ad.cif"),
    cif!("qmof-89d2f92.cif"),
    cif!("qmof-9247cfe.cif"),
    cif!("qmof-b6b9d01.cif"),
    cif!("qmof-c3fa563.cif"),
    cif!("qmof-e1456f8.cif"),
    cif!("qmof-ea5d153.cif"),
];

fn text(name: &str) -> &'static str {
    TEXT.iter()
        .find(|(file, _)| *file == name)
        .map(|(_, text)| *text)
        .expect("fixture")
}

fn cif_formula_counts(text: &str) -> BTreeMap<String, usize> {
    let formula = text
        .lines()
        .find_map(|line| line.trim().strip_prefix("_chemical_formula_sum"))
        .map(|rest| rest.trim().trim_matches('\'').to_owned())
        .expect("formula sum");
    let mut counts = BTreeMap::new();
    let mut chars = formula.chars().peekable();
    while let Some(letter) = chars.next() {
        if !letter.is_ascii_uppercase() {
            continue;
        }
        let mut element = letter.to_string();
        if let Some(lower) = chars.peek()
            && lower.is_ascii_lowercase()
        {
            element.push(chars.next().unwrap_or_default());
        }
        let mut digits = String::new();
        while let Some(digit) = chars.peek()
            && digit.is_ascii_digit()
        {
            digits.push(chars.next().unwrap_or_default());
        }
        *counts.entry(element).or_insert(0) += digits.parse::<usize>().unwrap_or(1);
    }
    counts
}

#[test]
fn every_qmof_file_matches_the_tabulated_values() {
    assert_eq!(FILES.len(), 14);
    for (name, atoms, volume, density) in FILES {
        let structure = read_cif(text(name)).unwrap_or_else(|error| panic!("{name}: {error}"));
        let summary = structure.summary().unwrap();
        let p1 = structure.p1().unwrap();
        assert_eq!(p1.atoms.len(), *atoms, "{name} atom count");
        assert!(
            (summary.volume - volume).abs() < 0.05,
            "{name}: volume {} vs {volume}",
            summary.volume
        );
        assert!(
            (summary.density - density).abs() < 0.002,
            "{name}: density {} vs {density}",
            summary.density
        );
        let mut actual: BTreeMap<String, usize> = BTreeMap::new();
        for site in structure.sites() {
            *actual.entry(site.element().to_owned()).or_insert(0) += 1;
        }
        assert_eq!(actual, cif_formula_counts(text(name)), "{name} formula");
    }
}

#[test]
fn qmof_cif_round_trip_preserves_the_cell() {
    for (name, _, volume, _) in FILES {
        let structure = read_cif(text(name)).unwrap();
        let written = write_cif(&structure).unwrap();
        let again = read_cif(&written).unwrap_or_else(|error| panic!("{name}: {error}"));
        assert_eq!(again.sites().len(), structure.sites().len(), "{name}");
        assert!(
            (again.summary().unwrap().volume - volume).abs() < 0.05,
            "{name}"
        );
    }
}

#[test]
fn symmetry_detection_recovers_the_qmof_space_group() {
    // Detection tolerances are in Angstrom; the relaxed cells are slightly
    // distorted, so 0.01 A is the value that recovers the tabulated group.
    let expected = [
        ("qmof-ea5d153.cif", 167_u32),
        ("qmof-832a1ad.cif", 12),
        ("qmof-40f119f.cif", 12),
        ("qmof-36ade61.cif", 51),
    ];
    for (name, number) in expected {
        let structure = read_cif(text(name)).unwrap();
        let dataset = xtalforge_symmetry::detect_symmetry(
            &structure,
            xtalforge_symmetry::SymmetryOptions {
                symprec: 0.01,
                ..Default::default()
            },
        )
        .unwrap_or_else(|error| panic!("{name}: {error}"));
        assert_eq!(dataset.number, number, "{name} space group");
    }
}
