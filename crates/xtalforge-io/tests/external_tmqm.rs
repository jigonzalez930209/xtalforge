//! External validation against the tmQM subset (Task 0.9).
//!
//! The 18 files under `tests/external/tmqm/` are unmodified tmQM records (MIT,
//! see that directory's README): one transition-metal complex per metal across
//! the 3d, 4d and 5d series. The header states the atom count and the
//! stoichiometry, which makes it an external check of the XYZ reader.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::collections::BTreeMap;

use xtalforge_io::{read_xyz, write_xyz};

macro_rules! xyz {
    ($name:literal) => {
        (
            $name,
            include_str!(concat!("../../../tests/external/tmqm/", $name)),
        )
    };
}

const FILES: &[(&str, &str)] = &[
    xyz!("ANAXOE_tmQM.xyz"),
    xyz!("CEYJEZ_tmQM.xyz"),
    xyz!("CILXOO_tmQM.xyz"),
    xyz!("CPMINC_tmQM.xyz"),
    xyz!("DAQHAI_tmQM.xyz"),
    xyz!("DASFIQ_tmQM.xyz"),
    xyz!("HEWTEM_tmQM.xyz"),
    xyz!("IPOSAK_tmQM.xyz"),
    xyz!("LAXLAB_tmQM.xyz"),
    xyz!("PITRAM_tmQM.xyz"),
    xyz!("QOSJIR_tmQM.xyz"),
    xyz!("REJFUL_tmQM.xyz"),
    xyz!("UMAHUO_tmQM.xyz"),
    xyz!("UTANIP_tmQM.xyz"),
    xyz!("VUDZEC_tmQM.xyz"),
    xyz!("WUGVEZ_tmQM.xyz"),
    xyz!("XARWUJ_tmQM.xyz"),
    xyz!("ZOYPEL_tmQM.xyz"),
];

fn header(text: &str) -> (&str, &str) {
    let mut lines = text.lines();
    let count = lines.next().expect("atom count");
    let comment = lines.next().expect("comment");
    (count, comment)
}

fn stoichiometry(comment: &str) -> &str {
    comment
        .split('|')
        .find_map(|part| part.trim().strip_prefix("Stoichiometry = "))
        .expect("stoichiometry in the header")
}

fn element_counts(formula: &str) -> BTreeMap<String, usize> {
    let main = formula.split('(').next().unwrap_or(formula);
    let mut counts = BTreeMap::new();
    let mut chars = main.chars().peekable();
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
        let count = digits.parse::<usize>().unwrap_or(1);
        *counts.entry(element).or_insert(0) += count;
    }
    counts
}

#[test]
fn every_tmqm_record_parses_with_its_header_counts() {
    assert_eq!(FILES.len(), 18);
    for (name, text) in FILES {
        let (count, comment) = header(text);
        let declared: usize = count.parse().expect("atom count");
        let structure = read_xyz(text).unwrap_or_else(|error| panic!("{name}: {error}"));
        assert_eq!(structure.sites().len(), declared, "{name} atom count");

        let expected = element_counts(stoichiometry(comment));
        let mut actual: BTreeMap<String, usize> = BTreeMap::new();
        for site in structure.sites() {
            *actual.entry(site.element().to_owned()).or_insert(0) += 1;
        }
        assert_eq!(actual, expected, "{name} stoichiometry");

        let graph = structure
            .bond_graph()
            .unwrap_or_else(|error| panic!("{name}: {error}"));
        let metal = structure.sites().first().expect("first atom").element();
        assert!(
            graph.degree(0) > 0,
            "{name}: the {metal} atom has no perceived bond"
        );
    }
}

#[test]
fn xyz_round_trip_preserves_the_atoms() {
    for (name, text) in FILES {
        let structure = read_xyz(text).unwrap();
        let written = write_xyz(&structure).unwrap();
        let again = read_xyz(&written).unwrap_or_else(|error| panic!("{name}: {error}"));
        assert_eq!(again.sites().len(), structure.sites().len(), "{name}");
        let first: Vec<&str> = structure
            .sites()
            .iter()
            .map(|site| site.element())
            .collect();
        let second: Vec<&str> = again.sites().iter().map(|site| site.element()).collect();
        assert_eq!(first, second, "{name}");
    }
}
