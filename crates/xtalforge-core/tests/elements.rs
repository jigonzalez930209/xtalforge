//! Integration tests for the periodic-table data and symbol parsing.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::elements::{
    all, atomic_number, canonical_symbol, data_source, element_by_symbol, element_by_z,
    element_info, is_symbol, isotope, isotopes, mass_of, parse_symbol, schema_version,
    symbol_from_z,
};

#[test]
fn table_has_every_element_plus_the_species() {
    let table = all().unwrap();
    assert_eq!(table.len(), 120);
    assert_eq!(table[0].symbol, "H");
    assert_eq!(table[117].symbol, "Og");
    assert_eq!(table[118].symbol, "D");
    assert_eq!(table[119].symbol, "X");
    assert_eq!(schema_version().unwrap(), 1);
    assert!(data_source().unwrap().contains("Crystal-Builder"));
    assert_eq!(isotopes().unwrap().len(), 1);
}

#[test]
fn iron_reference_values() {
    let iron = element_by_symbol("Fe").unwrap();
    assert_eq!(iron.z, 26);
    assert_eq!(iron.name, "Iron");
    assert!((iron.mass - 55.845).abs() < 1e-12);
    assert!((iron.covalent_radius - 1.32).abs() < 1e-6);
    assert!((iron.vdw_radius - 2.44).abs() < 1e-12);
    assert_eq!(iron.color, [224, 102, 51]);
    assert!(iron.is_metal);
    assert_eq!(iron.valence, None);
    assert!((iron.electronegativity.unwrap() - 1.83).abs() < 1e-12);
    assert!((iron.pyykko_single.unwrap() - 1.16).abs() < 1e-12);
    assert!((iron.pyykko_double.unwrap() - 1.09).abs() < 1e-12);
    assert!((iron.pyykko_aromatic.unwrap() - 1.125).abs() < 1e-12);
}

#[test]
fn pyykko_radii_and_isotopes() {
    let carbon = element_by_symbol("C").unwrap();
    assert!((carbon.pyykko_single.unwrap() - 0.75).abs() < 1e-12);
    assert!((carbon.pyykko_double.unwrap() - 0.67).abs() < 1e-12);
    assert!((carbon.pyykko_aromatic.unwrap() - 0.71).abs() < 1e-12);

    let hydrogen = element_by_symbol("H").unwrap();
    assert_eq!(hydrogen.pyykko_double, None);
    assert_eq!(hydrogen.pyykko_aromatic, hydrogen.pyykko_single);

    let deuterium = element_by_symbol("D").unwrap();
    assert!(deuterium.is_isotope());
    assert!((deuterium.mass - 2.0141).abs() < 1e-12);
    assert_eq!(isotope("D").unwrap().unwrap().element, "H");
    assert!((mass_of("D").unwrap() - 2.0141).abs() < 1e-12);
    assert!((mass_of("d").unwrap() - 2.0141).abs() < 1e-12);
}

#[test]
fn dummy_species_is_explicit() {
    let dummy = element_by_z(0).unwrap();
    assert_eq!(dummy.symbol, "X");
    assert!(dummy.is_dummy());
    assert_eq!(dummy.z, 0);
    assert_eq!(dummy.color, [255, 20, 147]);
    assert_eq!(dummy.pyykko_single, None);
}

#[test]
fn every_atomic_number_round_trips() {
    for z in 1..=118 {
        let symbol = symbol_from_z(z).unwrap();
        let element = element_by_z(z).unwrap();
        assert_eq!(element.symbol, symbol);
        assert_eq!(atomic_number(&symbol).unwrap(), z);
        assert_eq!(element_by_symbol(&symbol).unwrap().z, z);
        assert!(is_symbol(&symbol).unwrap());
        assert_eq!(parse_symbol(&symbol).unwrap(), symbol);
        assert_eq!(parse_symbol(&symbol.to_lowercase()).unwrap(), symbol);
        assert_eq!(parse_symbol(&symbol.to_uppercase()).unwrap(), symbol);
    }
    assert!(element_by_z(119).is_err());
    assert!(symbol_from_z(0).is_err());
}

#[test]
fn file_shaped_symbols_are_normalized() {
    assert_eq!(parse_symbol("Fe").unwrap(), "Fe");
    assert_eq!(parse_symbol("FE").unwrap(), "Fe");
    assert_eq!(parse_symbol("fe").unwrap(), "Fe");
    assert_eq!(parse_symbol("Fe2+").unwrap(), "Fe");
    assert_eq!(parse_symbol("Fe1").unwrap(), "Fe");
    assert_eq!(parse_symbol("O_1").unwrap(), "O");
    assert_eq!(parse_symbol("C1A").unwrap(), "C");
    assert_eq!(parse_symbol("Ow").unwrap(), "O");
    assert_eq!(parse_symbol("Co").unwrap(), "Co");
    assert_eq!(parse_symbol("CO").unwrap(), "Co");
    assert_eq!(parse_symbol("C a").unwrap(), "C");
    assert_eq!(parse_symbol("D").unwrap(), "D");
    assert_eq!(parse_symbol("X").unwrap(), "X");
    assert_eq!(parse_symbol("Kryptonite").unwrap(), "Kr");
}

#[test]
fn strict_symbols_reject_labels() {
    assert_eq!(canonical_symbol("fe").unwrap().as_deref(), Some("Fe"));
    assert_eq!(canonical_symbol("FE").unwrap().as_deref(), Some("Fe"));
    assert_eq!(canonical_symbol("D").unwrap().as_deref(), Some("D"));
    assert_eq!(canonical_symbol("Fe2+").unwrap(), None);
    assert_eq!(canonical_symbol("C1").unwrap(), None);
    assert_eq!(canonical_symbol("Kryptonite").unwrap(), None);
    assert_eq!(canonical_symbol("").unwrap(), None);
    assert_eq!(canonical_symbol("Zz").unwrap(), None);
}

#[test]
fn invalid_symbols_return_errors() {
    for text in ["", "   ", "123", "!", "Zz"] {
        let error = element_by_symbol(text).expect_err("invalid symbol must fail");
        assert_eq!(error.code(), "invalid_input");
    }
}

#[test]
fn element_info_accepts_numbers_and_symbols() {
    assert_eq!(element_info("26").unwrap().symbol, "Fe");
    assert_eq!(element_info("fe").unwrap().symbol, "Fe");
    assert_eq!(element_info("Fe2+").unwrap().symbol, "Fe");
    assert_eq!(element_info("0").unwrap().symbol, "X");
    assert!(element_info("119").is_err());
    assert!(element_info("Zz").is_err());
}

#[test]
fn shared_element_vectors() {
    let vectors: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/vectors/elements.json")).unwrap();
    let cases = vectors["elements"].as_array().unwrap();
    assert!(!cases.is_empty());
    for case in cases {
        let input = case["input"].as_str().unwrap();
        let element = element_info(input).unwrap();
        assert_eq!(element.symbol, case["symbol"].as_str().unwrap());
        assert_eq!(element.z, case["z"].as_u64().unwrap() as u8);
        let tolerance = case["tolerance"].as_f64().unwrap();
        for (field, value) in [
            ("mass", element.mass),
            ("covalent_radius", element.covalent_radius),
            ("vdw_radius", element.vdw_radius),
        ] {
            let expected = case[field].as_f64().unwrap();
            assert!(
                (value - expected).abs() <= tolerance,
                "{input}: {field} {value} != {expected}"
            );
        }
    }
}

#[test]
fn element_is_serializable() {
    let iron = element_by_symbol("Fe").unwrap();
    let json = serde_json::to_string(iron).unwrap();
    let back: xtalforge_core::Element = serde_json::from_str(&json).unwrap();
    assert_eq!(&back, iron);
}
