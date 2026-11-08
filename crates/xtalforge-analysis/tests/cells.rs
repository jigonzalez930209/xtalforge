//! Integration tests for free cell parameters and extinction symbols
//! (Tasks 7.5.3/7.5.4).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_analysis::{
    TrigonalSetting, extinction_symbol, free_cell_parameters, free_cell_parameters_with_setting,
};

#[test]
fn each_crystal_system_reports_its_free_numbers() {
    assert_eq!(
        free_cell_parameters(1).unwrap(),
        vec!["a", "b", "c", "alpha", "beta", "gamma"]
    );
    assert_eq!(
        free_cell_parameters(14).unwrap(),
        vec!["a", "b", "c", "beta"]
    );
    assert_eq!(free_cell_parameters(62).unwrap(), vec!["a", "b", "c"]);
    assert_eq!(free_cell_parameters(139).unwrap(), vec!["a", "c"]);
    assert_eq!(free_cell_parameters(143).unwrap(), vec!["a", "c"]);
    assert_eq!(free_cell_parameters(146).unwrap(), vec!["a", "alpha"]);
    assert_eq!(free_cell_parameters(168).unwrap(), vec!["a", "c"]);
    assert_eq!(free_cell_parameters(194).unwrap(), vec!["a", "c"]);
    assert_eq!(free_cell_parameters(227).unwrap(), vec!["a"]);
    assert_eq!(free_cell_parameters(230).unwrap(), vec!["a"]);
}

#[test]
fn trigonal_hexagonal_setting_is_available_explicitly() {
    assert_eq!(
        free_cell_parameters_with_setting(146, TrigonalSetting::HexagonalAxes).unwrap(),
        vec!["a", "c"]
    );
    assert_eq!(
        free_cell_parameters_with_setting(146, TrigonalSetting::RhombohedralAxes).unwrap(),
        vec!["a", "alpha"]
    );
    assert!(free_cell_parameters_with_setting(143, TrigonalSetting::RhombohedralAxes).is_err());
}

#[test]
fn extinction_symbols_follow_the_standard_settings() {
    assert_eq!(extinction_symbol(1).unwrap(), "P");
    assert_eq!(extinction_symbol(5).unwrap(), "C");
    assert_eq!(extinction_symbol(14).unwrap(), "P");
    assert_eq!(extinction_symbol(20).unwrap(), "C");
    assert_eq!(extinction_symbol(22).unwrap(), "F");
    assert_eq!(extinction_symbol(23).unwrap(), "I");
    assert_eq!(extinction_symbol(38).unwrap(), "A");
    assert_eq!(extinction_symbol(62).unwrap(), "P");
    assert_eq!(extinction_symbol(69).unwrap(), "F");
    assert_eq!(extinction_symbol(79).unwrap(), "I");
    assert_eq!(extinction_symbol(99).unwrap(), "P");
    assert_eq!(extinction_symbol(146).unwrap(), "R");
    assert_eq!(extinction_symbol(166).unwrap(), "R");
    assert_eq!(extinction_symbol(168).unwrap(), "P");
    assert_eq!(extinction_symbol(196).unwrap(), "F");
    assert_eq!(extinction_symbol(197).unwrap(), "I");
    assert_eq!(extinction_symbol(221).unwrap(), "P");
    assert_eq!(extinction_symbol(225).unwrap(), "F");
    assert_eq!(extinction_symbol(229).unwrap(), "I");
    assert_eq!(extinction_symbol(230).unwrap(), "I");
}

#[test]
fn invalid_space_group_numbers_are_rejected() {
    assert!(free_cell_parameters(0).is_err());
    assert!(free_cell_parameters(231).is_err());
    assert!(extinction_symbol(0).is_err());
    assert!(extinction_symbol(231).is_err());
}
