//! Cell properties are counted on the whole cell, not on the asymmetric unit.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_io::read_cif;

const NACL: &str = include_str!("../../../tests/corpus/simple/NaCl.cif");
const SILICON: &str = include_str!("../../../tests/corpus/simple/Si.cif");

#[test]
fn rock_salt_has_four_formula_units_and_its_real_density() {
    let structure = read_cif(NACL).unwrap();
    assert_eq!(structure.sites().len(), 2, "two asymmetric sites");
    assert_eq!(structure.formula(), "Na4 Cl4");
    assert_eq!(structure.z(), 4);
    // 4 NaCl (58.443 g/mol) in a = 5.62 Angstrom: 2.19 g/cm3 (halite measures
    // 2.165 at its slightly larger room-temperature cell).
    let density = structure.density().unwrap();
    assert!((density - 2.19).abs() < 0.03, "density {density}");
    assert!((structure.cell_mass().unwrap() - 4.0 * 58.44).abs() < 0.1);
    assert!(structure.total_charge().abs() < 1e-9);
    let summary = structure.summary().unwrap();
    assert!((summary.density - density).abs() < 1e-12);
    assert_eq!(summary.formula, "Na4 Cl4");
}

#[test]
fn diamond_silicon_has_eight_atoms() {
    let structure = read_cif(SILICON).unwrap();
    assert_eq!(structure.formula(), "Si8");
    assert_eq!(structure.z(), 8);
    // 2.329 g/cm3 for silicon.
    let density = structure.density().unwrap();
    assert!((density - 2.33).abs() < 0.02, "density {density}");
}
