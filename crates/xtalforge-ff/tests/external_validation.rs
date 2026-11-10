//! External force-field validation against sci-form (Task 4.2.5).
//!
//! sci-form is a separate MIT implementation of UFF by the same maintainer.
//! The same Cartesian coordinates are evaluated by both; the energies must
//! agree to numerical precision, which catches a wrong parameter, a missing
//! term or a different van der Waals combination.

#![cfg(feature = "sci-form")]
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::Structure;
use xtalforge_ff::uff::params;
use xtalforge_ff::{CalculatorOptions, UffModel};

/// One eV in kcal/mol, the unit sci-form reports.
const EV_TO_KCAL: f64 = 23.060_548;

fn molecule(name: &str, atoms: &[(&str, [f64; 3])]) -> Structure {
    let points: Vec<(String, [f64; 3])> = atoms
        .iter()
        .map(|(element, point)| ((*element).to_owned(), *point))
        .collect();
    Structure::molecule(name, &points, 8.0).expect("molecule")
}

fn coordinates(structure: &Structure) -> Vec<f64> {
    structure
        .sites()
        .iter()
        .flat_map(|site| {
            let cart = structure.lattice().to_cart(site.frac());
            [cart[0], cart[1], cart[2]]
        })
        .collect()
}

fn energy_kcal(structure: &Structure) -> f64 {
    let model = UffModel::build(structure, CalculatorOptions::default()).expect("UFF model");
    model.energy().expect("UFF energy") * EV_TO_KCAL
}

fn tetrahedral_methane(scale: f64) -> Structure {
    let length = xtalforge_ff::uff::natural_bond_length("C_3", "H_", 1.0).expect("C-H");
    let a = length / 3.0_f64.sqrt();
    molecule(
        "methane",
        &[
            ("C", [0.0, 0.0, 0.0]),
            ("H", [a * scale, a, a]),
            ("H", [a, -a, -a]),
            ("H", [-a, a, -a]),
            ("H", [-a, -a, a]),
        ],
    )
}

#[test]
fn uff_matches_sci_form_on_methane() {
    let methane = tetrahedral_methane(1.0);
    let ours = energy_kcal(&methane);
    let theirs = sci_form::compute_uff_energy("C", &coordinates(&methane)).expect("sci-form");
    assert!(
        (ours - theirs).abs() < 1e-5,
        "methane: ours {ours}, sci-form {theirs}"
    );
}

#[test]
fn uff_matches_sci_form_on_a_distorted_methane() {
    let methane = tetrahedral_methane(1.10);
    let ours = energy_kcal(&methane);
    let theirs = sci_form::compute_uff_energy("C", &coordinates(&methane)).expect("sci-form");
    assert!(
        (ours - theirs).abs() < 1e-4,
        "distorted methane: ours {ours}, sci-form {theirs}"
    );
    assert!(ours > 0.0, "a stretched bond must cost energy");
}

#[test]
fn uff_matches_sci_form_on_stretched_water() {
    let length = xtalforge_ff::uff::natural_bond_length("O_3", "H_", 1.0).expect("O-H");
    let angle = params::get("O_3").expect("O_3").theta0.to_radians();
    let water = molecule(
        "water",
        &[
            ("O", [0.0, 0.0, 0.0]),
            ("H", [length * 1.05, 0.0, 0.0]),
            ("H", [length * angle.cos(), length * angle.sin(), 0.0]),
        ],
    );
    let ours = energy_kcal(&water);
    let theirs = sci_form::compute_uff_energy("O", &coordinates(&water)).expect("sci-form");
    assert!(
        (ours - theirs).abs() < 1e-4,
        "water: ours {ours}, sci-form {theirs}"
    );
}
