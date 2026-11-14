//! Integration tests for phonon DOS and inelastic-neutron weights (Task 19.4.3).
//!
//! The weights are checked against the closed form
//! `sum_i (b_i^2 / m_i) |u_i|^2 / sum_i |u_i|^2` for a translation and for an
//! opposite-displacement optical mode.  The DOS checks use the unit-area
//! normalisation of the histogram.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::collections::BTreeMap;

use xtalforge_core::elements::mass_of;
use xtalforge_core::{Lattice, Site, Structure, StructureBuilder};
use xtalforge_ff::{
    DosOptions, Error, NeutronWeights, PhononResult, VibrationalModes, inelastic_neutron_dos,
    mode_neutron_weights, phonon_dos,
};

const CARBON_LENGTH: f64 = 6.646;
const HYDROGEN_LENGTH: f64 = -3.7406;

fn carbon_hydrogen() -> Structure {
    StructureBuilder::new("carbon-hydrogen", Lattice::cubic(10.0).unwrap())
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("H1", "H", [0.25, 0.0, 0.0]).unwrap())
        .build()
        .unwrap()
}

fn vibration_modes() -> VibrationalModes {
    let units = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let mut vectors = units
        .iter()
        .map(|unit| vec![*unit, *unit])
        .collect::<Vec<_>>();
    vectors.extend(units.iter().map(|unit| {
        vec![
            [0.5 * unit[0], 0.5 * unit[1], 0.5 * unit[2]],
            [-unit[0], -unit[1], -unit[2]],
        ]
    }));
    VibrationalModes {
        frequencies: vec![0.0, 0.0, 0.0, 300.0, 400.0, 500.0],
        modes: vectors,
        zero_modes: vec![0, 1, 2],
        imaginary: Vec::new(),
    }
}

fn phonon_result() -> PhononResult {
    PhononResult {
        q_points: vec![[0.0; 3], [0.5, 0.0, 0.0]],
        frequencies: vec![
            vec![0.0, 0.0, 0.0, 300.0, 400.0, 500.0],
            vec![10.0, 10.0, 10.0, 320.0, 420.0, 520.0],
        ],
        modes: vibration_modes(),
        dos: None,
        note: "integration fixture".to_owned(),
    }
}

fn weights() -> NeutronWeights {
    NeutronWeights::from_lengths(vec![CARBON_LENGTH, HYDROGEN_LENGTH]).unwrap()
}

fn dos_options() -> DosOptions {
    DosOptions {
        bins: 120,
        sigma: 8.0,
        min: Some(0.0),
        max: Some(600.0),
    }
}

fn integral(dos: &xtalforge_ff::PhononDos) -> f64 {
    let step = dos.frequencies[1] - dos.frequencies[0];
    dos.density.iter().sum::<f64>() * step
}

fn bin_near(dos: &xtalforge_ff::PhononDos, frequency: f64) -> f64 {
    let index = dos
        .frequencies
        .iter()
        .enumerate()
        .min_by(|(_, first), (_, second)| {
            (*first - frequency)
                .abs()
                .partial_cmp(&(*second - frequency).abs())
                .unwrap()
        })
        .map(|(index, _)| index)
        .unwrap();
    dos.density[index]
}

#[test]
fn neutron_weights_follow_scattering_length_squared_over_mass() {
    let structure = carbon_hydrogen();
    let modes = vibration_modes();
    let found = mode_neutron_weights(&structure, &modes, &weights()).unwrap();
    assert_eq!(found.len(), 6);

    let mass_c = mass_of("C").unwrap();
    let mass_h = mass_of("H").unwrap();
    let translation =
        (CARBON_LENGTH * CARBON_LENGTH / mass_c + HYDROGEN_LENGTH * HYDROGEN_LENGTH / mass_h) / 2.0;
    // The optical modes move carbon by 0.5 and hydrogen by 1.0, so the
    // hydrogen contribution enters with a factor of four.
    let optical = (0.25 * CARBON_LENGTH * CARBON_LENGTH / mass_c
        + HYDROGEN_LENGTH * HYDROGEN_LENGTH / mass_h)
        / 1.25;
    for (mode, weight) in found.iter().take(3).enumerate() {
        assert!(
            (weight - translation).abs() < 1e-12 * translation,
            "translation {mode}: {weight} != {translation}"
        );
    }
    for (mode, weight) in found.iter().skip(3).enumerate() {
        assert!(
            (weight - optical).abs() < 1e-12 * optical,
            "optical {mode}: {weight} != {optical}"
        );
    }
    assert!(
        optical > translation,
        "hydrogen should dominate the optical modes"
    );
}

#[test]
fn phonon_dos_is_normalised_to_unit_area() {
    let result = phonon_result();
    let dos = phonon_dos(&result, &dos_options()).unwrap();
    assert_eq!(dos.frequencies.len(), 120);
    assert!((integral(&dos) - 1.0).abs() < 1e-9, "{}", integral(&dos));
    assert!(dos.density.iter().all(|value| *value >= 0.0));
    assert_eq!(dos.sigma, 8.0);
}

#[test]
fn inelastic_neutron_dos_exposes_both_curves_and_the_branch_weights() {
    let structure = carbon_hydrogen();
    let result = phonon_result();
    let weighted = inelastic_neutron_dos(&structure, &result, &weights(), &dos_options()).unwrap();
    assert_eq!(
        weighted.dos.frequencies.len(),
        weighted.weighted.frequencies.len()
    );
    assert!((integral(&weighted.dos) - 1.0).abs() < 1e-9);
    assert!((integral(&weighted.weighted) - 1.0).abs() < 1e-9);
    assert_eq!(weighted.mode_weights.len(), 6);
    let total = weighted.mode_weights.iter().sum::<f64>();
    assert!((weighted.total_weight - total).abs() < 1e-15);
    assert!(total > 0.0);
    assert!(
        weighted
            .dos
            .density
            .iter()
            .zip(&weighted.weighted.density)
            .any(|(plain, heavy)| (plain - heavy).abs() > 1e-6),
        "the weighted curve should differ from the plain one"
    );
    let plain_peak = bin_near(&weighted.dos, 500.0);
    let weighted_peak = bin_near(&weighted.weighted, 500.0);
    assert!(
        weighted_peak > 1.1 * plain_peak,
        "hydrogen modes should dominate: {weighted_peak} vs {plain_peak}"
    );
}

#[test]
fn neutron_lengths_can_come_from_a_table_or_a_callback() {
    let structure = carbon_hydrogen();
    let mut table = BTreeMap::new();
    table.insert("C".to_owned(), CARBON_LENGTH);
    table.insert("H".to_owned(), HYDROGEN_LENGTH);
    let from_table = NeutronWeights::from_table(&structure, &table).unwrap();
    assert_eq!(from_table, weights());
    let from_callback = NeutronWeights::from_callback(&structure, |symbol| match symbol {
        "C" => Ok(CARBON_LENGTH),
        "H" => Ok(HYDROGEN_LENGTH),
        other => Err(Error::InvalidInput(format!("no length for {other}"))),
    })
    .unwrap();
    assert_eq!(from_callback, weights());
    assert_eq!(from_callback.len(), 2);
    assert!(!from_callback.is_empty());

    let mut incomplete = BTreeMap::new();
    incomplete.insert("C".to_owned(), CARBON_LENGTH);
    let error = NeutronWeights::from_table(&structure, &incomplete).unwrap_err();
    assert_eq!(error.code(), "invalid_input");
}

#[test]
fn invalid_phonon_inputs_are_rejected() {
    let structure = carbon_hydrogen();

    let mut empty = phonon_result();
    empty.frequencies.clear();
    let error = phonon_dos(&empty, &dos_options()).unwrap_err();
    assert_eq!(error.code(), "invalid_input");

    let error = mode_neutron_weights(
        &structure,
        &vibration_modes(),
        &NeutronWeights::from_lengths(vec![1.0]).unwrap(),
    )
    .unwrap_err();
    assert_eq!(error.code(), "invalid_input");

    let error = NeutronWeights::from_lengths(vec![f64::NAN, 1.0]).unwrap_err();
    assert_eq!(error.code(), "invalid_input");

    let mut mismatched = phonon_result();
    mismatched.frequencies[1].pop();
    let error =
        inelastic_neutron_dos(&structure, &mismatched, &weights(), &dos_options()).unwrap_err();
    assert_eq!(error.code(), "invalid_input");

    let zero_lengths = NeutronWeights::from_lengths(vec![0.0, 0.0]).unwrap();
    let error = inelastic_neutron_dos(&structure, &phonon_result(), &zero_lengths, &dos_options())
        .unwrap_err();
    assert_eq!(error.code(), "invalid_input");

    let mut malformed = vibration_modes();
    malformed.frequencies.pop();
    let error = mode_neutron_weights(&structure, &malformed, &weights()).unwrap_err();
    assert_eq!(error.code(), "invalid_input");

    let error = phonon_dos(
        &phonon_result(),
        &DosOptions {
            bins: 0,
            ..DosOptions::default()
        },
    )
    .unwrap_err();
    assert_eq!(error.code(), "invalid_input");
}
