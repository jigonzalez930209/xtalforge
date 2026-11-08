//! Integration tests for charges, EQeq and Ewald summation (Task 4.3).
//!
//! The QEq checks use the sign and neutrality a chemical reader expects; the
//! Ewald check compares against a spherical direct lattice sum, which is the
//! same boundary convention as Ewald and converges to it (the cubic direct sum
//! does not, because of its surface term).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::{Lattice, Site, Structure, StructureBuilder};
use xtalforge_ff::charges::{COULOMB_EV_ANGSTROM, ChargeOptions, EwaldOptions};
use xtalforge_ff::{eqeq_charges, ewald_energy_forces, qeq_charges};

fn molecule(name: &str, atoms: &[(&str, [f64; 3])]) -> Structure {
    let points: Vec<(String, [f64; 3])> = atoms
        .iter()
        .map(|(element, point)| ((*element).to_owned(), *point))
        .collect();
    Structure::molecule(name, &points, 8.0).expect("molecule")
}

fn water() -> Structure {
    molecule(
        "water",
        &[
            ("O", [0.0, 0.0, 0.0]),
            ("H", [0.96, 0.0, 0.0]),
            ("H", [-0.24, 0.93, 0.0]),
        ],
    )
}

fn nacl_molecule() -> Structure {
    molecule("nacl", &[("Na", [0.0, 0.0, 0.0]), ("Cl", [2.36, 0.0, 0.0])])
}

fn rock_salt() -> Structure {
    let lattice = Lattice::cubic(4.0).unwrap();
    let frac = [
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.5],
        [0.5, 0.5, 0.5],
        [0.0, 0.0, 0.5],
        [0.0, 0.5, 0.0],
        [0.5, 0.0, 0.0],
    ];
    let mut builder = StructureBuilder::new("rock-salt", lattice);
    for (index, point) in frac.iter().enumerate() {
        let element = if index < 4 { "Na" } else { "Cl" };
        builder = builder.site(Site::new(format!("{element}{index}"), element, *point).unwrap());
    }
    builder.build().unwrap()
}

fn rock_salt_charges() -> Vec<f64> {
    let mut charges = vec![1.0; 4];
    charges.extend(vec![-1.0; 4]);
    charges
}

fn direct_coulomb(structure: &Structure, charges: &[f64], radius: i32) -> f64 {
    let lattice = structure.lattice();
    let cart = structure.cart_sites();
    let matrix = lattice.matrix();
    let mut energy = 0.0;
    for na in -radius..=radius {
        for nb in -radius..=radius {
            for nc in -radius..=radius {
                if na * na + nb * nb + nc * nc > radius * radius {
                    continue;
                }
                let shift = [
                    na as f64 * matrix[0][0] + nb as f64 * matrix[1][0] + nc as f64 * matrix[2][0],
                    na as f64 * matrix[0][1] + nb as f64 * matrix[1][1] + nc as f64 * matrix[2][1],
                    na as f64 * matrix[0][2] + nb as f64 * matrix[1][2] + nc as f64 * matrix[2][2],
                ];
                for (i, first) in cart.iter().enumerate() {
                    for (j, second) in cart.iter().enumerate() {
                        if i == j && na == 0 && nb == 0 && nc == 0 {
                            continue;
                        }
                        let delta = [
                            second[0] + shift[0] - first[0],
                            second[1] + shift[1] - first[1],
                            second[2] + shift[2] - first[2],
                        ];
                        let distance =
                            (delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2])
                                .sqrt();
                        energy += 0.5 * COULOMB_EV_ANGSTROM * charges[i] * charges[j] / distance;
                    }
                }
            }
        }
    }
    energy
}

fn shifted(structure: &Structure, index: usize, cart_shift: [f64; 3]) -> Structure {
    let lattice = structure.lattice();
    let mut frac: Vec<[f64; 3]> = structure.sites().iter().map(Site::frac).collect();
    let delta = lattice.to_frac(cart_shift);
    if let Some(point) = frac.get_mut(index) {
        *point = [
            point[0] + delta[0],
            point[1] + delta[1],
            point[2] + delta[2],
        ];
    }
    let sites: Vec<Site> = structure
        .sites()
        .iter()
        .zip(&frac)
        .map(|(site, point)| site.clone().with_frac(*point).unwrap())
        .collect();
    let mut builder = StructureBuilder::new(structure.name(), *lattice).hall(structure.hall());
    for site in sites {
        builder = builder.site(site);
    }
    builder.build().unwrap()
}

#[test]
fn qeq_water_is_neutral_and_polar() {
    let charges = qeq_charges(&water(), &ChargeOptions::default()).unwrap();
    assert_eq!(charges.per_atom.len(), 3);
    assert!(charges.total.abs() < 1e-10);
    assert!(charges.per_atom[0] < -0.2, "oxygen should be negative");
    assert!(charges.per_atom[1] > 0.1 && charges.per_atom[2] > 0.1);
    assert!(charges.iterations >= 1);
}

#[test]
fn qeq_nacl_gives_ions() {
    let charges = qeq_charges(&nacl_molecule(), &ChargeOptions::default()).unwrap();
    assert!(charges.total.abs() < 1e-10);
    assert!(charges.per_atom[0] > 0.5, "sodium should be positive");
    assert!(charges.per_atom[1] < -0.5, "chloride should be negative");
}

#[test]
fn eqeq_nacl_gives_ions_with_the_reference_dielectric() {
    let charges = eqeq_charges(&nacl_molecule(), &ChargeOptions::eqeq()).unwrap();
    assert!(charges.total.abs() < 1e-9);
    assert!(charges.per_atom[0] > 0.5);
    assert!(charges.per_atom[1] < -0.5);
}

#[test]
fn charges_are_stored_per_site() {
    let structure = water();
    let charges = qeq_charges(&structure, &ChargeOptions::default()).unwrap();
    let per_site = charges.per_site(&structure).unwrap();
    assert_eq!(per_site.len(), structure.sites().len());
    let with_charges = charges.apply_to(&structure).unwrap();
    assert_eq!(with_charges.sites().len(), structure.sites().len());
    assert_eq!(with_charges.name(), structure.name());
    let total: f64 = with_charges
        .sites()
        .iter()
        .map(|site| site.charge().unwrap_or(0.0))
        .sum();
    assert!(total.abs() < 1e-9);
    for site in with_charges.sites() {
        assert!(site.charge().is_some());
    }
}

#[test]
fn ewald_matches_the_spherical_direct_sum() {
    let structure = rock_salt();
    let charges = rock_salt_charges();
    let options = EwaldOptions {
        accuracy: 1e-8,
        ..EwaldOptions::default()
    };
    let result = ewald_energy_forces(&structure, &charges, &options).unwrap();
    let direct = direct_coulomb(&structure, &charges, 10);
    assert!(
        (result.energy - direct).abs() < 1e-3,
        "Ewald {} vs direct {}",
        result.energy,
        direct
    );
    assert_eq!(result.forces.len(), structure.sites().len());
    let net: [f64; 3] = result.forces.iter().fold([0.0; 3], |sum, force| {
        [sum[0] + force[0], sum[1] + force[1], sum[2] + force[2]]
    });
    assert!(net[0].abs() < 1e-8 && net[1].abs() < 1e-8 && net[2].abs() < 1e-8);
}

#[test]
fn ewald_forces_match_finite_differences() {
    let structure = shifted(&rock_salt(), 0, [0.03, 0.02, -0.01]);
    let charges = rock_salt_charges();
    let options = EwaldOptions {
        accuracy: 1e-9,
        ..EwaldOptions::default()
    };
    let result = ewald_energy_forces(&structure, &charges, &options).unwrap();
    let step = 1e-5;
    for (index, force) in result.forces.iter().enumerate() {
        for axis in 0..3 {
            let mut shift = [0.0; 3];
            shift[axis] = step;
            let high = shifted(&structure, index, shift);
            let low = shifted(&structure, index, [-shift[0], -shift[1], -shift[2]]);
            let energy_high = ewald_energy_forces(&high, &charges, &options)
                .unwrap()
                .energy;
            let energy_low = ewald_energy_forces(&low, &charges, &options)
                .unwrap()
                .energy;
            let numeric = -(energy_high - energy_low) / (2.0 * step);
            assert!(
                (force[axis] - numeric).abs() < 2e-3,
                "atom {index} axis {axis}: analytic {} numeric {}",
                force[axis],
                numeric
            );
        }
    }
}

#[test]
fn ewald_rejects_bad_counts_and_options() {
    let structure = rock_salt();
    let charges = rock_salt_charges();
    assert!(ewald_energy_forces(&structure, &charges[..3], &EwaldOptions::default()).is_err());
    let bad = EwaldOptions {
        accuracy: 0.0,
        ..EwaldOptions::default()
    };
    assert!(ewald_energy_forces(&structure, &charges, &bad).is_err());
    let charges = qeq_charges(
        &water(),
        &ChargeOptions {
            total_charge: f64::NAN,
            ..ChargeOptions::default()
        },
    );
    assert!(charges.is_err());
}
