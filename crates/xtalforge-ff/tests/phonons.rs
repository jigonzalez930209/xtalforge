//! Phonon q-mesh and DOS against an analytic nearest-neighbor spring model.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::{Lattice, Site, Structure, StructureBuilder};
use xtalforge_ff::vibrations::{DosOptions, PhononOptions, phonons};
use xtalforge_ff::{Calculator, Error, UffCalculator};

/// A simple cubic lattice joined by nearest-neighbour springs.
struct SpringCalculator {
    spring: f64,
    equilibrium: f64,
}

impl Calculator for SpringCalculator {
    fn name(&self) -> &'static str {
        "spring"
    }

    fn energy(&self, structure: &Structure) -> Result<f64, Error> {
        Ok(self.evaluate(structure)?.0)
    }

    fn energy_forces(&self, structure: &Structure) -> Result<(f64, Vec<[f64; 3]>), Error> {
        self.evaluate(structure)
    }
}

impl SpringCalculator {
    fn evaluate(&self, structure: &Structure) -> Result<(f64, Vec<[f64; 3]>), Error> {
        let p1 = structure.p1()?;
        let matrix = p1.lattice.matrix();
        let mut forces = vec![[0.0; 3]; p1.atoms.len()];
        let mut energy = 0.0;
        for (i, atom) in p1.atoms.iter().enumerate() {
            for (j, other) in p1.atoms.iter().enumerate().skip(i + 1) {
                let delta = [
                    other.frac[0] - atom.frac[0],
                    other.frac[1] - atom.frac[1],
                    other.frac[2] - atom.frac[2],
                ];
                let delta = [
                    delta[0] - delta[0].round(),
                    delta[1] - delta[1].round(),
                    delta[2] - delta[2].round(),
                ];
                let cart = [
                    delta[0] * matrix[0][0] + delta[1] * matrix[1][0] + delta[2] * matrix[2][0],
                    delta[0] * matrix[0][1] + delta[1] * matrix[1][1] + delta[2] * matrix[2][1],
                    delta[0] * matrix[0][2] + delta[1] * matrix[1][2] + delta[2] * matrix[2][2],
                ];
                let length = (cart[0] * cart[0] + cart[1] * cart[1] + cart[2] * cart[2]).sqrt();
                if (length - self.equilibrium).abs() > 0.1 {
                    continue;
                }
                let extension = length - self.equilibrium;
                energy += 0.5 * self.spring * extension * extension;
                if length > 0.0 {
                    let unit = [cart[0] / length, cart[1] / length, cart[2] / length];
                    for axis in 0..3 {
                        forces[i][axis] += self.spring * extension * unit[axis];
                        forces[j][axis] -= self.spring * extension * unit[axis];
                    }
                }
            }
        }
        Ok((energy, forces))
    }
}

fn spring_structure() -> Structure {
    StructureBuilder::new("spring", Lattice::cubic(4.0).unwrap())
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .build()
        .unwrap()
}

fn expected(k: f64, mass: f64, q: [f64; 3]) -> Vec<f64> {
    let factor = 521.470_898_372_506_6;
    let mut values = [q[0], q[1], q[2]]
        .iter()
        .map(|value| {
            let omega = (4.0 * k / mass).sqrt() * (std::f64::consts::PI * value).sin().abs();
            omega * factor
        })
        .collect::<Vec<_>>();
    values.sort_by(f64::total_cmp);
    values
}

#[test]
fn q_mesh_matches_the_analytic_spring_dispersion() {
    let calculator = SpringCalculator {
        spring: 10.0,
        equilibrium: 4.0,
    };
    let structure = spring_structure();
    let options = PhononOptions {
        q_grid: [2, 2, 2],
        supercell: [3, 3, 3],
        dos: Some(DosOptions::default()),
        ..PhononOptions::default()
    };
    let result = phonons(&calculator, &structure, &options).unwrap();
    assert_eq!(result.q_points.len(), 8);
    let mass = 12.011;
    for (q, found) in result.q_points.iter().zip(result.frequencies.iter()) {
        let mut found = found.clone();
        found.sort_by(f64::total_cmp);
        let expected = expected(10.0, mass, *q);
        for (got, want) in found.iter().zip(expected.iter()) {
            assert!((got - want).abs() < 0.1, "q={q:?}: {got} != {want}");
        }
    }
    let dos = result.dos.expect("DOS");
    assert_eq!(dos.frequencies.len(), 200);
    assert!(dos.density.iter().all(|value| value.is_finite()));
}

#[test]
fn a_q_mesh_needs_a_supercell() {
    let options = PhononOptions {
        q_grid: [2, 2, 2],
        ..PhononOptions::default()
    };
    assert!(options.validate().is_err());
}

#[test]
fn gamma_remains_the_primitive_hessian() {
    let calculator = UffCalculator::with_defaults().unwrap();
    let structure = StructureBuilder::new("methane", Lattice::cubic(10.0).unwrap())
        .site(Site::new("C1", "C", [0.5, 0.5, 0.5]).unwrap())
        .build()
        .unwrap();
    let result = phonons(&calculator, &structure, &PhononOptions::default()).unwrap();
    assert_eq!(result.q_points, vec![[0.0, 0.0, 0.0]]);
    assert_eq!(result.frequencies.len(), 1);
}

#[test]
fn q_mesh_matches_the_phonopy_spring_golden() {
    let golden: serde_json::Value = serde_json::from_str(include_str!(
        "../../../tests/golden/phase-4/phonons/spring-sc.json"
    ))
    .unwrap();
    let spring = golden["spring_ev_per_a2"].as_f64().unwrap();
    let equilibrium = golden["equilibrium_a"].as_f64().unwrap();
    let calculator = SpringCalculator {
        spring,
        equilibrium,
    };
    let structure = spring_structure();
    let options = PhononOptions {
        q_grid: [2, 2, 2],
        supercell: [3, 3, 3],
        dos: None,
        ..PhononOptions::default()
    };
    let result = phonons(&calculator, &structure, &options).unwrap();
    let points = golden["q_points"].as_array().unwrap();
    let expected = golden["frequencies_cm1"].as_array().unwrap();
    for (q, found) in result.q_points.iter().zip(result.frequencies.iter()) {
        let index = points
            .iter()
            .position(|point| {
                let values = point.as_array().unwrap();
                (0..3).all(|axis| (values[axis].as_f64().unwrap() - q[axis]).abs() < 1e-12)
            })
            .expect("golden q point");
        let mut found = found.clone();
        found.sort_by(f64::total_cmp);
        let mut want = expected[index]
            .as_array()
            .unwrap()
            .iter()
            .map(|value| value.as_f64().unwrap())
            .collect::<Vec<_>>();
        want.sort_by(f64::total_cmp);
        for (got, reference) in found.iter().zip(want.iter()) {
            assert!(
                (got - reference).abs() < 0.5,
                "q={q:?}: {got} != {reference}"
            );
        }
    }
}

#[test]
fn q_path_replaces_the_monkhorst_pack_grid() {
    let calculator = SpringCalculator {
        spring: 10.0,
        equilibrium: 4.0,
    };
    let structure = spring_structure();
    let options = PhononOptions {
        supercell: [3, 3, 3],
        q_path: vec![[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]],
        ..PhononOptions::default()
    };
    let result = phonons(&calculator, &structure, &options).unwrap();
    assert_eq!(result.q_points, options.q_path);
    assert_eq!(result.frequencies.len(), 2);
}
