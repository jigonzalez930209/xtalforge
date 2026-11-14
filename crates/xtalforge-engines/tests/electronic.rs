//! Integration tests for electronic-structure analysis (Task 18.3).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_engines::electronic::{
    BandEdge, DensityOfStates, DosOptions, DosProjectionKind, ELECTRON_MASS_FACTOR,
    ElectronicStructure, KPointLabel, Orbital, OrbitalWeight, Projection, Spin, band_edge_masses,
    band_gap, band_structure, d_band_centre, density_of_states, effective_mass, fit_parabolic_band,
    path_distances, smear_dos, smear_levels, spin_band_edge_masses, spin_band_structures,
    spin_gaps,
};

/// Uniform grid from `min` to `max` with `points` nodes.
fn grid(min: f64, max: f64, points: usize) -> Vec<f64> {
    (0..points)
        .map(|index| min + (max - min) * index as f64 / (points - 1) as f64)
        .collect()
}

/// Trapezoid integral of `values` on `energies`.
fn trapezoid(energies: &[f64], values: &[f64]) -> f64 {
    energies
        .windows(2)
        .zip(values.windows(2))
        .map(|(energy, value)| 0.5 * (value[0] + value[1]) * (energy[1] - energy[0]))
        .sum()
}

/// Builds two bands with explicit occupations and a Fermi energy.
fn two_bands(
    valence: Vec<f64>,
    conduction: Vec<f64>,
    with_occupations: bool,
    with_fermi: bool,
) -> ElectronicStructure {
    let kpoints = (0..valence.len())
        .map(|index| [index as f64 * 0.1, 0.0, 0.0])
        .collect();
    let occupations = vec![vec![2.0; valence.len()], vec![0.0; conduction.len()]];
    let mut input = ElectronicStructure::new(kpoints, vec![valence, conduction]);
    if with_occupations {
        input = input.with_occupations(occupations);
    }
    if with_fermi {
        input = input.with_fermi(0.0);
    }
    input
}

#[test]
fn direct_gap_of_two_bands() {
    let input = two_bands(
        vec![-0.5, -1.0, -1.5, -2.0, -2.5],
        vec![0.5, 1.0, 1.5, 2.0, 2.5],
        true,
        true,
    );
    let structure = band_structure(&input).unwrap();
    assert_eq!(structure.spin, Spin::Up);
    assert_eq!(structure.band_count(), 2);
    assert_eq!(structure.kpoint_count(), 5);
    let gap = structure.gap().unwrap();
    assert!((gap.energy - 1.0).abs() < 1e-12);
    assert!((gap.direct_energy - 1.0).abs() < 1e-12);
    assert!(gap.direct);
    assert!(!gap.metallic);
    assert_eq!(gap.vbm_kpoint, 0);
    assert_eq!(gap.cbm_kpoint, 0);
    assert!(gap.vbm_label.is_none());
    assert_eq!(band_gap(&input).unwrap(), gap);
}

#[test]
fn indirect_gap_of_two_bands() {
    let input = two_bands(
        vec![-0.5, -0.7, -0.9, -1.1, -1.3],
        vec![1.5, 1.0, 0.5, 1.0, 1.5],
        true,
        false,
    );
    let gap = band_gap(&input).unwrap();
    assert!((gap.energy - 1.0).abs() < 1e-12);
    assert!((gap.direct_energy - 1.4).abs() < 1e-12);
    assert!(!gap.direct);
    assert!(!gap.metallic);
    assert!((gap.vbm + 0.5).abs() < 1e-12);
    assert!((gap.cbm - 0.5).abs() < 1e-12);
    assert_eq!(gap.vbm_kpoint, 0);
    assert_eq!(gap.cbm_kpoint, 2);
}

#[test]
fn degenerate_direct_and_metallic_gaps() {
    let direct = two_bands(vec![-1.0, -0.5, -0.5], vec![0.5, 1.0, 0.5], true, false);
    let gap = band_gap(&direct).unwrap();
    assert!(!gap.metallic);
    assert!((gap.energy - 1.0).abs() < 1e-12);
    assert!((gap.direct_energy - 1.0).abs() < 1e-12);
    assert!(gap.direct);

    let metal = two_bands(vec![0.0, 0.0, 0.0], vec![-0.5, -0.5, -0.5], true, false);
    let gap = band_gap(&metal).unwrap();
    assert!(gap.metallic);
    assert_eq!(gap.energy, 0.0);
    assert!(!gap.direct);
}

#[test]
fn spin_resolved_gaps() {
    let input = two_bands(vec![-1.0, -1.0, -1.0], vec![1.0, 1.0, 1.0], true, false)
        .with_spin_down(vec![vec![-1.0, -1.0, -1.0], vec![0.5, 0.5, 0.5]])
        .with_spin_down_occupations(vec![vec![2.0, 2.0, 2.0], vec![0.0, 0.0, 0.0]]);
    assert_eq!(input.spin_count(), 2);
    let structures = spin_band_structures(&input).unwrap();
    assert_eq!(structures.len(), 2);
    assert_eq!(structures[0].spin, Spin::Up);
    assert_eq!(structures[1].spin, Spin::Down);
    assert_eq!(Spin::Up.as_str(), "up");
    assert_eq!(Spin::Down.as_str(), "down");
    let gaps = spin_gaps(&input).unwrap();
    assert_eq!(gaps.len(), 2);
    assert_eq!(gaps[0].spin, Spin::Up);
    assert!((gaps[0].gap.energy - 2.0).abs() < 1e-12);
    assert_eq!(gaps[1].spin, Spin::Down);
    assert!((gaps[1].gap.energy - 1.5).abs() < 1e-12);
    assert_eq!(input.energies_for(Spin::Down).unwrap().len(), 2);
    assert_eq!(input.occupations_for(Spin::Down).unwrap().len(), 2);
}

#[test]
fn effective_mass_of_parabolic_band() {
    let expected = 0.5;
    let coefficient = ELECTRON_MASS_FACTOR / (2.0 * expected);
    let k: Vec<f64> = (0..=12).map(|index| index as f64 * 0.05).collect();
    let conduction: Vec<f64> = k.iter().map(|value| coefficient * value * value).collect();
    let mass = effective_mass(&k, &conduction, 6, 3).unwrap();
    assert!((mass - expected).abs() < 1e-10, "mass = {mass}");

    let fit = fit_parabolic_band(&k, &conduction, 6, 3, BandEdge::Conduction).unwrap();
    assert_eq!(fit.points, 7);
    assert_eq!(fit.edge, BandEdge::Conduction);
    assert!((fit.k - 0.3).abs() < 1e-12);
    assert!((fit.curvature - ELECTRON_MASS_FACTOR / expected).abs() < 1e-9);
    assert!((fit.magnitude() - expected).abs() < 1e-10);

    let valence: Vec<f64> = k.iter().map(|value| -coefficient * value * value).collect();
    let hole = fit_parabolic_band(&k, &valence, 6, 3, BandEdge::Valence).unwrap();
    assert!((hole.mass + expected).abs() < 1e-10, "hole = {}", hole.mass);
    assert_eq!(BandEdge::Valence.as_str(), "valence");
}

#[test]
fn effective_masses_at_band_edges() {
    let coefficient = ELECTRON_MASS_FACTOR;
    let kpoints: Vec<[f64; 3]> = (0..=10)
        .map(|index| [index as f64 * 0.1, 0.0, 0.0])
        .collect();
    let offset = |index: usize| index as f64 * 0.1 - 0.5;
    let conduction: Vec<f64> = (0..=10)
        .map(|index| 0.5 + coefficient * offset(index).powi(2))
        .collect();
    let valence: Vec<f64> = (0..=10)
        .map(|index| -0.5 - coefficient * offset(index).powi(2))
        .collect();
    let occupations = vec![vec![2.0; 11], vec![0.0; 11]];
    let input = ElectronicStructure::new(kpoints, vec![valence.clone(), conduction.clone()])
        .with_occupations(occupations.clone())
        .with_spin_down(vec![valence, conduction])
        .with_spin_down_occupations(occupations);
    let structure = band_structure(&input).unwrap();
    let masses = band_edge_masses(&structure, 3).unwrap();
    assert_eq!(masses.len(), 2);
    assert!((masses[0].mass - 0.5).abs() < 1e-9, "{:?}", masses[0]);
    assert!((masses[1].mass + 0.5).abs() < 1e-9, "{:?}", masses[1]);
    assert_eq!(masses[0].edge, BandEdge::Conduction);
    assert_eq!(masses[1].edge, BandEdge::Valence);

    let spin_masses = spin_band_edge_masses(&input, 3).unwrap();
    assert_eq!(spin_masses.len(), 2);
    assert_eq!(spin_masses[0].spin, Spin::Up);
    assert_eq!(spin_masses[1].spin, Spin::Down);
    assert!((spin_masses[1].electron.mass - 0.5).abs() < 1e-9);
    assert!((spin_masses[1].hole.mass + 0.5).abs() < 1e-9);
}

#[test]
fn gaussian_smearing_normalization() {
    let energies = grid(-10.0, 10.0, 2001);
    let levels = smear_levels(&[0.0], &energies, 0.5).unwrap();
    let area = trapezoid(&energies, &levels);
    assert!((area - 1.0).abs() < 1e-9, "level area = {area}");

    let step = energies[1] - energies[0];
    let center = energies
        .iter()
        .position(|value| value.abs() < 1e-12)
        .unwrap();
    let mut sampled = vec![0.0; energies.len()];
    sampled[center] = 1.0 / step;
    let smoothed = smear_dos(&energies, &sampled, 0.5).unwrap();
    let smoothed_area = trapezoid(&energies, &smoothed);
    assert!(
        (smoothed_area - 1.0).abs() < 1e-7,
        "smoothed area = {smoothed_area}"
    );
    let empty = smear_levels(&[], &energies, 0.5).unwrap();
    assert!(empty.iter().all(|value| *value == 0.0));
}

#[test]
fn d_band_centre_of_known_dos() {
    let energies = vec![-3.0, -2.0, -1.0, 0.0];
    let values = vec![0.0, 1.0, 3.0, 0.0];
    let centre = d_band_centre(&energies, &values).unwrap();
    assert!((centre + 1.25).abs() < 1e-12, "centre = {centre}");

    let fine = grid(-6.0, 6.0, 1201);
    let gaussian = smear_levels(&[-1.5], &fine, 0.5).unwrap();
    let centre = d_band_centre(&fine, &gaussian).unwrap();
    assert!((centre + 1.5).abs() < 1e-6, "gaussian centre = {centre}");
}

#[test]
fn projection_grouping() {
    let input = ElectronicStructure::new(vec![[0.0, 0.0, 0.0]], vec![vec![-1.0]])
        .with_occupations(vec![vec![2.0]])
        .with_projections(vec![Projection {
            band: 0,
            kpoint: 0,
            weights: vec![
                OrbitalWeight {
                    site: 0,
                    element: "Fe".to_owned(),
                    orbital: Orbital::D,
                    weight: 0.6,
                },
                OrbitalWeight {
                    site: 0,
                    element: "Fe".to_owned(),
                    orbital: Orbital::S,
                    weight: 0.2,
                },
                OrbitalWeight {
                    site: 1,
                    element: "O".to_owned(),
                    orbital: Orbital::P,
                    weight: 0.2,
                },
            ],
        }]);
    let options = DosOptions {
        sigma: 0.2,
        points: 2001,
        energy_range: Some((-5.0, 5.0)),
    };
    let dos = density_of_states(&input, &options).unwrap();
    assert!((dos.total_states().unwrap() - 1.0).abs() < 1e-9);
    assert!((trapezoid(&dos.energies, &dos.total) - 1.0).abs() < 1e-9);

    let iron = dos.element("Fe").unwrap();
    assert_eq!(iron.kind, DosProjectionKind::Element);
    assert_eq!(iron.label, "Fe");
    let oxygen = dos.element("O").unwrap();
    assert_eq!(oxygen.kind, DosProjectionKind::Element);
    assert!((trapezoid(&dos.energies, &iron.values) - 0.8).abs() < 1e-9);
    assert!((trapezoid(&dos.energies, &oxygen.values) - 0.2).abs() < 1e-9);

    let d = dos.orbital(Orbital::D).unwrap();
    assert_eq!(d.kind, DosProjectionKind::Orbital);
    assert_eq!(d.orbital.as_ref(), Some(&Orbital::D));
    assert!((trapezoid(&dos.energies, &d.values) - 0.6).abs() < 1e-9);
    assert!(
        (trapezoid(&dos.energies, &dos.orbital(Orbital::S).unwrap().values) - 0.2).abs() < 1e-9
    );
    assert!(
        (trapezoid(&dos.energies, &dos.orbital(Orbital::P).unwrap().values) - 0.2).abs() < 1e-9
    );
    assert!((dos.d_band_centre().unwrap() + 1.0).abs() < 1e-9);
    assert!(dos.projection("Fe").is_some());
    assert!(dos.projection("Cu").is_none());
    assert_eq!(dos.total_buffer().values.len(), dos.energies.len());
    assert_eq!(d.buffer(&dos.energies).energies, dos.energies);
    assert_eq!(Orbital::parse("d"), Orbital::D);
    assert_eq!(Orbital::parse("dxy").as_str(), "dxy");

    let json = serde_json::to_string(&dos).unwrap();
    let decoded: DensityOfStates = serde_json::from_str(&json).unwrap();
    assert_eq!(decoded.energies.len(), dos.energies.len());
    for (decoded_energy, energy) in decoded.energies.iter().zip(&dos.energies) {
        assert!((decoded_energy - energy).abs() < 1e-12);
    }
    assert_eq!(decoded.projections.len(), dos.projections.len());
    for (decoded_projection, original) in decoded.projections.iter().zip(&dos.projections) {
        assert_eq!(decoded_projection.label, original.label);
        assert_eq!(decoded_projection.kind, original.kind);
        for (decoded_value, value) in decoded_projection.values.iter().zip(&original.values) {
            assert!((decoded_value - value).abs() < 1e-12);
        }
    }
}

#[test]
fn band_path_metadata_and_serialization() {
    let kpoints = vec![[0.0, 0.0, 0.0], [0.25, 0.0, 0.0], [0.5, 0.0, 0.0]];
    let input = ElectronicStructure::new(kpoints.clone(), vec![vec![-1.0, -0.8, -0.6]])
        .with_labels(vec![
            KPointLabel {
                index: 0,
                label: "Gamma".to_owned(),
            },
            KPointLabel {
                index: 2,
                label: "X".to_owned(),
            },
        ])
        .with_distances(vec![0.0, 1.0, 2.0])
        .with_fermi(-0.7);
    assert_eq!(input.kpoint_count(), 3);
    assert_eq!(input.band_count(), 1);
    assert_eq!(input.spin_count(), 1);
    let structure = band_structure(&input).unwrap();
    structure.validate().unwrap();
    assert_eq!(structure.label_at(0), Some("Gamma"));
    assert_eq!(structure.label_at(1), None);
    assert_eq!(structure.energy(0, 1), Some(-0.8));
    assert_eq!(structure.energy(1, 0), None);
    assert_eq!(structure.distances, vec![0.0, 1.0, 2.0]);
    assert_eq!(structure.fermi, Some(-0.7));

    let json = serde_json::to_string(&structure).unwrap();
    let decoded: xtalforge_engines::electronic::BandStructure =
        serde_json::from_str(&json).unwrap();
    assert_eq!(decoded.spin, structure.spin);
    assert_eq!(decoded.labels, structure.labels);
    assert_eq!(decoded.kpoints, structure.kpoints);
    assert_eq!(decoded.distances, structure.distances);
    assert_eq!(decoded.fermi, structure.fermi);
    for (decoded_band, band) in decoded.bands.iter().zip(&structure.bands) {
        for (decoded_energy, energy) in decoded_band.energies.iter().zip(&band.energies) {
            assert!((decoded_energy - energy).abs() < 1e-12);
        }
    }

    let automatic = path_distances(&[[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [0.5, 0.5, 0.0]]);
    assert!((automatic[1] - 0.5).abs() < 1e-12);
    assert!((automatic[2] - 1.0).abs() < 1e-12);
}

#[test]
fn invalid_inputs_are_rejected() {
    let mismatched =
        ElectronicStructure::new(vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], vec![vec![0.0]]);
    assert!(mismatched.validate().is_err());
    assert!(band_structure(&mismatched).is_err());

    let nan = ElectronicStructure::new(vec![[0.0, 0.0, 0.0]], vec![vec![f64::NAN]]);
    assert!(nan.validate().is_err());

    let bad_occupations = ElectronicStructure::new(vec![[0.0, 0.0, 0.0]], vec![vec![0.0]])
        .with_occupations(vec![vec![2.5]]);
    assert!(bad_occupations.validate().is_err());

    let bad_labels =
        ElectronicStructure::new(vec![[0.0, 0.0, 0.0]], vec![vec![0.0]]).with_labels(vec![
            KPointLabel {
                index: 4,
                label: "X".to_owned(),
            },
        ]);
    assert!(bad_labels.validate().is_err());

    let bad_distances =
        ElectronicStructure::new(vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], vec![vec![0.0, 1.0]])
            .with_distances(vec![1.0, 0.0]);
    assert!(bad_distances.validate().is_err());

    let bad_projection = ElectronicStructure::new(vec![[0.0, 0.0, 0.0]], vec![vec![0.0]])
        .with_projections(vec![Projection {
            band: 3,
            kpoint: 0,
            weights: Vec::new(),
        }]);
    assert!(bad_projection.validate().is_err());

    let negative_weight = ElectronicStructure::new(vec![[0.0, 0.0, 0.0]], vec![vec![0.0]])
        .with_projections(vec![Projection {
            band: 0,
            kpoint: 0,
            weights: vec![OrbitalWeight {
                site: 0,
                element: "Fe".to_owned(),
                orbital: Orbital::D,
                weight: -1.0,
            }],
        }]);
    assert!(negative_weight.validate().is_err());

    let empty_orbital = ElectronicStructure::new(vec![[0.0, 0.0, 0.0]], vec![vec![0.0]])
        .with_projections(vec![Projection {
            band: 0,
            kpoint: 0,
            weights: vec![OrbitalWeight {
                site: 0,
                element: "Fe".to_owned(),
                orbital: Orbital::Other(String::new()),
                weight: 1.0,
            }],
        }]);
    assert!(empty_orbital.validate().is_err());

    let spin_without_channel = ElectronicStructure::new(vec![[0.0, 0.0, 0.0]], vec![vec![0.0]])
        .with_spin_down_occupations(vec![vec![0.0]]);
    assert!(spin_without_channel.validate().is_err());

    let unclassified = ElectronicStructure::new(
        vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        vec![vec![-1.0, 1.0]],
    );
    assert!(unclassified.validate().is_ok());
    assert!(band_gap(&unclassified).is_err());

    let all_occupied = ElectronicStructure::new(
        vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        vec![vec![-1.0, -0.5]],
    )
    .with_occupations(vec![vec![2.0, 2.0]]);
    assert!(band_gap(&all_occupied).is_err());

    let k = [0.0, 0.1, 0.2];
    assert!(effective_mass(&k, &[0.0, 0.1], 0, 1).is_err());
    assert!(effective_mass(&k[..2], &[0.0, 0.1], 0, 1).is_err());
    assert!(effective_mass(&k, &[0.0, 0.1, 0.2], 3, 1).is_err());
    assert!(effective_mass(&k, &[0.0, 0.1, 0.2], 1, 0).is_err());
    assert!(effective_mass(&[0.0, 0.2, 0.1], &[0.0, 0.1, 0.2], 1, 1).is_err());
    assert!(effective_mass(&k, &[1.0, 1.0, 1.0], 1, 1).is_err());

    let fine = grid(-2.0, 2.0, 41);
    assert!(smear_levels(&[0.0], &fine, 0.0).is_err());
    assert!(smear_levels(&[0.0], &[0.0, 0.0], 0.1).is_err());
    assert!(smear_dos(&fine, &vec![0.0; 40], 0.1).is_err());
    assert!(smear_dos(&fine, &vec![0.0; 41], f64::NAN).is_err());
    assert!(d_band_centre(&fine, &vec![0.0; 41]).is_err());
    assert!(d_band_centre(&[0.0, 1.0], &[-1.0, 1.0]).is_err());

    assert!(DosOptions::new(0.0).validate().is_err());
    assert!(
        DosOptions {
            sigma: 0.1,
            points: 1,
            energy_range: None,
        }
        .validate()
        .is_err()
    );
    assert!(
        DosOptions {
            sigma: 0.1,
            points: 10,
            energy_range: Some((1.0, 0.0)),
        }
        .validate()
        .is_err()
    );

    let empty = DensityOfStates {
        energies: fine.clone(),
        total: vec![0.0; fine.len()],
        projections: Vec::new(),
        fermi: None,
    };
    assert!(empty.d_band_centre().is_err());
    let broken = DensityOfStates {
        energies: fine,
        total: vec![0.0; 3],
        projections: Vec::new(),
        fermi: None,
    };
    assert!(broken.total_states().is_err());
}
