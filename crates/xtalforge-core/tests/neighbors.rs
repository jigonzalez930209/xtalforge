//! Integration tests for the periodic neighbor search (Task 1.10).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::{Lattice, Site, StructureBuilder, neighbors, neighbors_of};

fn structure(sites: &[(&str, [f64; 3])]) -> xtalforge_core::Structure {
    let mut builder = StructureBuilder::new("test", Lattice::cubic(4.0).unwrap());
    for (index, (element, frac)) in sites.iter().enumerate() {
        builder = builder.site(Site::new(format!("{element}{index}"), element, *frac).unwrap());
    }
    builder.build().unwrap()
}

#[test]
fn nearest_pair_and_image() {
    let structure = structure(&[("C", [0.0, 0.0, 0.0]), ("C", [0.5, 0.0, 0.0])]);
    let found = neighbors(&structure, 2.1).unwrap();
    assert_eq!(found.len(), 2);
    assert!(
        found
            .iter()
            .all(|neighbor| neighbor.i == 0 && neighbor.j == 1)
    );
    assert!(found.iter().any(|neighbor| neighbor.image == [0, 0, 0]));
    assert!(found.iter().any(|neighbor| neighbor.image == [-1, 0, 0]));
    assert!(
        found
            .iter()
            .all(|neighbor| (neighbor.distance - 2.0).abs() < 1e-12)
    );
}

#[test]
fn self_images_are_reported_once() {
    let structure = structure(&[("C", [0.0, 0.0, 0.0])]);
    let found = neighbors(&structure, 4.1).unwrap();
    assert_eq!(found.len(), 6);
    for neighbor in &found {
        assert_eq!(neighbor.i, neighbor.j);
        assert!((neighbor.distance - 4.0).abs() < 1e-12);
    }
}

#[test]
fn periodic_image_between_boundary_atoms() {
    let structure = structure(&[("C", [0.9, 0.0, 0.0]), ("C", [0.1, 0.0, 0.0])]);
    let found = neighbors(&structure, 1.0).unwrap();
    assert_eq!(found.len(), 1);
    assert_eq!(found[0].image, [1, 0, 0]);
    assert!((found[0].distance - 0.8).abs() < 1e-12);
}

#[test]
fn results_are_sorted_and_filterable() {
    let structure = structure(&[
        ("C", [0.0, 0.0, 0.0]),
        ("C", [0.25, 0.0, 0.0]),
        ("C", [0.5, 0.0, 0.0]),
    ]);
    let found = neighbors(&structure, 1.1).unwrap();
    assert!(found.windows(2).all(|pair| {
        (pair[0].i, pair[0].j, pair[0].image) <= (pair[1].i, pair[1].j, pair[1].image)
    }));
    let first = neighbors_of(&structure, 0, 1.1).unwrap();
    assert!(
        first
            .iter()
            .all(|neighbor| neighbor.i == 0 || neighbor.j == 0)
    );
}

#[test]
fn cutoff_is_validated() {
    let structure = structure(&[("C", [0.0, 0.0, 0.0])]);
    assert!(neighbors(&structure, 0.0).is_err());
    assert!(neighbors(&structure, f64::NAN).is_err());
}

#[test]
fn linked_cells_match_the_direct_search_on_a_triclinic_cell() {
    let lattice = Lattice::new([[6.1, 0.0, 0.0], [1.3, 5.2, 0.0], [0.8, 1.1, 4.9]]).unwrap();
    let mut seed = 0x5eed_u64;
    let mut builder = StructureBuilder::new("random", lattice);
    for index in 0..32 {
        let frac = [
            pseudo_random(&mut seed) * 1.8 - 0.4,
            pseudo_random(&mut seed) * 1.8 - 0.4,
            pseudo_random(&mut seed) * 1.8 - 0.4,
        ];
        builder = builder.site(Site::new(format!("C{index}"), "C", frac).unwrap());
    }
    let structure = builder.build().unwrap();

    for cutoff in [1.7, 2.4, 3.6] {
        let linked = neighbors(&structure, cutoff).unwrap();
        let direct = direct_neighbors(&structure, cutoff);
        assert_eq!(linked.len(), direct.len(), "cutoff {cutoff}");
        for (found, expected) in linked.iter().zip(direct.iter()) {
            assert_eq!(
                (found.i, found.j, found.image),
                (expected.i, expected.j, expected.image),
                "cutoff {cutoff}"
            );
            assert!(
                (found.distance - expected.distance).abs() < 1e-12,
                "cutoff {cutoff}"
            );
        }
    }
}

fn pseudo_random(seed: &mut u64) -> f64 {
    *seed = seed
        .wrapping_mul(6_364_136_223_846_793_005)
        .wrapping_add(1_442_695_040_888_963_407);
    ((*seed >> 11) as f64) / ((1_u64 << 53) as f64)
}

fn direct_neighbors(
    structure: &xtalforge_core::Structure,
    cutoff: f64,
) -> Vec<xtalforge_core::Neighbor> {
    let p1 = structure.p1().unwrap();
    let reciprocal = p1.lattice.reciprocal(false).unwrap();
    let minimum_spacing = reciprocal
        .matrix()
        .iter()
        .map(|row| 1.0 / row.iter().map(|value| value * value).sum::<f64>().sqrt())
        .fold(f64::INFINITY, f64::min);
    let range = (cutoff / minimum_spacing).ceil() as i32 + 1;
    let matrix = p1.lattice.matrix();
    let mut found = Vec::new();
    for (i, atom_i) in p1.atoms.iter().enumerate() {
        for (j, atom_j) in p1.atoms.iter().enumerate().skip(i) {
            for x in -range..=range {
                for y in -range..=range {
                    for z in -range..=range {
                        if i == j && x == 0 && y == 0 && z == 0 {
                            continue;
                        }
                        let delta = [
                            atom_j.frac[0] + f64::from(x) - atom_i.frac[0],
                            atom_j.frac[1] + f64::from(y) - atom_i.frac[1],
                            atom_j.frac[2] + f64::from(z) - atom_i.frac[2],
                        ];
                        let distance = norm(multiply_row(delta, matrix));
                        if distance <= cutoff + 1e-12 {
                            found.push(xtalforge_core::Neighbor {
                                i,
                                j,
                                image: [x, y, z],
                                distance,
                            });
                        }
                    }
                }
            }
        }
    }
    found.sort_by(|first, second| {
        (first.i, first.j, first.image, first.distance.to_bits()).cmp(&(
            second.i,
            second.j,
            second.image,
            second.distance.to_bits(),
        ))
    });
    found
}

fn multiply_row(vector: [f64; 3], matrix: [[f64; 3]; 3]) -> [f64; 3] {
    let [x, y, z] = vector;
    let [[a, b, c], [d, e, f], [g, h, i]] = matrix;
    [
        x * a + y * d + z * g,
        x * b + y * e + z * h,
        x * c + y * f + z * i,
    ]
}

fn norm(vector: [f64; 3]) -> f64 {
    vector.iter().map(|value| value * value).sum::<f64>().sqrt()
}
