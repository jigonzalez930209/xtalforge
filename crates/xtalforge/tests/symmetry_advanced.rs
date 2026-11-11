//! Integration tests for the symmetry workflow and analysis (Task 2.3-2.10).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_io::read_cif;
use xtalforge_symmetry::{SpaceGroup, SymmetryOptions, detect_symmetry};

const SI: &str = include_str!("../../../tests/corpus/simple/Si.cif");

#[test]
fn detect_asymmetrize_and_p1_round_trip_silicon() {
    let structure = read_cif(SI).unwrap();
    let result =
        xtalforge_symmetry::detect_and_symmetrize(&structure, SymmetryOptions::default()).unwrap();
    assert_eq!(result.structure.sites().len(), 1);
    assert_eq!(result.structure.hall(), "F 4d 2 3 -1d");
    assert_eq!(result.structure.operations().len(), 192);
    let expanded = result.structure.p1().unwrap();
    assert_eq!(expanded.n_atoms(), 8);
    let p1 = xtalforge_symmetry::to_p1(&structure).unwrap();
    assert!(p1.is_p1());
    assert_eq!(p1.sites().len(), 8);
}

#[test]
fn impose_averages_a_distorted_nacl_cell() {
    let lattice = xtalforge_core::Lattice::cubic(5.64).unwrap();
    let distorted = xtalforge_core::StructureBuilder::new("NaCl", lattice)
        .site(xtalforge_core::Site::new("Na1", "Na", [0.0, 0.0, 0.0]).unwrap())
        .site(xtalforge_core::Site::new("Cl1", "Cl", [0.51, 0.5, 0.5]).unwrap())
        .build()
        .unwrap();
    let group = SpaceGroup::from_number(225).unwrap();
    let report = xtalforge_symmetry::set_space_group(
        &distorted,
        &group,
        xtalforge_symmetry::SetMode::Impose,
    )
    .unwrap();
    assert_eq!(report.n_input, 2);
    assert_eq!(report.n_output, 2);
    assert!(report.max_displacement > 0.0);
    let dataset = detect_symmetry(&report.structure, SymmetryOptions::default()).unwrap();
    assert_eq!(dataset.number, 225);
    let summary = report.summary(&group);
    assert_eq!(summary.number, 225);
    assert_eq!(summary.hall, "-F 4 2 3");
}

#[test]
fn reinterpret_keeps_the_sites() {
    let structure = read_cif(SI).unwrap();
    let asymmetric =
        xtalforge_symmetry::asymmetrize(&structure, SymmetryOptions::default()).unwrap();
    let group = SpaceGroup::from_number(227).unwrap();
    let report = xtalforge_symmetry::set_space_group(
        &asymmetric,
        &group,
        xtalforge_symmetry::SetMode::Reinterpret,
    )
    .unwrap();
    assert_eq!(report.n_output, 1);
    assert_eq!(report.max_displacement, 0.0);
    assert_eq!(report.structure.sites().len(), 1);
}

#[test]
fn wyckoff_assignments_carry_multiplicity_and_parameters() {
    let structure = read_cif(SI).unwrap();
    let assignments = xtalforge_symmetry::wyckoff_assignments(
        &structure,
        SymmetryOptions {
            symprec: 0.001,
            ..SymmetryOptions::default()
        },
    )
    .unwrap();
    assert_eq!(assignments.len(), 8);
    assert!(
        assignments
            .iter()
            .all(|assignment| assignment.letter == 'a' && assignment.coordinates == "0,0,0")
    );
    assert!(
        assignments
            .iter()
            .all(|assignment| assignment.multiplicity == 8)
    );
    assert!(
        assignments
            .iter()
            .all(|assignment| !assignment.coordinates.is_empty())
    );
    let parameters = xtalforge_symmetry::free_parameters("x,x,z", [0.1, 0.2, 0.3]);
    assert_eq!(parameters, vec![('x', 0.1), ('z', 0.3)]);
    let parameters = xtalforge_symmetry::free_parameters("1/4,y,1/4", [0.0, 0.2, 0.0]);
    assert_eq!(parameters, vec![('y', 0.2)]);
}

#[test]
fn standardize_and_primitive_silicon() {
    let structure = read_cif(SI).unwrap();
    let standard = xtalforge_symmetry::standardize(&structure, SymmetryOptions::default()).unwrap();
    assert_eq!(standard.structure.sites().len(), 1);
    assert_eq!(standard.structure.operations().len(), 192);
    assert!((standard.structure.lattice().volume() - 160.186_257_298_018_75).abs() < 1e-6);
    assert_eq!(standard.structure.p1().unwrap().n_atoms(), 8);
    assert_eq!(standard.dataset.number, 227);

    let primitive =
        xtalforge_symmetry::primitive_cell(&structure, SymmetryOptions::default()).unwrap();
    assert!((primitive.structure.lattice().volume() - 40.046_564_324_504_69).abs() < 1e-6);
    assert_eq!(primitive.structure.p1().unwrap().n_atoms(), 2);
    assert!(primitive.structure.sites().len() <= 2);
    let ideal = xtalforge_symmetry::idealize(&structure, SymmetryOptions::default()).unwrap();
    assert_eq!(ideal.structure.sites().len(), 1);
}

#[test]
fn reductions_preserve_volume_and_atoms() {
    let lattice =
        xtalforge_core::Lattice::new([[3.0, 0.0, 0.0], [2.8, 1.0, 0.0], [2.8, 0.0, 1.0]]).unwrap();
    let structure = xtalforge_core::StructureBuilder::new("skew", lattice)
        .site(xtalforge_core::Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .build()
        .unwrap();
    for reduction in [
        xtalforge_symmetry::Reduction::Niggli,
        xtalforge_symmetry::Reduction::Delaunay,
        xtalforge_symmetry::Reduction::Minkowski,
    ] {
        let reduced = xtalforge_symmetry::reduce_structure(&structure, reduction).unwrap();
        assert!((reduced.structure.lattice().volume() - structure.lattice().volume()).abs() < 1e-9);
        assert_eq!(reduced.structure.sites().len(), 1);
        assert!((reduced.determinant.abs() - 1.0).abs() < 1e-9);
    }
    let equivalent = xtalforge_symmetry::lattices_equivalent(
        structure.lattice(),
        &structure
            .lattice()
            .transform([[1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
            .unwrap(),
        1e-6,
    )
    .unwrap();
    assert!(equivalent);
}

#[test]
fn subgroup_enumeration_and_identification() {
    let primitive = xtalforge_symmetry::primitive_operations(525).unwrap();
    assert_eq!(primitive.len(), 48);
    let (number, hall) = xtalforge_symmetry::identify_space_group(&primitive).unwrap();
    assert_eq!(number, 227);
    assert!(matches!(hall, 525 | 526));

    let groups = xtalforge_symmetry::translationengleiche_subgroups(525, 1e-5).unwrap();
    assert!(!groups.is_empty());
    assert!(
        groups
            .iter()
            .any(|group| group.index == 1 && group.depth == 0)
    );
    assert!(
        groups
            .iter()
            .all(|group| group.kind == xtalforge_symmetry::SubgroupKind::Translationengleiche)
    );

    let classes =
        xtalforge_symmetry::klassengleiche_subgroups(525, [[1, 0, 0], [0, 1, 0], [0, 0, 1]], 1e-5)
            .unwrap();
    assert!(!classes.is_empty());
    assert!(classes.iter().any(|group| group.index == 1));
}

#[test]
fn tolerance_ladder_and_diagnostics() {
    let structure = read_cif(SI).unwrap();
    let ladder =
        xtalforge_symmetry::tolerance_ladder(&structure, &xtalforge_symmetry::DEFAULT_LADDER)
            .unwrap();
    assert!(!ladder.is_empty());
    assert!(ladder.iter().any(|range| range.number == 227));
    assert!(ladder.iter().all(|range| range.min <= range.max));

    let group = SpaceGroup::from_number(227).unwrap();
    let diagnostics = xtalforge_symmetry::diagnose(&structure, &group).unwrap();
    assert!(diagnostics.max_displacement < 1e-6);
    assert_eq!(diagnostics.atoms.len(), 8);
}

#[test]
fn fcc_brillouin_zone_is_a_truncated_octahedron() {
    let fcc =
        xtalforge_core::Lattice::new([[0.0, 2.0, 2.0], [2.0, 0.0, 2.0], [2.0, 2.0, 0.0]]).unwrap();
    let zone = xtalforge_symmetry::brillouin_zone(&fcc).unwrap();
    assert_eq!(zone.n_faces(), 14);
    assert_eq!(zone.n_vertices(), 24);
    assert_eq!(zone.edges.len(), 36);
    let distances: Vec<f64> = zone
        .vertices
        .iter()
        .map(|vertex| {
            (vertex[0] * vertex[0] + vertex[1] * vertex[1] + vertex[2] * vertex[2]).sqrt()
        })
        .collect();
    assert!(distances.iter().all(|distance| *distance > 0.0));
}

#[test]
fn lowering_to_a_subgroup_reduces_the_asymmetric_unit() {
    let structure = read_cif(SI).unwrap();
    let lowered = xtalforge_symmetry::lower_to_subgroup(&structure, 216, 0.5).unwrap();
    assert_eq!(lowered.operations().len(), 96);
    assert!(!lowered.sites().is_empty());
    assert!(lowered.sites().len() <= 8);
    assert!(lowered.p1().unwrap().n_atoms() >= 8);
    assert!(xtalforge_symmetry::lower_to_subgroup(&structure, 216, 0.0).is_err());
}

fn pseudo_random(seed: &mut u64) -> f64 {
    *seed = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    ((*seed >> 11) as f64) / ((1u64 << 53) as f64)
}

#[test]
fn generated_structures_recover_their_group() {
    let cases = [
        (
            1u32,
            xtalforge_core::Lattice::orthorhombic(5.0, 6.0, 7.0).unwrap(),
        ),
        (
            2,
            xtalforge_core::Lattice::orthorhombic(5.0, 6.0, 7.0).unwrap(),
        ),
        (
            3,
            xtalforge_core::Lattice::from_parameters(5.0, 6.0, 7.0, 90.0, 100.0, 90.0).unwrap(),
        ),
        (
            16,
            xtalforge_core::Lattice::orthorhombic(5.0, 6.0, 7.0).unwrap(),
        ),
        (
            75,
            xtalforge_core::Lattice::from_parameters(4.0, 4.0, 6.0, 90.0, 90.0, 90.0).unwrap(),
        ),
        (
            143,
            xtalforge_core::Lattice::from_parameters(4.0, 4.0, 6.0, 90.0, 90.0, 120.0).unwrap(),
        ),
        (
            168,
            xtalforge_core::Lattice::from_parameters(4.0, 4.0, 6.0, 90.0, 90.0, 120.0).unwrap(),
        ),
        (195, xtalforge_core::Lattice::cubic(6.0).unwrap()),
    ];
    let mut seed = 7u64;
    for (number, lattice) in cases {
        let group = SpaceGroup::from_number(number).unwrap();
        let mut builder = xtalforge_core::StructureBuilder::new("generated", lattice)
            .hall(group.hall_symbol.clone());
        for operation in &group.operations {
            builder = builder.operation(*operation);
        }
        for (index, element) in ["C", "N", "O"].iter().enumerate() {
            let frac = [
                0.13 + 0.07 * f64::from(index as u32) + 0.03 * pseudo_random(&mut seed),
                0.21 + 0.05 * f64::from(index as u32) + 0.02 * pseudo_random(&mut seed),
                0.31 + 0.09 * f64::from(index as u32) + 0.04 * pseudo_random(&mut seed),
            ];
            builder = builder.site(
                xtalforge_core::Site::new(format!("{element}{index}"), element, frac).unwrap(),
            );
        }
        let structure = builder.build().unwrap();
        let p1 = structure.p1().unwrap();
        assert_eq!(p1.n_atoms(), 3 * group.operations.len(), "group {number}");
        let dataset = detect_symmetry(&structure, SymmetryOptions::default()).unwrap();
        assert_eq!(dataset.number, number, "detected group for {number}");
        let reduced =
            xtalforge_symmetry::asymmetrize(&structure, SymmetryOptions::default()).unwrap();
        assert_eq!(
            reduced.p1().unwrap().n_atoms(),
            p1.n_atoms(),
            "round trip {number}"
        );
    }
}

#[test]
fn standardization_is_idempotent() {
    let structure = read_cif(SI).unwrap();
    let first = xtalforge_symmetry::standardize(&structure, SymmetryOptions::default()).unwrap();
    let second =
        xtalforge_symmetry::standardize(&first.structure, SymmetryOptions::default()).unwrap();
    assert!(
        (first.structure.lattice().volume() - second.structure.lattice().volume()).abs() < 1e-9
    );
    assert_eq!(
        first.structure.sites().len(),
        second.structure.sites().len()
    );
}

#[test]
fn layer_group_detection_for_a_square_layer() {
    let lattice =
        xtalforge_core::Lattice::from_parameters(4.0, 4.0, 20.0, 90.0, 90.0, 90.0).unwrap();
    let structure = xtalforge_core::StructureBuilder::new("layer", lattice)
        .site(xtalforge_core::Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .site(xtalforge_core::Site::new("C2", "C", [0.5, 0.5, 0.0]).unwrap())
        .build()
        .unwrap();
    let layer = xtalforge_symmetry::detect_layer_symmetry(&structure, 1e-3).unwrap();
    assert!((1..=80).contains(&layer.number));
    assert!(!layer.operations.is_empty());
    assert_eq!(layer.wyckoffs.len(), 2);
    assert!(xtalforge_symmetry::detect_layer_symmetry(&structure, 0.0).is_err());
}

#[test]
fn collinear_magnetic_search_returns_a_uni_number() {
    let lattice = xtalforge_core::Lattice::cubic(4.0).unwrap();
    let structure = xtalforge_core::StructureBuilder::new("magnetic", lattice)
        .site(
            xtalforge_core::Site::new("Fe1", "Fe", [0.0, 0.0, 0.0])
                .unwrap()
                .with_prop("magmom", serde_json::json!(2.0)),
        )
        .build()
        .unwrap();
    let dataset = xtalforge_symmetry::detect_magnetic_symmetry(&structure, 1e-3).unwrap();
    assert!(dataset.uni_number >= 1);
    assert!(!dataset.operations.is_empty());
    assert!(xtalforge_symmetry::detect_magnetic_symmetry(&structure, 0.0).is_err());
}

#[test]
fn change_of_setting_uses_the_requested_hall_number() {
    let structure = read_cif(SI).unwrap();
    let changed = xtalforge_symmetry::change_setting(&structure, 526).unwrap();
    assert_eq!(changed.structure.hall(), "-F 4vw 2vw 3");
    assert!((changed.structure.lattice().volume() - 160.186_257_298_018_75).abs() < 1e-6);
    assert_eq!(changed.structure.p1().unwrap().n_atoms(), 8);
    assert!(xtalforge_symmetry::change_setting(&structure, 0).is_err());
}

#[test]
fn bravais_classification_and_band_paths() {
    let fcc =
        xtalforge_core::Lattice::new([[0.0, 2.0, 2.0], [2.0, 0.0, 2.0], [2.0, 2.0, 0.0]]).unwrap();
    assert_eq!(xtalforge_symmetry::classify_bravais(&fcc), "cF");
    let path = xtalforge_symmetry::band_path(&fcc).unwrap();
    assert_eq!(path.bravais, "cF");
    for label in ["G", "X", "W", "K", "L", "U"] {
        assert!(path.special_points.contains_key(label), "missing {label}");
    }
    let points = path.cartesian_points(&fcc).unwrap();
    assert!(points["G"].iter().all(|value| value.abs() < 1e-12));
    assert!(points["X"][0].abs() > 0.0 || points["X"][1].abs() > 0.0);
    let vertices = path.vertices(&fcc).unwrap();
    assert!(vertices.len() >= 4);
    assert!(vertices.last().expect("last").2 > 0.0);

    let cubic = xtalforge_core::Lattice::cubic(4.0).unwrap();
    assert_eq!(xtalforge_symmetry::classify_bravais(&cubic), "cP");
    let hex = xtalforge_core::Lattice::from_parameters(3.0, 3.0, 5.0, 90.0, 90.0, 120.0).unwrap();
    assert_eq!(xtalforge_symmetry::classify_bravais(&hex), "hP");
    let tetragonal =
        xtalforge_core::Lattice::from_parameters(3.0, 3.0, 5.0, 90.0, 90.0, 90.0).unwrap();
    assert_eq!(xtalforge_symmetry::classify_bravais(&tetragonal), "tP");
    let orthorhombic = xtalforge_core::Lattice::orthorhombic(3.0, 4.0, 5.0).unwrap();
    assert_eq!(xtalforge_symmetry::classify_bravais(&orthorhombic), "oP");
}

fn molecule(points: &[(&str, [f64; 3])]) -> Vec<(String, [f64; 3])> {
    points
        .iter()
        .map(|(element, point)| ((*element).to_owned(), *point))
        .collect()
}

#[test]
fn molecular_point_groups() {
    let water = molecule(&[
        ("O", [0.0, 0.0, 0.0]),
        ("H", [0.757, 0.586, 0.0]),
        ("H", [-0.757, 0.586, 0.0]),
    ]);
    assert_eq!(
        xtalforge_symmetry::molecular_point_group(&water, 1e-3)
            .unwrap()
            .schoenflies,
        "C2v"
    );

    let methane = molecule(&[
        ("C", [0.0, 0.0, 0.0]),
        ("H", [1.0, 1.0, 1.0]),
        ("H", [1.0, -1.0, -1.0]),
        ("H", [-1.0, 1.0, -1.0]),
        ("H", [-1.0, -1.0, 1.0]),
    ]);
    assert_eq!(
        xtalforge_symmetry::molecular_point_group(&methane, 1e-3)
            .unwrap()
            .schoenflies,
        "Td"
    );

    let ammonia = molecule(&[
        ("N", [0.0, 0.0, 0.3]),
        ("H", [0.94, 0.0, 0.0]),
        ("H", [-0.47, 0.814, 0.0]),
        ("H", [-0.47, -0.814, 0.0]),
    ]);
    assert_eq!(
        xtalforge_symmetry::molecular_point_group(&ammonia, 1e-3)
            .unwrap()
            .schoenflies,
        "C3v"
    );

    assert!(xtalforge_symmetry::molecular_point_group(&water, 0.0).is_err());
}
