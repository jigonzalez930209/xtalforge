//! Integration tests for moyo conversion and symmetry detection (Task 2.1).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::Structure;
use xtalforge_io::read_cif;
use xtalforge_symmetry::{AtomType, SpaceGroup, SymmetryOptions, detect_symmetry};

const SI: &str = include_str!("../../../tests/corpus/simple/Si.cif");
const MOF5: &str = include_str!("../../../tests/corpus/cod/MOF-5.cif");

#[test]
fn silicon_is_fd3m_and_matches_the_golden_number() {
    let structure = read_cif(SI).unwrap();
    let dataset = detect_symmetry(&structure, SymmetryOptions::default()).unwrap();
    assert_eq!(dataset.number, 227);
    assert_eq!(dataset.hm_symbol.replace(' ', ""), "Fd-3m");
    assert_eq!(dataset.operations.len(), 192);
    assert_eq!(dataset.wyckoffs.len(), 8);
    assert_eq!(dataset.orbits.len(), 8);
    assert_eq!(dataset.pearson_symbol, "cF8");
}

#[test]
fn mof5_is_fm3m_and_matches_the_golden_number() {
    let structure = read_cif(MOF5).unwrap();
    let dataset = detect_symmetry(&structure, SymmetryOptions::default()).unwrap();
    assert_eq!(dataset.number, 225);
    assert_eq!(dataset.hm_symbol.replace(' ', ""), "Fm-3m");
    assert_eq!(dataset.operations.len(), 192);
    assert_eq!(dataset.wyckoffs.len(), structure.p1().unwrap().n_atoms());
}

#[test]
fn space_group_lookup_accepts_symbols_and_numbers() {
    assert_eq!(
        SpaceGroup::from_number(225)
            .unwrap()
            .hm_short
            .replace(' ', ""),
        "Fm-3m"
    );
    assert_eq!(SpaceGroup::parse("Fm-3m").unwrap().number, 225);
    assert_eq!(SpaceGroup::parse("F m -3 m").unwrap().number, 225);
    assert_eq!(SpaceGroup::parse("225").unwrap().number, 225);
    assert_eq!(SpaceGroup::parse("F d -3 m:1").unwrap().number, 227);
    assert_eq!(SpaceGroup::parse("-F 4 2 3").unwrap().number, 225);
    assert!(SpaceGroup::parse("not a group").is_err());
    assert!(SpaceGroup::from_number(0).is_err());
    assert!(SpaceGroup::from_hall_number(531).is_err());
    assert_eq!(SpaceGroup::from_number(225).unwrap().order(), 192);
    assert!(SpaceGroup::from_number(225).unwrap().is_centrosymmetric());
    assert!(!SpaceGroup::from_number(225).unwrap().is_chiral());
    assert!(SpaceGroup::from_number(1).unwrap().is_chiral());
    assert_eq!(xtalforge_symmetry::all_settings().unwrap().len(), 530);
}

#[test]
fn wyckoff_positions_are_listed() {
    let positions = xtalforge_symmetry::WyckoffPosition::list(526).unwrap();
    assert!(positions.iter().any(|position| position.letter == 'a'));
    assert!(
        positions
            .iter()
            .all(|position| !position.coordinates.is_empty())
    );
    assert!(xtalforge_symmetry::WyckoffPosition::list(0).is_err());
}

#[test]
fn atom_types_include_occupancy_and_dummy_species() {
    let carbon = AtomType::new("C", 1.0).unwrap();
    let half = AtomType::new("C", 0.5).unwrap();
    assert_ne!(carbon, half);
    assert_eq!(AtomType::new("X", 1.0).unwrap().z, 0);
    assert!(AtomType::new("Zz", 1.0).is_err());
    assert!(AtomType::new("C", f64::NAN).is_err());
}

#[test]
fn invalid_options_and_empty_structures_are_refused() {
    let structure = read_cif(SI).unwrap();
    let options = SymmetryOptions {
        symprec: 0.0,
        ..SymmetryOptions::default()
    };
    assert_eq!(
        detect_symmetry(&structure, options).unwrap_err().code(),
        "symmetry"
    );
    let empty = Structure::molecule("empty", &[], 5.0);
    assert!(empty.is_err());
}

/// Golden files whose header group is not recovered at 0.001 A with the
/// documented atom-type policy (element plus occupancy). Each one is listed in
/// the Phase 2 algorithm page with its reason.
const KNOWN_DETECTION_EXCEPTIONS: [&str; 6] = [
    "CFA1.cif",
    "cod/MIL-101.cif",
    "cod/MIL-100.cif",
    "cod/PCN-222.cif",
    "cod/MIL-88B.cif",
    "edge/special_position.cif",
];

#[test]
#[ignore = "slow full corpus sweep; run with --ignored"]
fn corpus_space_group_numbers_match_the_golden_data() {
    let golden_root =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/phase-1");
    let corpus_root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/corpus");
    let mut stack = vec![golden_root];
    let mut checked = 0;
    while let Some(directory) = stack.pop() {
        for entry in std::fs::read_dir(&directory).unwrap() {
            let path = entry.unwrap().path();
            if path.is_dir() {
                stack.push(path);
                continue;
            }
            if path.extension().and_then(|value| value.to_str()) != Some("json") {
                continue;
            }
            let golden: serde_json::Value =
                serde_json::from_str(&std::fs::read_to_string(&path).unwrap()).unwrap();
            let (Some(source), Some(number)) = (
                golden["source"].as_str(),
                golden["space_group"]["number"].as_u64(),
            ) else {
                continue;
            };
            if number == 1 || KNOWN_DETECTION_EXCEPTIONS.contains(&source) {
                continue;
            }
            let text = std::fs::read_to_string(corpus_root.join(source)).unwrap();
            let structure = read_cif(&text).unwrap();
            let symprec = golden["options"]["symprec"].as_f64().unwrap_or(1e-5);
            let dataset = detect_symmetry(
                &structure,
                SymmetryOptions {
                    symprec,
                    ..SymmetryOptions::default()
                },
            )
            .unwrap();
            assert_eq!(
                u64::from(dataset.number),
                number,
                "space group mismatch for {source}"
            );
            checked += 1;
        }
    }
    assert!(
        checked >= 40,
        "expected the Phase 1 corpus, checked {checked}"
    );
}

#[test]
fn small_corpus_space_groups_match_the_golden_data() {
    let golden_root =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/phase-1");
    let corpus_root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/corpus");
    let mut stack = vec![golden_root];
    let mut checked = 0;
    while let Some(directory) = stack.pop() {
        for entry in std::fs::read_dir(&directory).unwrap() {
            let path = entry.unwrap().path();
            if path.is_dir() {
                stack.push(path);
                continue;
            }
            if path.extension().and_then(|value| value.to_str()) != Some("json") {
                continue;
            }
            let golden: serde_json::Value =
                serde_json::from_str(&std::fs::read_to_string(&path).unwrap()).unwrap();
            let (Some(source), Some(number)) = (
                golden["source"].as_str(),
                golden["space_group"]["number"].as_u64(),
            ) else {
                continue;
            };
            if number == 1 || KNOWN_DETECTION_EXCEPTIONS.contains(&source) {
                continue;
            }
            if !source.starts_with("simple/") {
                continue;
            }
            let text = std::fs::read_to_string(corpus_root.join(source)).unwrap();
            let structure = read_cif(&text).unwrap();
            let symprec = golden["options"]["symprec"].as_f64().unwrap_or(1e-5);
            let dataset = detect_symmetry(
                &structure,
                SymmetryOptions {
                    symprec,
                    ..SymmetryOptions::default()
                },
            )
            .unwrap();
            assert_eq!(
                u64::from(dataset.number),
                number,
                "space group mismatch for {source}"
            );
            checked += 1;
        }
    }
    assert!(checked >= 10, "expected a small corpus, checked {checked}");
}

#[test]
fn space_group_value_object_has_the_full_identity() {
    let group = SpaceGroup::from_number(225).unwrap();
    assert_eq!(group.crystal_system, "cubic");
    assert_eq!(group.lattice_system, "cubic");
    assert_eq!(group.centering, 'F');
    assert_eq!(group.bravais, "cF");
    assert_eq!(group.point_group, "m-3m");
    assert_eq!(group.schoenflies, "Oh");
    assert_eq!(group.laue_class, "m-3m");
    assert!(group.is_centrosymmetric());
    assert!(!group.is_chiral());
    let trigonal = SpaceGroup::from_number(150).unwrap();
    assert_eq!(trigonal.point_group, "32");
    assert_eq!(trigonal.laue_class, "-3m");
}

#[test]
fn catalogue_filters_by_system_centering_and_point_group() {
    let cubic = xtalforge_symmetry::catalogue(xtalforge_symmetry::SpaceGroupFilter {
        crystal_system: Some("cubic"),
        ..Default::default()
    })
    .unwrap();
    assert!(cubic.len() >= 36);
    assert!(cubic.iter().all(|group| group.crystal_system == "cubic"));
    let face_centered = xtalforge_symmetry::catalogue(xtalforge_symmetry::SpaceGroupFilter {
        centering: Some('F'),
        ..Default::default()
    })
    .unwrap();
    assert!(face_centered.iter().all(|group| group.centering == 'F'));
    assert!(face_centered.iter().any(|group| group.number == 225));
    let cubic_m3m = xtalforge_symmetry::catalogue(xtalforge_symmetry::SpaceGroupFilter {
        point_group: Some("m-3m"),
        ..Default::default()
    })
    .unwrap();
    assert!(!cubic_m3m.is_empty());
    assert!(cubic_m3m.iter().all(|group| group.point_group == "m-3m"));
}

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
            .all(|assignment| assignment.letter == 'b')
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
