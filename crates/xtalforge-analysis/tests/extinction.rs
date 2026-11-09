//! Integration tests for the systematic-absence screen (Task 7.10.1).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_analysis::{
    LatticeSystem, PxrdOptions, Radiation, calculate_pxrd, is_systematically_absent,
    space_group_candidates,
};
use xtalforge_core::{Lattice, Site, Structure, StructureBuilder};

/// Reflections of a structure above a fraction of the strongest line.
fn reflections(structure: &Structure, threshold: f64) -> Vec<[i32; 3]> {
    let pattern = calculate_pxrd(
        structure,
        &PxrdOptions {
            radiation: Radiation::CuKa,
            two_theta_range: [10.0, 90.0],
            step: 0.02,
            fwhm: 0.1,
        },
    )
    .unwrap();
    let maximum = pattern
        .reflections
        .iter()
        .map(|reflection| reflection.intensity)
        .fold(0.0_f64, f64::max);
    pattern
        .reflections
        .iter()
        .filter(|reflection| reflection.intensity > threshold * maximum)
        .map(|reflection| reflection.hkl)
        .collect()
}

fn rutile() -> Structure {
    let cif = include_str!("../../../tests/corpus/simple/TiO2.cif");
    let base = xtalforge_io::cif::read_cif(cif).unwrap();
    let mut builder = StructureBuilder::new(
        "rutile",
        Lattice::from_parameters(4.5937, 4.5937, 2.9587, 90.0, 90.0, 90.0).unwrap(),
    )
    .hall(base.hall());
    builder = builder.site(Site::new("Ti", "Ti", [0.0, 0.0, 0.0]).unwrap());
    builder = builder.site(Site::new("O", "O", [0.30478, 0.30478, 0.0]).unwrap());
    for operation in base.operations() {
        builder = builder.operation(*operation);
    }
    builder.build().unwrap()
}

/// Diamond in the conventional cubic cell (eight silicon atoms).
fn silicon() -> Structure {
    let mut builder = StructureBuilder::new("Si", Lattice::cubic(5.430941).unwrap());
    for (index, frac) in [
        [0.0, 0.0, 0.0],
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.5],
        [0.5, 0.5, 0.0],
        [0.25, 0.25, 0.25],
        [0.25, 0.75, 0.75],
        [0.75, 0.25, 0.75],
        [0.75, 0.75, 0.25],
    ]
    .into_iter()
    .enumerate()
    {
        builder = builder.site(Site::new(format!("Si{}", index + 1), "Si", frac).unwrap());
    }
    builder.build().unwrap()
}

#[test]
fn rutile_reflections_admit_its_space_group() {
    let hkl = reflections(&rutile(), 1e-3);
    let candidates = space_group_candidates(LatticeSystem::Tetragonal, &hkl).unwrap();
    println!(
        "rutile {} reflections -> {:?}",
        hkl.len(),
        candidates
            .iter()
            .map(|entry| (entry.number, entry.symbol.as_str(), entry.order))
            .collect::<Vec<_>>()
    );
    assert!(
        candidates.iter().any(|entry| entry.number == 136),
        "P4_2/mnm is not among {candidates:?}"
    );
}

#[test]
fn a_forbidden_reflection_removes_the_group() {
    let mut hkl = reflections(&rutile(), 1e-3);
    // 001 violates the 4_2 screw of P4_2/mnm.
    hkl.push([0, 0, 1]);
    let candidates = space_group_candidates(LatticeSystem::Tetragonal, &hkl).unwrap();
    println!(
        "with 001 -> {:?}",
        candidates
            .iter()
            .map(|entry| (entry.number, entry.symbol.as_str(), entry.order))
            .collect::<Vec<_>>()
    );
    assert!(
        !candidates.iter().any(|entry| entry.number == 136),
        "P4_2/mnm survived an impossible 001"
    );
}

#[test]
fn silicon_reflections_admit_the_diamond_group() {
    let hkl = reflections(&silicon(), 1e-3);
    let candidates = space_group_candidates(LatticeSystem::Cubic, &hkl).unwrap();
    println!(
        "silicon {} reflections -> {:?}",
        hkl.len(),
        candidates
            .iter()
            .map(|entry| (entry.number, entry.symbol.as_str(), entry.order))
            .collect::<Vec<_>>()
    );
    assert!(
        candidates.iter().any(|entry| entry.number == 227),
        "Fd-3m is not among {candidates:?}"
    );
    // 200 is the diamond-glide absence that separates Fd-3m from Fm-3m.
    assert!(hkl.contains(&[1, 1, 1]));
    let with_200 = {
        let mut list = hkl.clone();
        list.push([2, 0, 0]);
        space_group_candidates(LatticeSystem::Cubic, &list).unwrap()
    };
    assert!(
        !with_200.iter().any(|entry| entry.number == 227),
        "Fd-3m survived an impossible 200"
    );
}

#[test]
fn absence_test_matches_known_conditions() {
    let group = xtalforge_symmetry::spacegroup::SpaceGroup::from_number(136).unwrap();
    // Rutile P4_2/mnm: 00l even; h0l and 0kl have even index sums; hk0 has
    // an even h + k; hhl is not restricted.
    assert!(is_systematically_absent(&group.operations, [0, 0, 1]).unwrap());
    assert!(!is_systematically_absent(&group.operations, [0, 0, 2]).unwrap());
    assert!(is_systematically_absent(&group.operations, [1, 0, 0]).unwrap());
    assert!(!is_systematically_absent(&group.operations, [1, 0, 1]).unwrap());
    assert!(is_systematically_absent(&group.operations, [1, 0, 2]).unwrap());
    assert!(is_systematically_absent(&group.operations, [0, 1, 0]).unwrap());
    assert!(!is_systematically_absent(&group.operations, [0, 1, 1]).unwrap());
    assert!(is_systematically_absent(&group.operations, [0, 1, 2]).unwrap());
    assert!(!is_systematically_absent(&group.operations, [1, 1, 0]).unwrap());
    assert!(!is_systematically_absent(&group.operations, [1, 1, 1]).unwrap());
    assert!(!is_systematically_absent(&group.operations, [1, 1, 2]).unwrap());
    // Diamond Fd-3m: the group itself allows 222; its absence in diamond is
    // a special-position extinction, not a general-position one.
    let diamond = xtalforge_symmetry::spacegroup::SpaceGroup::from_number(227).unwrap();
    assert!(is_systematically_absent(&diamond.operations, [2, 0, 0]).unwrap());
    assert!(is_systematically_absent(&diamond.operations, [4, 2, 0]).unwrap());
    assert!(is_systematically_absent(&diamond.operations, [0, 0, 2]).unwrap());
    assert!(is_systematically_absent(&diamond.operations, [0, 1, 1]).unwrap());
    assert!(!is_systematically_absent(&diamond.operations, [1, 1, 1]).unwrap());
    assert!(!is_systematically_absent(&diamond.operations, [2, 2, 0]).unwrap());
    assert!(!is_systematically_absent(&diamond.operations, [2, 2, 2]).unwrap());
    assert!(!is_systematically_absent(&diamond.operations, [4, 0, 0]).unwrap());
    let empty = is_systematically_absent(&[], [1, 0, 0]);
    assert!(empty.is_err());
}
