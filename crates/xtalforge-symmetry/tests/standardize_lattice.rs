//! Regression tests for the standardized lattice orientation.
//!
//! moyo stores the basis as columns; the wrapper used to transpose it into a
//! row-major array, which is correct for orthogonal cells but shears every
//! non-orthogonal standardized cell. The quartz cell below is the smallest
//! case that exposes it: the standardized cell must stay hexagonal.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::{Lattice, Site, StructureBuilder};
use xtalforge_symmetry::{
    Reduction, SymmetryOptions, detect_symmetry, reduce_lattice, reduce_structure, standardize,
};

fn quartz() -> xtalforge_core::Structure {
    let lattice = Lattice::new([
        [4.912_39, 0.0, 0.0],
        [-2.456_195, 4.254_254_533_296_639, 0.0],
        [0.0, 0.0, 5.403_85],
    ])
    .expect("hexagonal lattice");
    let positions = [
        [0.4701, 0.0, 0.666_666_666_666_666_6],
        [0.0, 0.4701, 0.333_333_333_333_333_3],
        [0.5299, 0.5299, 0.0],
        [0.4139, 0.2674, 0.7856],
        [0.732_6, 0.146_5, 0.452_266_666_666_666_6],
        [0.853_5, 0.586_1, 0.118_933_333_333_333_3],
        [0.2674, 0.4139, 0.2144],
        [0.146_5, 0.732_6, 0.547_733_333_333_333_4],
        [0.586_1, 0.853_5, 0.881_066_666_666_666_7],
    ];
    let elements = ["Si", "Si", "Si", "O", "O", "O", "O", "O", "O"];
    let mut builder = StructureBuilder::new("quartz", lattice);
    for (index, (position, element)) in positions.iter().zip(elements).enumerate() {
        builder = builder.site(
            Site::new(format!("{element}{}", index + 1), element, *position).expect("quartz site"),
        );
    }
    builder.build().expect("quartz structure")
}

#[test]
fn standardized_hexagonal_cell_keeps_its_basis() {
    let structure = quartz();
    let options = SymmetryOptions {
        symprec: 1e-3,
        ..SymmetryOptions::default()
    };
    let dataset = detect_symmetry(&structure, options).expect("detect quartz");
    assert_eq!(dataset.number, 154);

    let lattice = dataset.std_lattice;
    let a = (lattice[0][0].powi(2) + lattice[0][1].powi(2) + lattice[0][2].powi(2)).sqrt();
    let b = (lattice[1][0].powi(2) + lattice[1][1].powi(2) + lattice[1][2].powi(2)).sqrt();
    let c = (lattice[2][0].powi(2) + lattice[2][1].powi(2) + lattice[2][2].powi(2)).sqrt();
    assert!((a - 4.912_39).abs() < 1e-9, "a = {a}");
    assert!((b - 4.912_39).abs() < 1e-9, "b = {b}");
    assert!((c - 5.403_85).abs() < 1e-9, "c = {c}");
    let cosine = (lattice[0][0] * lattice[1][0]
        + lattice[0][1] * lattice[1][1]
        + lattice[0][2] * lattice[1][2])
        / (a * b);
    assert!((cosine + 0.5).abs() < 1e-9, "gamma cosine = {cosine}");

    let standard = standardize(&structure, options).expect("standardize quartz");
    assert_eq!(standard.structure.sites().len(), 2);
    assert_eq!(standard.structure.p1().expect("P1").n_atoms(), 9);
    assert!((standard.structure.lattice().volume() - 112.932_669_550_927).abs() < 1e-9);
}

#[test]
fn niggli_reduction_keeps_the_hexagonal_basis() {
    let structure = quartz();
    let (reduced, _) = reduce_lattice(structure.lattice(), Reduction::Niggli).expect("Niggli");
    let lattice = reduced.matrix();
    let a = (lattice[0][0].powi(2) + lattice[0][1].powi(2) + lattice[0][2].powi(2)).sqrt();
    let b = (lattice[1][0].powi(2) + lattice[1][1].powi(2) + lattice[1][2].powi(2)).sqrt();
    let c = (lattice[2][0].powi(2) + lattice[2][1].powi(2) + lattice[2][2].powi(2)).sqrt();
    assert!((a - 4.912_39).abs() < 1e-9, "a = {a}");
    assert!((b - 4.912_39).abs() < 1e-9, "b = {b}");
    assert!((c - 5.403_85).abs() < 1e-9, "c = {c}");
    let cosine = (lattice[0][0] * lattice[1][0]
        + lattice[0][1] * lattice[1][1]
        + lattice[0][2] * lattice[1][2])
        / (a * b);
    assert!((cosine + 0.5).abs() < 1e-9, "gamma cosine = {cosine}");

    let reduced_structure = reduce_structure(&structure, Reduction::Niggli).expect("reduce");
    assert_eq!(reduced_structure.structure.p1().expect("P1").n_atoms(), 9);
    assert!((reduced_structure.structure.lattice().volume() - 112.932_669_550_927).abs() < 1e-9);
}
