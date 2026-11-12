//! Integration tests for strain, rattle, editing and optimal supercells
//! (Task 9.2).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::collections::BTreeMap;

use xtalforge_core::{
    Lattice, RattleOptions, Site, SiteSort, SpeciesReplacement, Structure, StructureBuilder,
    add_sites, deform, optimal_supercell, rattle, remove_sites, replace_site, replace_species,
    scale_to_density, scale_to_volume, sort_sites, strain, voigt_to_tensor, wrap_sites,
};

fn diamond() -> Structure {
    StructureBuilder::new("Si", Lattice::cubic(5.430941).unwrap())
        .site(Site::new("Si1", "Si", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("Si2", "Si", [0.25, 0.25, 0.25]).unwrap())
        .build()
        .unwrap()
}

#[test]
fn isotropic_strain_scales_the_volume() {
    let structure = diamond();
    let strained = strain(
        &structure,
        voigt_to_tensor([0.01, 0.01, 0.01, 0.0, 0.0, 0.0]).unwrap(),
    )
    .unwrap();
    let ratio = strained.lattice().volume() / structure.lattice().volume();
    assert!((ratio - 1.030_301).abs() < 1e-4);
}

#[test]
fn deformation_gradient_and_volume_helpers_agree() {
    let structure = diamond();
    let scaled = scale_to_volume(&structure, 2.0 * structure.lattice().volume()).unwrap();
    let deform_gradient = deform(
        &structure,
        [
            [2.0_f64.cbrt(), 0.0, 0.0],
            [0.0, 2.0_f64.cbrt(), 0.0],
            [0.0, 0.0, 2.0_f64.cbrt()],
        ],
    )
    .unwrap();
    assert!((scaled.lattice().volume() - deform_gradient.lattice().volume()).abs() < 1e-9);
}

#[test]
fn density_scaling_uses_the_cell_mass() {
    let structure = diamond();
    let target = 3.0;
    let scaled = scale_to_density(&structure, target).unwrap();
    assert!((scaled.density().unwrap() - target).abs() < 1e-6);
}

#[test]
fn rattle_is_deterministic_for_a_seed() {
    let options = RattleOptions {
        amplitude: 0.1,
        seed: 42,
        ..RattleOptions::default()
    };
    let first = rattle(&diamond(), &options).unwrap();
    let second = rattle(&diamond(), &options).unwrap();
    assert_eq!(first.structure, second.structure);
    assert!(first.max_displacement > 0.0);
    let other = rattle(
        &diamond(),
        &RattleOptions {
            seed: 43,
            ..options
        },
    )
    .unwrap();
    assert_ne!(first.structure, other.structure);
}

#[test]
fn species_replacement_can_be_partial_and_seeded() {
    let structure = StructureBuilder::new("alloy", Lattice::cubic(4.0).unwrap())
        .site(Site::new("Si1", "Si", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("Si2", "Si", [0.5, 0.5, 0.5]).unwrap())
        .build()
        .unwrap();
    let mut map = BTreeMap::new();
    map.insert("Si".to_owned(), "Ge".to_owned());
    let replaced = replace_species(
        &structure,
        &SpeciesReplacement {
            map,
            fraction: 0.5,
            seed: 7,
        },
    )
    .unwrap();
    let germanium = replaced
        .sites()
        .iter()
        .filter(|site| site.element() == "Ge")
        .count();
    assert_eq!(germanium, 1);
}

#[test]
fn sorting_and_wrapping_keep_the_content() {
    let structure = StructureBuilder::new("mix", Lattice::cubic(4.0).unwrap())
        .site(Site::new("O1", "O", [1.2, -0.2, 0.0]).unwrap())
        .site(Site::new("Si1", "Si", [0.0, 0.0, 0.0]).unwrap())
        .build()
        .unwrap();
    let sorted = sort_sites(&structure, SiteSort::Element).unwrap();
    assert_eq!(sorted.sites()[0].element(), "O");
    assert_eq!(sorted.sites()[1].element(), "Si");
    let wrapped = wrap_sites(&structure).unwrap();
    for site in wrapped.sites() {
        for value in site.frac() {
            assert!((0.0..1.0).contains(&value));
        }
    }
}

#[test]
fn site_editing_adds_removes_and_replaces() {
    let structure = diamond();
    let added = add_sites(&structure, &[("Ge".to_owned(), [0.5, 0.5, 0.5])]).unwrap();
    assert_eq!(added.sites().len(), 3);
    let removed = remove_sites(&added, &[2]).unwrap();
    assert_eq!(removed.sites().len(), 2);
    let replaced = replace_site(&removed, 0, "Ge").unwrap();
    assert_eq!(replaced.sites()[0].element(), "Ge");
}

#[test]
fn optimal_supercell_approaches_a_cube() {
    let choice = optimal_supercell(&diamond(), 8).unwrap();
    assert_eq!(choice.cells, 8);
    assert_eq!(choice.structure.sites().len(), 16);
    assert!(choice.score >= 0.0);
    let determinant = choice.matrix[0][0]
        * (choice.matrix[1][1] * choice.matrix[2][2] - choice.matrix[1][2] * choice.matrix[2][1])
        - choice.matrix[0][1]
            * (choice.matrix[1][0] * choice.matrix[2][2]
                - choice.matrix[1][2] * choice.matrix[2][0])
        + choice.matrix[0][2]
            * (choice.matrix[1][0] * choice.matrix[2][1]
                - choice.matrix[1][1] * choice.matrix[2][0]);
    assert!((determinant.abs() - 8.0).abs() < 1e-9);
}
