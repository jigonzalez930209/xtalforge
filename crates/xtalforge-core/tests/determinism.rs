//! Determinism across repeated and concurrent runs (Task 8.6.2).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::sync::Arc;

use xtalforge_core::{Lattice, Site, Structure, StructureBuilder, neighbors};

/// A 3x3x3 grid of carbon atoms in a 6 A cell.
fn carbon_cell() -> Structure {
    let mut builder = StructureBuilder::new("grid", Lattice::cubic(6.0).unwrap());
    let mut index = 0;
    for x in 0..3 {
        for y in 0..3 {
            for z in 0..3 {
                let frac = [f64::from(x) / 3.0, f64::from(y) / 3.0, f64::from(z) / 3.0];
                builder = builder.site(Site::new(format!("C{index}"), "C", frac).unwrap());
                index += 1;
            }
        }
    }
    builder.build().unwrap()
}

fn structure() -> Arc<Structure> {
    Arc::new(carbon_cell().supercell(2, 2, 2).unwrap())
}

fn signature(structure: &Structure) -> String {
    let graph = structure.bond_graph().unwrap();
    let neighbors = neighbors(structure, 3.0).unwrap();
    let mut text = format!(
        "{} sites={} neighbors={}\n",
        structure.summary().unwrap().formula,
        structure.sites().len(),
        neighbors.len()
    );
    for bond in &graph.bonds {
        text.push_str(&format!(
            "{} {} {:?} {:.12} {:.12}\n",
            bond.i, bond.j, bond.image, bond.distance, bond.order
        ));
    }
    text
}

#[test]
fn repeated_runs_give_identical_results() {
    let structure = structure();
    let first = signature(&structure);
    for _ in 0..4 {
        assert_eq!(signature(&structure), first);
    }
}

#[test]
fn concurrent_runs_give_identical_results() {
    let structure = structure();
    let expected = signature(&structure);
    let mut handles = Vec::new();
    for _ in 0..4 {
        let structure = Arc::clone(&structure);
        handles.push(std::thread::spawn(move || signature(&structure)));
    }
    for handle in handles {
        assert_eq!(handle.join().unwrap(), expected);
    }
}

#[test]
fn a_supercell_is_the_same_whatever_the_thread_count() {
    let structure = carbon_cell();
    let first = structure.supercell(3, 2, 1).unwrap().summary().unwrap();
    let second = structure.supercell(3, 2, 1).unwrap().summary().unwrap();
    assert_eq!(first, second);
    assert_eq!(
        structure.supercell(3, 2, 1).unwrap().sites().len(),
        structure.sites().len() * 6
    );
}
