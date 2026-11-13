//! Integration tests for the structure-of-arrays views (Task 15.2.1).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::{
    Lattice, NeighborArrays, P1Arrays, Site, StructureBuilder, SymOp, neighbors, neighbors_arrays,
};

fn structure() -> xtalforge_core::Structure {
    StructureBuilder::new("soa", Lattice::cubic(5.0).unwrap())
        .hall("P -1")
        .operation(SymOp::identity())
        .operation(SymOp::parse("-x,-y,-z").unwrap())
        .site(Site::new("C1", "C", [0.1, 0.2, 0.3]).unwrap())
        .site(Site::new("O1", "O", [0.0, 0.0, 0.0]).unwrap())
        .build()
        .unwrap()
}

#[test]
fn p1_columns_match_the_cell() {
    let structure = structure();
    let cell = structure.p1().unwrap();
    let view = P1Arrays::from_structure_default(&structure).unwrap();
    assert_eq!(view.len(), cell.n_atoms());
    assert_eq!(view.lattice(), &cell.lattice);
    for (index, atom) in cell.atoms.iter().enumerate() {
        assert_eq!(view.elements()[index], atom.element);
        assert_eq!(view.labels()[index], atom.label);
        assert_eq!(view.frac(index), Some(atom.frac));
        assert_eq!(view.occupancy()[index], atom.occupancy);
        assert_eq!(view.site_indices()[index], atom.site_idx);
        assert_eq!(view.operation_indices()[index], atom.op_idx);
        assert_eq!(view.taus()[index], atom.tau);
        assert_eq!(view.disorder_groups()[index], atom.disorder_group);
    }
    let (x, y, z) = view.columns();
    assert_eq!(
        (x.len(), y.len(), z.len()),
        (cell.n_atoms(), cell.n_atoms(), cell.n_atoms())
    );
}

#[test]
fn p1_view_round_trips() {
    let cell = structure().p1().unwrap();
    let view = P1Arrays::from_p1(&cell);
    assert_eq!(view.to_p1(), cell);
    assert_eq!(view.frac(view.len()), None);
}

#[test]
fn explicit_tolerance_is_validated() {
    let structure = structure();
    assert!(P1Arrays::from_structure(&structure, -0.1).is_err());
    assert!(P1Arrays::from_structure(&structure, f64::NAN).is_err());
}

#[test]
fn neighbor_columns_match_the_neighbor_list() {
    let structure = structure();
    let found = neighbors(&structure, 1.2).unwrap();
    let view = NeighborArrays::new(&found, structure.p1().unwrap().n_atoms()).unwrap();
    assert_eq!(view.len(), found.len());
    assert_eq!(view.iter().collect::<Vec<_>>(), found);
    for (index, neighbor) in found.iter().enumerate() {
        assert_eq!(view.get(index), Some(*neighbor));
        assert_eq!(view.i()[index], neighbor.i);
        assert_eq!(view.j()[index], neighbor.j);
        assert_eq!(view.distances()[index], neighbor.distance);
    }
    let (x, y, z) = view.images();
    assert!(
        found
            .iter()
            .enumerate()
            .all(|(index, neighbor)| { [x[index], y[index], z[index]] == neighbor.image })
    );
}

#[test]
fn rows_group_every_neighbor_of_one_atom() {
    let structure = structure();
    let found = neighbors(&structure, 2.4).unwrap();
    let n_atoms = structure.p1().unwrap().n_atoms();
    let view = NeighborArrays::new(&found, n_atoms).unwrap();
    for atom in 0..n_atoms {
        let range = view.row_range(atom).unwrap();
        assert!(
            (range.start..range.end).all(|position| view.i()[position] == atom),
            "row {atom} is not grouped"
        );
        let expected = found.iter().filter(|neighbor| neighbor.i == atom).count();
        assert_eq!(range.len(), expected);
    }
    assert_eq!(view.row_range(n_atoms), None);
}

#[test]
fn out_of_range_neighbors_are_rejected() {
    let found = [xtalforge_core::Neighbor {
        i: 2,
        j: 0,
        image: [0, 0, 0],
        distance: 1.0,
    }];
    let error = NeighborArrays::new(&found, 1).expect_err("an index outside the cell must fail");
    assert!(error.message().contains("outside a cell"));
}

#[test]
fn the_shortcut_matches_the_direct_path() {
    let structure = structure();
    let direct = neighbors_arrays(&structure, 1.2).unwrap();
    let view = P1Arrays::from_structure_default(&structure).unwrap();
    assert_eq!(view.neighbor_arrays(1.2).unwrap(), direct);
    assert_eq!(
        direct.iter().collect::<Vec<_>>(),
        neighbors(&structure, 1.2).unwrap()
    );
}

#[test]
fn empty_cells_are_supported() {
    let empty = StructureBuilder::new("empty", Lattice::cubic(3.0).unwrap())
        .build()
        .unwrap();
    let view = P1Arrays::from_structure_default(&empty).unwrap();
    assert!(view.is_empty());
    assert_eq!(view.to_p1().atoms.len(), 0);
    assert!(view.neighbor_arrays(1.0).unwrap().is_empty());
    assert!(NeighborArrays::new(&[], 0).unwrap().is_empty());
}
