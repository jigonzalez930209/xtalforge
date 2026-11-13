//! Integration tests for mmCIF/PDBx, the CIF Cartesian path and the PDB
//! reader (Task 13.2.3).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::{Lattice, Site, StructureBuilder};
use xtalforge_io::registry::{find, read_with, sniff};
use xtalforge_io::{read_cif, read_mmcif, read_pdb, write_pdb};

const MMCIF: &str = "\
data_1ABC
_cell.length_a 10.0
_cell.length_b 10.0
_cell.length_c 10.0
_cell.angle_alpha 90.0
_cell.angle_beta 90.0
_cell.angle_gamma 90.0
_symmetry.space_group_name_H-M 'P 1'
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_asym_id
_atom_site.label_comp_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_PDB_model_num
ATOM 1 C C1 A GLY 1 1.0 2.0 3.0 1.0 8.0 1
ATOM 2 O O1 A GLY 1 4.0 5.0 6.0 0.5 16.0 1
ATOM 3 C C2 B GLY 2 7.0 8.0 9.0 1.0 0.0 2
";

const CARTESIAN_CIF: &str = "\
data_cart
_cell_length_a 10.0
_cell_length_b 10.0
_cell_length_c 10.0
_cell_angle_alpha 90.0
_cell_angle_beta 90.0
_cell_angle_gamma 90.0
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_Cartn_x
_atom_site_Cartn_y
_atom_site_Cartn_z
C1 C 1.0 2.0 3.0
";

const PDB: &str = "\
HEADER    TEST STRUCTURE
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1     1
ATOM      1  C1  GLY A   1       1.000   2.000   3.000  1.00  8.00           C
ATOM      2  O1  GLY A   1       4.000   5.000   6.000  0.50 16.00           O
CONECT    1    2
END
";

fn silicon() -> xtalforge_core::Structure {
    StructureBuilder::new("silicon", Lattice::cubic(5.43).unwrap())
        .site(Site::new("Si1", "Si", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("Si2", "Si", [0.25, 0.25, 0.25]).unwrap())
        .build()
        .unwrap()
}

#[test]
fn mmcif_reads_only_the_first_model_and_keeps_hierarchy() {
    let structure = read_mmcif(MMCIF).unwrap();
    assert_eq!(structure.sites().len(), 2);
    assert_eq!(structure.sites().len(), 2);
    assert_eq!(structure.sites()[0].element(), "C");
    assert_eq!(structure.sites()[1].element(), "O");
    assert_eq!(structure.sites()[0].label(), "C1");
    assert_eq!(structure.sites()[0].element(), "C");
    assert_eq!(structure.sites()[0].props().get("chain_id").unwrap(), "A");
    assert_eq!(
        structure.sites()[0].props().get("residue_name").unwrap(),
        "GLY"
    );
    assert_eq!(
        structure.sites()[0].props().get("residue_number").unwrap(),
        "1"
    );
    let [x, y, z] = structure.sites()[0].frac();
    assert!((x - 0.1).abs() < 1e-9);
    assert!((y - 0.2).abs() < 1e-9);
    assert!((z - 0.3).abs() < 1e-9);
    assert!((structure.sites()[1].occupancy() - 0.5).abs() < 1e-9);
    let u_iso = structure.sites()[0].u_iso().unwrap();
    assert!((u_iso - 8.0 / (8.0 * std::f64::consts::PI * std::f64::consts::PI)).abs() < 1e-9);
}

#[test]
fn cif_accepts_cartesian_atom_columns() {
    let structure = read_cif(CARTESIAN_CIF).unwrap();
    assert_eq!(structure.formula(), "C");
    let [x, y, z] = structure.sites()[0].frac();
    assert!((x - 0.1).abs() < 1e-9);
    assert!((y - 0.2).abs() < 1e-9);
    assert!((z - 0.3).abs() < 1e-9);
    assert_eq!(sniff(MMCIF).unwrap().name, "mmcif");
    assert_eq!(sniff(CARTESIAN_CIF).unwrap().name, "cif");
}

#[test]
fn mmcif_rejects_a_missing_atom_loop_and_missing_cell() {
    assert!(read_mmcif("data_x\n_cell.length_a 10\n").is_err());
    let no_cell = MMCIF.replace("_cell.length_a 10.0\n", "");
    assert!(read_mmcif(&no_cell).is_err());
    assert!(read_mmcif("").is_err());
}

#[test]
fn pdb_reads_cell_symmetry_occupancy_and_bonds() {
    let structure = read_pdb(PDB).unwrap();
    assert_eq!(structure.sites().len(), 2);
    assert_eq!(structure.sites()[0].element(), "C");
    assert!((structure.lattice().lengths()[0] - 10.0).abs() < 1e-6);
    assert!((structure.sites()[1].occupancy() - 0.5).abs() < 1e-9);
    assert!((structure.sites()[0].u_iso().unwrap() - 0.101_321).abs() < 1e-4);
    assert_eq!(structure.sites()[0].props().get("chain_id").unwrap(), "A");
    assert_eq!(structure.explicit_bonds().len(), 1);
    assert!(structure.is_p1());
}

#[test]
fn pdb_without_a_cell_becomes_a_molecule() {
    let text = "\
ATOM      1  O1  HOH A   1       0.000   0.000   0.000  1.00  0.00           O
ATOM      2  H1  HOH A   1       1.000   0.000   0.000  1.00  0.00           H
END
";
    let structure = read_pdb(text).unwrap();
    assert_eq!(structure.sites().len(), 2);
    let elements: Vec<&str> = structure
        .sites()
        .iter()
        .map(xtalforge_core::Site::element)
        .collect();
    assert!(elements.contains(&"H") && elements.contains(&"O"));
    assert!(structure.lattice().volume() > 0.0);
}

#[test]
fn pdb_writer_preserves_the_cell_on_a_round_trip() {
    let text = write_pdb(&silicon()).unwrap();
    assert!(text.starts_with("CRYST1"));
    assert!(text.ends_with("END\n"));
    let back = read_pdb(&text).unwrap();
    assert_eq!(back.formula(), "Si2");
    assert!((back.lattice().lengths()[0] - 5.43).abs() < 1e-3);
    assert!((back.lattice().lengths()[1] - 5.43).abs() < 1e-3);
}

#[test]
fn pdb_rejects_malformed_atom_records() {
    let short = "\
ATOM      1  C1  GLY A   1       1.000   nope    3.000  1.00  0.00           C
";
    assert!(read_pdb(short).is_err());
    let bad_cell = "CRYST1   10.000   10.000\n";
    assert!(read_pdb(bad_cell).is_err());
    assert!(read_pdb("").is_err());
}

#[test]
fn mmcif_and_pdb_are_registered() {
    let mmcif = find(".mmcif").unwrap();
    assert_eq!(mmcif.name, "mmcif");
    assert!(mmcif.read.is_some());
    assert!(mmcif.write.is_none());
    assert!(mmcif.capabilities.cell && mmcif.capabilities.labels);
    assert_eq!(find(".pdbx").unwrap().name, "mmcif");
    assert_eq!(find(".pdb").unwrap().name, "pdb");
    assert!(find(".pdb").unwrap().read.is_some());
    assert_eq!(read_with(PDB, "pdb").unwrap().sites().len(), 2);
    assert_eq!(read_with(MMCIF, "mmcif").unwrap().sites().len(), 2);
}
