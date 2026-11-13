//! Integration tests for the plane-wave and electronic-structure formats:
//! CASTEP, FHI-aims, CP2K, ABINIT and CRYSTAL (Task 13.2.2).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use serde_json::Value;
use xtalforge_core::{Lattice, Site, StructureBuilder};
use xtalforge_io::registry::{find, read_with, sniff, write_with};
use xtalforge_io::{
    read_abinit, read_aims, read_castep, read_cp2k, read_crystal, write_abinit, write_aims,
    write_castep, write_cp2k, write_crystal,
};

const CASTEP: &str = "\
%BLOCK LATTICE_CART
  5.43 0.0 0.0
  0.0 5.43 0.0
  0.0 0.0 5.43
%ENDBLOCK LATTICE_CART

%BLOCK POSITIONS_FRAC
  Si 0.0 0.0 0.0
  Si 0.5 0.5 0.5
%ENDBLOCK POSITIONS_FRAC
";

const AIMS: &str = "\
lattice_vector 5.43 0.0 0.0
lattice_vector 0.0 5.43 0.0
lattice_vector 0.0 0.0 5.43
atom_frac 0.0 0.0 0.0 Si
atom_frac 0.25 0.25 0.25 Si
constrain_relaxation .true.
initial_charge 0.5
";

const CP2K: &str = "\
&GLOBAL
  PROJECT cp2k-si
&END GLOBAL
&FORCE_EVAL
  &SUBSYS
    &CELL
      ABC 5.43 5.43 5.43
    &END CELL
    &COORD
      SCALED
      Si 0.0 0.0 0.0
      Si 0.25 0.25 0.25
    &END COORD
    &KIND Si
      ELEMENT Si
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
";

const ABINIT: &str = "\
natom 2
ntypat 1
typat 1 1
znucl 14
acell 5.43 5.43 5.43 angstrom
rprim 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0
xred 0.0 0.0 0.0 0.25 0.25 0.25
";

const CRYSTAL: &str = "\
Silicon
CRYSTAL
0 0 0
227
5.43 5.43 5.43 90.0 90.0 90.0
2
14 0.0 0.0 0.0
14 0.25 0.25 0.25
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
fn castep_round_trip_keeps_the_cell() {
    let structure = read_castep(CASTEP).unwrap();
    assert_eq!(structure.formula(), "Si2");
    assert!((structure.lattice().volume() - 160.103_007).abs() < 1e-4);
    let text = write_castep(&structure).unwrap();
    assert!(text.contains("%BLOCK LATTICE_CART"));
    assert!(text.contains("%BLOCK POSITIONS_FRAC"));
    assert!(text.contains("%BLOCK SPECIES_POT"));
    let back = read_castep(&text).unwrap();
    assert_eq!(back.formula(), "Si2");
    assert!((back.lattice().volume() - 160.103_007).abs() < 1e-4);
}

#[test]
fn castep_reads_abc_cells_and_rejects_missing_blocks() {
    let abc = "\
%BLOCK LATTICE_ABC
  4.0 5.0 6.0 90.0 90.0 90.0
%ENDBLOCK LATTICE_ABC
%BLOCK POSITIONS_CART
  C 0.0 0.0 0.0
%ENDBLOCK POSITIONS_CART
";
    let structure = read_castep(abc).unwrap();
    assert_eq!(structure.formula(), "C");
    assert!((structure.lattice().volume() - 120.0).abs() < 1e-9);
    assert!(read_castep("%BLOCK POSITIONS_FRAC\n C 0 0 0\n%ENDBLOCK POSITIONS_FRAC").is_err());
    assert!(read_castep("%BLOCK LATTICE_CART\n5 0 0\n0 5 0\n0 0 5\n%ENDBLOCK LATTICE_CART\n%BLOCK POSITIONS_FRAC\nSi 0 0 0\n").is_err());
}

#[test]
fn fhi_aims_round_trip_keeps_constraints_and_charges() {
    let structure = read_aims(AIMS).unwrap();
    assert_eq!(structure.formula(), "Si2");
    let last = &structure.sites()[1];
    assert_eq!(
        last.props().get("constrain_relaxation"),
        Some(&Value::Bool(true))
    );
    assert_eq!(last.props().get("initial_charge"), Some(&Value::from(0.5)));
    let text = write_aims(&structure).unwrap();
    assert!(text.contains("lattice_vector"));
    assert!(text.contains("atom_frac"));
    assert!(text.contains("constrain_relaxation .true."));
    let back = read_aims(&text).unwrap();
    assert_eq!(back.sites().len(), 2);
    assert!((back.lattice().volume() - 160.103_007).abs() < 1e-4);
    assert_eq!(
        back.sites()[1].props().get("constrain_relaxation"),
        Some(&Value::Bool(true))
    );
}

#[test]
fn fhi_aims_without_a_lattice_becomes_a_molecule() {
    let text = "atom 0.0 0.0 0.0 O\natom 1.0 0.0 0.0 H\n";
    let structure = read_aims(text).unwrap();
    assert_eq!(structure.sites().len(), 2);
    let mut elements: Vec<&str> = structure
        .sites()
        .iter()
        .map(xtalforge_core::Site::element)
        .collect();
    elements.sort_unstable();
    assert_eq!(elements, vec!["H", "O"]);
    assert!(structure.lattice().volume() > 0.0);
    assert!(read_aims("atom_frac 0.0 0.0 0.0 Si\n").is_err());
}

#[test]
fn cp2k_reads_scaled_coordinates_and_round_trips() {
    let structure = read_cp2k(CP2K).unwrap();
    assert_eq!(structure.formula(), "Si2");
    assert_eq!(structure.name(), "cp2k-si");
    assert!((structure.lattice().volume() - 160.103_007).abs() < 1e-4);
    let text = write_cp2k(&structure).unwrap();
    assert!(text.contains("&FORCE_EVAL"));
    assert!(text.contains("ABC"));
    assert!(text.contains("SCALED"));
    let back = read_cp2k(&text).unwrap();
    assert_eq!(back.formula(), "Si2");
    assert!((back.sites()[1].frac()[0] - 0.25).abs() < 1e-9);
}

#[test]
fn cp2k_rejects_unclosed_sections_and_bad_vectors() {
    let unclosed = "&FORCE_EVAL\n &SUBSYS\n &END SUBSYS\n";
    assert!(read_cp2k(unclosed).is_err());
    let bad_cell = CP2K.replace("ABC 5.43 5.43 5.43", "ABC 5.43 nope 5.43");
    assert!(read_cp2k(&bad_cell).is_err());
    assert!(read_cp2k("&GLOBAL\n&END GLOBAL\n").is_err());
}

#[test]
fn abinit_round_trip_keeps_cell_and_types() {
    let structure = read_abinit(ABINIT).unwrap();
    assert_eq!(structure.formula(), "Si2");
    assert!((structure.lattice().volume() - 160.103_007).abs() < 1e-4);
    let text = write_abinit(&structure).unwrap();
    assert!(text.contains("acell 1.0 1.0 1.0 angstrom"));
    assert!(text.contains("rprim"));
    assert!(text.contains("xred"));
    let back = read_abinit(&text).unwrap();
    assert_eq!(back.formula(), "Si2");
    assert!((back.lattice().volume() - 160.103_007).abs() < 1e-4);
    assert!((back.sites()[1].frac()[0] - 0.25).abs() < 1e-9);
}

#[test]
fn abinit_supports_repetition_and_xcart() {
    let text = "\
natom 2
ntypat 2
typat 1 2
znucl 6 8
acell 3*4.0 angstrom
rprim 1 0 0 0 1 0 0 0 1
xcart 0.0 0.0 0.0 1.2 0.0 0.0
";
    let structure = read_abinit(text).unwrap();
    assert_eq!(structure.formula(), "C O");
    assert!((structure.lattice().lengths()[0] - 4.0).abs() < 1e-9);
    assert!(read_abinit("natom 1\nntypat 1\ntypat 1\nznucl 6\n").is_err());
}

#[test]
fn crystal_reads_space_groups_and_round_trips() {
    let structure = read_crystal(CRYSTAL).unwrap();
    assert_eq!(structure.sites().len(), 2);
    assert_eq!(structure.operations().len(), 192);
    assert_eq!(structure.p1().unwrap().n_atoms(), 16);
    let text = write_crystal(&structure).unwrap();
    assert!(text.contains("CRYSTAL"));
    assert!(text.contains("227"));
    let back = read_crystal(&text).unwrap();
    assert_eq!(back.sites().len(), 2);
    assert_eq!(back.operations().len(), 192);
    assert_eq!(back.p1().unwrap().n_atoms(), 16);
}

#[test]
fn crystal_rejects_slabs_and_bad_atomic_numbers() {
    let slab = CRYSTAL.replace("0 0 0", "1 0 0");
    assert!(read_crystal(&slab).is_err());
    let bad_z = CRYSTAL.replace("14 0.0 0.0 0.0", "200 0.0 0.0 0.0");
    assert!(read_crystal(&bad_z).is_err());
    assert!(read_crystal("not a crystal").is_err());
}

#[test]
fn the_new_formats_are_registered_and_sniffed() {
    for (extension, name) in [
        (".cell", "castep"),
        (".geometry.in", "fhi-aims"),
        (".aims", "fhi-aims"),
        (".cp2k", "cp2k"),
        (".inp", "cp2k"),
        (".abi", "abinit"),
        (".crystal", "crystal"),
    ] {
        let format = find(extension).unwrap();
        assert_eq!(format.name, name, "{extension}");
        assert!(format.read.is_some(), "{extension}");
        assert!(format.write.is_some(), "{extension}");
    }
    assert_eq!(sniff(CASTEP).unwrap().name, "castep");
    assert_eq!(sniff(AIMS).unwrap().name, "fhi-aims");
    assert_eq!(sniff(CP2K).unwrap().name, "cp2k");
    assert_eq!(sniff(ABINIT).unwrap().name, "abinit");
    assert_eq!(read_with(CASTEP, "castep").unwrap().formula(), "Si2");
}

#[test]
fn the_formats_report_round_trip_loss_consistently() {
    let castep = silicon();
    let text = write_with(&castep, "castep");
    assert!(text.is_ok());
    let back = read_with(&text.unwrap(), ".cell").unwrap();
    assert_eq!(back.sites().len(), 2);
}
