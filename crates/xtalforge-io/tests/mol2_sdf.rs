//! Integration tests for the MOL2 and SDF/MOL formats (Task 13.2.5).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_io::registry::{find, read_with, sniff};
use xtalforge_io::{read_mol2, read_sdf, read_sdf_all, write_mol2, write_sdf};

const MOL2: &str = "\
@<TRIPOS>MOLECULE
water
3 2 1 0 0
SMALL
USER_CHARGES

@<TRIPOS>ATOM
1 O1 0.0000 0.0000 0.0000 O.3 1 MOL -0.8000
2 H1 1.0000 0.0000 0.0000 H 1 MOL 0.4000
3 H2 0.0000 1.0000 0.0000 H 1 MOL 0.4000
@<TRIPOS>BOND
1 1 2 1
2 1 3 1
";

const SDF: &str = "\
water
  xtalforge
comment
  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0
    0.0000    1.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
M  CHG  3   1  -1   2   1   3   1
M  END
>  <ID>
water

$$$$
";

const SDF_V3000: &str = "\
water
  xtalforge
comment
  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0 0.0 0.0 0
M  V30 2 O 1.2 0.0 0.0 0 CHG=-1
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 1 2
M  V30 END BOND
M  V30 END CTAB
M  END
$$$$
";

#[test]
fn mol2_reads_atoms_charges_bonds_and_metadata() {
    let structure = read_mol2(MOL2).unwrap();
    assert_eq!(structure.sites().len(), 3);
    assert_eq!(structure.sites()[0].element(), "O");
    assert!((structure.sites()[0].charge().unwrap() + 0.8).abs() < 1e-9);
    assert!((structure.sites()[1].charge().unwrap() - 0.4).abs() < 1e-9);
    assert_eq!(structure.explicit_bonds().len(), 2);
    assert_eq!(
        structure.metadata().get("mol2_charge_type"),
        Some(&serde_json::Value::String("USER_CHARGES".to_owned()))
    );
    assert_eq!(
        structure
            .sites()
            .iter()
            .map(xtalforge_core::Site::element)
            .filter(|element| *element == "H")
            .count(),
        2
    );
}

#[test]
fn mol2_round_trip_keeps_charges_and_bonds() {
    let structure = read_mol2(MOL2).unwrap();
    let text = write_mol2(&structure).unwrap();
    assert!(text.starts_with("@<TRIPOS>MOLECULE"));
    assert!(text.contains("@<TRIPOS>ATOM"));
    assert!(text.contains("@<TRIPOS>BOND"));
    let back = read_mol2(&text).unwrap();
    assert_eq!(back.sites().len(), 3);
    assert!((back.sites()[0].charge().unwrap() + 0.8).abs() < 1e-3);
    assert_eq!(back.explicit_bonds().len(), 2);
}

#[test]
fn mol2_reports_malformed_sections_and_atom_types() {
    assert!(read_mol2("").is_err());
    assert!(read_mol2("@<TRIPOS>MOLECULE\nx\n").is_err());
    let bad_type = MOL2.replace("O.3", "Q");
    assert!(read_mol2(&bad_type).is_err());
    let bad_bond = MOL2.replace("1 1 2 1", "1 1 9 1");
    assert!(read_mol2(&bad_bond).is_err());
}

#[test]
fn sdf_v2000_reads_charges_bonds_and_data_fields() {
    let structure = read_sdf(SDF).unwrap();
    assert_eq!(structure.sites().len(), 3);
    assert_eq!(structure.sites()[0].element(), "O");
    assert!((structure.sites()[0].charge().unwrap() + 1.0).abs() < 1e-9);
    assert!((structure.sites()[1].charge().unwrap() - 1.0).abs() < 1e-9);
    assert_eq!(structure.explicit_bonds().len(), 2);
    let data = structure.metadata().get("sdf_data").unwrap();
    assert_eq!(data.get("ID").unwrap(), "water");
}

#[test]
fn sdf_round_trip_keeps_charges_bonds_and_data() {
    let structure = read_sdf(SDF).unwrap();
    let text = write_sdf(&structure).unwrap();
    assert!(text.contains("V2000"));
    assert!(text.contains("M  CHG"));
    assert!(text.contains("$$$$"));
    let back = read_sdf(&text).unwrap();
    assert_eq!(back.sites().len(), 3);
    assert!((back.sites()[0].charge().unwrap() + 1.0).abs() < 1e-9);
    assert_eq!(back.explicit_bonds().len(), 2);
    assert_eq!(
        back.metadata()
            .get("sdf_data")
            .and_then(|value| value.get("ID"))
            .and_then(serde_json::Value::as_str),
        Some("water")
    );
}

#[test]
fn sdf_reads_v3000_and_multiple_records() {
    let structure = read_sdf(SDF_V3000).unwrap();
    assert_eq!(structure.sites().len(), 2);
    assert_eq!(structure.sites()[0].element(), "C");
    assert!((structure.sites()[1].charge().unwrap() + 1.0).abs() < 1e-9);
    assert_eq!(structure.explicit_bonds().len(), 1);
    assert!((structure.explicit_bonds()[0].order - 2.0).abs() < 1e-9);
    let two = format!("{SDF}{SDF_V3000}");
    assert_eq!(read_sdf_all(&two).unwrap().len(), 2);
}

#[test]
fn sdf_reports_malformed_records() {
    assert!(read_sdf("").is_err());
    assert!(read_sdf("title\nprogram\ncomment\nxx yy\n").is_err());
    let short = SDF.replace(
        "    1.0000    0.0000    0.0000 H   0  0  0",
        "    1.0000    0.0000",
    );
    assert!(read_sdf(&short).is_err());
    let no_counts = SDF_V3000.replace("M  V30 COUNTS 2 1 0 0 0", "");
    assert!(read_sdf(&no_counts).is_err());
}

#[test]
fn molecular_formats_are_registered_and_sniffed() {
    for (extension, name) in [
        (".mol2", "mol2"),
        (".sdf", "sdf"),
        (".mol", "sdf"),
        (".mdl", "sdf"),
    ] {
        let format = find(extension).unwrap();
        assert_eq!(format.name, name, "{extension}");
        assert!(format.read.is_some(), "{extension}");
        assert!(format.write.is_some(), "{extension}");
    }
    assert_eq!(sniff(MOL2).unwrap().name, "mol2");
    assert_eq!(sniff(SDF).unwrap().name, "sdf");
    assert_eq!(read_with(MOL2, "mol2").unwrap().sites().len(), 3);
    assert_eq!(read_with(SDF, ".mol").unwrap().sites().len(), 3);
}
