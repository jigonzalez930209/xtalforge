//! Integration tests for the LAMMPS, SHELX and PDB formats (Tasks 3.3.4-3.3.6).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::{ExplicitBond, Lattice, Site, Structure, StructureBuilder};
use xtalforge_io::registry::{find, read_with};
use xtalforge_io::{read_lammps, read_shelx, write_lammps, write_pdb, write_shelx};

fn triclinic() -> Structure {
    StructureBuilder::new(
        "triclinic",
        Lattice::from_parameters(10.0, 11.0, 12.0, 90.0, 95.0, 100.0).unwrap(),
    )
    .site(
        Site::new("C1", "C", [0.1, 0.2, 0.3])
            .unwrap()
            .with_charge(-0.4)
            .unwrap()
            .with_disorder_group(2),
    )
    .site(
        Site::new("O1", "O", [0.4, 0.5, 0.6])
            .unwrap()
            .with_occupancy(0.5)
            .unwrap()
            .with_charge(0.2)
            .unwrap(),
    )
    .explicit_bond(ExplicitBond {
        i: 0,
        j: 1,
        image: [0, 0, 0],
        order: 1.0,
        remove: false,
    })
    .build()
    .unwrap()
}

#[test]
fn lammps_round_trip_keeps_cell_sites_charges_and_bonds() {
    let structure = triclinic();
    let text = write_lammps(&structure).unwrap();
    assert!(text.contains("2 atoms"));
    assert!(text.contains("Atoms  # full"));
    assert!(text.contains("Masses"));
    assert!(text.contains("Bonds"));
    assert!(text.contains("xy xz yz"));

    let back = read_lammps(&text).unwrap();
    assert_eq!(back.sites().len(), 2);
    assert_eq!(back.formula(), "C O");
    let [a, b, c, alpha, beta, gamma] = back.lattice().parameters();
    assert!((a - 10.0).abs() < 1e-8);
    assert!((b - 11.0).abs() < 1e-8);
    assert!((c - 12.0).abs() < 1e-8);
    assert!((alpha - 90.0).abs() < 1e-8);
    assert!((beta - 95.0).abs() < 1e-8);
    assert!((gamma - 100.0).abs() < 1e-8);
    assert!((back.lattice().volume() - structure.lattice().volume()).abs() < 1e-6);
    assert!((back.sites()[0].charge().unwrap() + 0.4).abs() < 1e-6);
    assert!((back.sites()[1].charge().unwrap() - 0.2).abs() < 1e-6);
    assert_eq!(back.explicit_bonds().len(), 1);
    assert_eq!(back.explicit_bonds()[0].i, 0);
    assert_eq!(back.explicit_bonds()[0].j, 1);
}

#[test]
fn lammps_mass_mapping_reports_ambiguity() {
    let text = "\
LAMMPS data file

1 atoms
1 atom types

0.0 10.0 xlo xhi
0.0 10.0 ylo yhi
0.0 10.0 zlo zhi

Masses

1 58.693

Atoms # full

1 1 1 0.0 1.0 1.0 1.0
";
    let structure = read_lammps(text).unwrap();
    assert_eq!(structure.sites()[0].element(), "Ni");
    let warnings = structure.metadata().get("warnings").unwrap();
    assert!(warnings.as_array().is_some_and(|items| !items.is_empty()));
}

#[test]
fn lammps_rejects_malformed_input() {
    assert!(read_lammps("").is_err());
    assert!(read_lammps("1 atoms\n0 1 xlo xhi\n0 1 ylo yhi\n0 1 zlo zhi\n").is_err());
    let short_atom = "\
1 atoms
1 atom types

0.0 10.0 xlo xhi
0.0 10.0 ylo yhi
0.0 10.0 zlo zhi

Masses

1 12.011

Atoms # full

1 1 1 0.0 1.0 1.0
";
    assert!(read_lammps(short_atom).is_err());

    let bad_bond = "\
2 atoms
1 atom types

0.0 10.0 xlo xhi
0.0 10.0 ylo yhi
0.0 10.0 zlo zhi

Masses

1 12.011

Atoms # full

1 1 1 0.0 1.0 1.0 1.0
2 1 1 0.0 2.0 1.0 1.0

Bonds

1 1 1 3
";
    assert!(read_lammps(bad_bond).is_err());
}

#[test]
fn shelx_round_trip_keeps_cell_occupancy_and_part() {
    let structure = triclinic();
    let text = write_shelx(&structure).unwrap();
    assert!(text.starts_with("TITL triclinic"));
    assert!(text.contains("CELL 0.71073"));
    assert!(text.contains("LATT 1"));
    assert!(text.contains("SFAC C O"));

    let back = read_shelx(&text).unwrap();
    assert_eq!(back.sites().len(), 2);
    let [a, b, c, alpha, beta, gamma] = back.lattice().parameters();
    assert!((a - 10.0).abs() < 1e-6);
    assert!((b - 11.0).abs() < 1e-6);
    assert!((c - 12.0).abs() < 1e-6);
    assert!((alpha - 90.0).abs() < 1e-6);
    assert!((beta - 95.0).abs() < 1e-6);
    assert!((gamma - 100.0).abs() < 1e-6);
    assert!((back.sites()[0].occupancy() - 1.0).abs() < 1e-6);
    assert!((back.sites()[1].occupancy() - 0.5).abs() < 1e-6);
    assert_eq!(back.sites()[0].disorder_group(), Some(2));
    assert_eq!(back.sites()[1].disorder_group(), None);
}

#[test]
fn shelx_resolves_free_variables_and_part_instructions() {
    let text = "\
TITL free
CELL 0.71073 10 10 10 90 90 90
LATT 1
SFAC C
FVAR 1.0 0.5
PART 3
C1 1 b 0.2 0.3 21.0 0.02
END
";
    let structure = read_shelx(text).unwrap();
    let site = &structure.sites()[0];
    assert!((site.frac()[0] - 0.5).abs() < 1e-9);
    assert!((site.frac()[1] - 0.2).abs() < 1e-9);
    assert!((site.occupancy() - 0.5).abs() < 1e-9);
    assert_eq!(site.disorder_group(), Some(3));
}

#[test]
fn shelx_reads_element_first_atom_rows() {
    let text = "\
CELL 0.71073 10 10 10 90 90 90
SFAC C O
C 0 0.1 0.3 10.5 0.02 2
END
";
    let structure = read_shelx(text).unwrap();
    let site = &structure.sites()[0];
    assert_eq!(site.element(), "C");
    assert!((site.frac()[0] - 0.0).abs() < 1e-9);
    assert!((site.frac()[1] - 0.1).abs() < 1e-9);
    assert!((site.occupancy() - 0.5).abs() < 1e-9);
    assert_eq!(site.disorder_group(), Some(2));
}

#[test]
fn shelx_rejects_malformed_input() {
    assert!(read_shelx("SFAC C\nC1 1 0.1 0.2 0.3 11 0.02\n").is_err());
    let no_sfac = "\
CELL 0.71073 10 10 10 90 90 90
SFAC C
C1 2 0.1 0.2 0.3 11 0.02
END
";
    assert!(read_shelx(no_sfac).is_err());
    let bad_occupancy = "\
CELL 0.71073 10 10 10 90 90 90
SFAC C
C1 1 0.1 0.2 0.3 5.0 0.02
END
";
    assert!(read_shelx(bad_occupancy).is_err());
    let missing_variable = "\
CELL 0.71073 10 10 10 90 90 90
SFAC C
FVAR 1.0
C1 1 z 0.2 0.3 11 0.02
END
";
    assert!(read_shelx(missing_variable).is_err());
}

#[test]
fn pdb_writes_hetatm_conect_and_end() {
    let structure = triclinic();
    let text = write_pdb(&structure).unwrap();
    let atoms: Vec<&str> = text
        .lines()
        .filter(|line| line.starts_with("HETATM"))
        .collect();
    assert_eq!(atoms.len(), 2);
    assert_eq!(atoms[0].len(), 78);
    assert_eq!(&atoms[0][76..78], " C");
    assert_eq!(&atoms[1][76..78], " O");
    let conect: Vec<&str> = text
        .lines()
        .filter(|line| line.starts_with("CONECT"))
        .collect();
    assert_eq!(conect.len(), 2);
    assert!(conect.iter().any(|line| line.starts_with("CONECT    1")));
    assert!(text.ends_with("END\n"));
}

#[test]
fn pdb_without_bonds_has_no_conect() {
    let structure = StructureBuilder::new("one", Lattice::cubic(5.0).unwrap())
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .build()
        .unwrap();
    let text = write_pdb(&structure).unwrap();
    assert!(!text.contains("CONECT"));
    assert_eq!(
        text.lines()
            .filter(|line| line.starts_with("HETATM"))
            .count(),
        1
    );
}

#[test]
fn registry_lists_the_new_formats_and_capabilities() {
    assert_eq!(find(".data").unwrap().name, "lammps-data");
    assert_eq!(find(".lmp").unwrap().name, "lammps-data");
    assert_eq!(find("lammps-data").unwrap().name, "lammps-data");
    assert_eq!(find(".res").unwrap().name, "shelx");
    assert_eq!(find(".ins").unwrap().name, "shelx");
    assert_eq!(find("shelx").unwrap().name, "shelx");
    assert_eq!(find(".pdb").unwrap().name, "pdb");
    assert!(find("pdb").unwrap().read.is_some());
    assert!(find("pdb").unwrap().write.is_some());
    assert!(find("pdb").unwrap().capabilities.cell);
    assert!(read_with("", "pdb").is_err());

    let lammps = find("lammps-data").unwrap().capabilities;
    assert!(lammps.cell && lammps.charges && lammps.bonds);
    assert!(!lammps.symmetry && !lammps.anisotropic_adp && !lammps.occupancy);
    let shelx = find("shelx").unwrap().capabilities;
    assert!(shelx.cell && shelx.occupancy && shelx.labels);
    assert!(!shelx.symmetry && !shelx.anisotropic_adp && !shelx.bonds);
    let pdb = find("pdb").unwrap().capabilities;
    assert!(pdb.cell && pdb.symmetry && pdb.bonds && pdb.labels);
    assert!(!pdb.anisotropic_adp);
}
