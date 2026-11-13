//! Integration tests for the Quantum ESPRESSO `pw.x` input and output readers
//! and the input writer (Task 13.2.1).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::{Lattice, Site, StructureBuilder};
use xtalforge_io::registry::{find, read_with, sniff};
use xtalforge_io::{read_pwscf, write_pwscf};

const INPUT: &str = "\
&CONTROL
  calculation = 'scf'
  prefix = 'si'
/
&SYSTEM
  ibrav = 0
  nat = 2
  ntyp = 1
  ecutwfc = 30.0
/
&ELECTRONS
  conv_thr = 1.0d-8
/
ATOMIC_SPECIES
  Si 28.0855 Si.UPF
CELL_PARAMETERS angstrom
  0.000000 2.715000 2.715000
  2.715000 0.000000 2.715000
  2.715000 2.715000 0.000000
ATOMIC_POSITIONS crystal
  Si 0.00 0.00 0.00
  Si 0.25 0.25 0.25
K_POINTS gamma
";

const OUTPUT: &str = "\
 Program PWSCF v.7.0 starts on 13Nov2026

     lattice parameter (alat)  =      10.2600  a.u.
     unit-cell volume          =     271.5000 (a.u.)^3
     number of atoms/cell      =          2
     number of atomic types    =          1

     crystal axes: (cart. coord. in units of alat)
       a(1) = (  -0.500000   0.000000   0.500000 )
       a(2) = (   0.000000   0.500000   0.500000 )
       a(3) = (  -0.500000   0.500000   0.000000 )

     site n.     atom                  positions (alat units)
         1           Si  tau(   1) = (    0.0000000   0.0000000   0.0000000  )
         2           Si  tau(   2) = (    0.2500000   0.2500000   0.2500000  )

!    total energy              =     -22.84600000 Ry
     convergence has been achieved in   4 iterations
";

fn silicon() -> xtalforge_core::Structure {
    StructureBuilder::new("silicon", Lattice::cubic(5.43).unwrap())
        .site(Site::new("Si1", "Si", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("Si2", "Si", [0.25, 0.25, 0.25]).unwrap())
        .build()
        .unwrap()
}

#[test]
fn pwscf_input_reads_cell_species_and_positions() {
    let structure = read_pwscf(INPUT).unwrap();
    assert_eq!(structure.formula(), "Si2");
    assert_eq!(structure.sites().len(), 2);
    assert_eq!(structure.metadata().get("pwscf_ibrav").unwrap(), "0");
    assert_eq!(structure.metadata().get("pwscf_prefix").unwrap(), "si");
    let [a, b, c, alpha, beta, gamma] = structure.lattice().parameters();
    assert!((a - 3.839_5).abs() < 1e-3);
    assert!((b - 3.839_5).abs() < 1e-3);
    assert!((c - 3.839_5).abs() < 1e-3);
    assert!((alpha - 60.0).abs() < 1e-3);
    assert!((beta - 60.0).abs() < 1e-3);
    assert!((gamma - 60.0).abs() < 1e-3);
    assert!(structure.is_p1());
}

#[test]
fn pwscf_output_reads_final_positions_and_energy() {
    let structure = read_pwscf(OUTPUT).unwrap();
    assert_eq!(structure.formula(), "Si2");
    let expected = 5.43 * std::f64::consts::SQRT_2 / 2.0;
    assert!((structure.lattice().lengths()[0] - expected).abs() < 1e-3);
    let energy = structure
        .metadata()
        .get("pwscf_energy_ev")
        .and_then(serde_json::Value::as_f64)
        .unwrap();
    assert!((energy + 310.84).abs() < 0.1);
    assert_eq!(
        structure.metadata().get("pwscf_output"),
        Some(&serde_json::Value::Bool(true))
    );
}

#[test]
fn pwscf_writer_round_trips_a_p1_cell() {
    let text = write_pwscf(&silicon()).unwrap();
    assert!(text.contains("ibrav = 0"));
    assert!(text.contains("CELL_PARAMETERS angstrom"));
    assert!(text.contains("ATOMIC_POSITIONS crystal"));
    assert!(text.contains("K_POINTS gamma"));
    let back = read_pwscf(&text).unwrap();
    assert_eq!(back.formula(), "Si2");
    assert!((back.lattice().volume() - 160.1).abs() < 1e-2);
}

#[test]
fn pwscf_input_supports_ibrav_and_bohr_units() {
    let text = "\
&SYSTEM
  ibrav = 1
  celldm(1) = 10.0
  nat = 1
  ntyp = 1
/
ATOMIC_SPECIES
  C 12.011 C.UPF
ATOMIC_POSITIONS bohr
  C 0.0 0.0 0.0
";
    let structure = read_pwscf(text).unwrap();
    assert_eq!(structure.formula(), "C");
    assert!((structure.lattice().lengths()[0] - 10.0 * 0.529_177_210_903).abs() < 1e-9);
}

#[test]
fn pwscf_rejects_unsupported_and_incomplete_input() {
    let unsupported = "\
&SYSTEM
  ibrav = 5
  celldm(1) = 10.0
  nat = 1
  ntyp = 1
/
ATOMIC_SPECIES
  C 12.011 C.UPF
ATOMIC_POSITIONS crystal
  C 0.0 0.0 0.0
";
    assert!(read_pwscf(unsupported).is_err());
    let no_positions = INPUT.replace("ATOMIC_POSITIONS crystal", "POSITIONS crystal");
    assert!(read_pwscf(&no_positions).is_err());
    let bad_atom = INPUT.replace("Si 0.25 0.25 0.25", "Si 0.25 nope 0.25");
    assert!(read_pwscf(&bad_atom).is_err());
    assert!(read_pwscf("").is_err());
}

#[test]
fn pwscf_is_registered_and_sniffed() {
    assert_eq!(find(".pwi").unwrap().name, "pwscf");
    assert_eq!(find("pwscf").unwrap().name, "pwscf");
    assert!(find("pwscf").unwrap().read.is_some());
    assert!(find("pwscf").unwrap().write.is_some());
    assert_eq!(sniff(INPUT).unwrap().name, "pwscf");
    assert_eq!(read_with(OUTPUT, "pwscf").unwrap().formula(), "Si2");
}
