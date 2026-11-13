//! Integration tests for the VASP `OUTCAR` and `vasprun.xml` readers and the
//! trajectory detection (Task 13.2.4).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_io::registry::{find, sniff};
use xtalforge_io::{
    read_outcar, read_outcar_structure, read_trajectory, read_vasprun, read_vasprun_structure,
};

const OUTCAR: &str = "\
 vasp.6.3.0 13Nov2026

 TITEL  = PAW_PBE Si 05Jan2001
 ions per type =               2
 direct lattice vectors                 reciprocal lattice vectors
     0.000000000  2.715000000  2.715000000     0.368324    0.368324    0.368324
     2.715000000  0.000000000  2.715000000     0.368324    0.368324    0.368324
     2.715000000  2.715000000  0.000000000     0.368324    0.368324    0.368324

 POSITION                                       TOTAL-FORCE (eV/Angst)
 -----------------------------------------------------------------------------------
      0.00000      0.00000      0.00000         0.100000      0.200000      0.300000
      1.35750      1.35750      1.35750         0.000000      0.000000      0.000000
 -----------------------------------------------------------------------------------

  FORCE on cell =-STRESS in cart. coord.  units (eV):
  Direction    XX          YY          ZZ          XY          YZ          ZX
  --------------------------------------------------------------------------------------
  Total     0.10 0.20 0.30 0.40 0.50 0.60
  in kB     1.00 2.00 3.00 4.00 5.00 6.00

 free  energy   TOTEN  =        -10.80000000 eV

 POSITION                                       TOTAL-FORCE (eV/Angst)
 -----------------------------------------------------------------------------------
      0.00000      0.00000      0.00000         0.050000      0.000000      0.000000
      1.35750      1.35750      1.35750         0.000000      0.000000      0.000000
 -----------------------------------------------------------------------------------

  in kB     2.00 3.00 4.00 5.00 6.00 7.00

 free  energy   TOTEN  =        -10.90000000 eV
";

const VASPRUN: &str = "\
<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>
<modeling>
 <generator>
  <i name=\"program\">vasp</i>
  <i name=\"version\">6.3.0</i>
 </generator>
 <atominfo>
  <array name=\"atoms\">
   <set>
    <rc><c>Si</c><c>1</c></rc>
    <rc><c>Si</c><c>2</c></rc>
   </set>
  </array>
 </atominfo>
 <structure name=\"initialpos\">
  <crystal>
   <varray name=\"basis\">
    <v>0.0 2.715 2.715</v>
    <v>2.715 0.0 2.715</v>
    <v>2.715 2.715 0.0</v>
   </varray>
  </crystal>
  <varray name=\"positions\">
   <v>0.0 0.0 0.0</v>
   <v>0.25 0.25 0.25</v>
  </varray>
 </structure>
 <calculation>
  <structure>
   <crystal>
    <varray name=\"basis\">
     <v>0.0 2.715 2.715</v>
     <v>2.715 0.0 2.715</v>
     <v>2.715 2.715 0.0</v>
    </varray>
   </crystal>
   <varray name=\"positions\">
    <v>0.0 0.0 0.0</v>
    <v>0.25 0.25 0.25</v>
   </varray>
  </structure>
  <energy>
   <i name=\"e_fr_energy\">-10.8</i>
  </energy>
  <varray name=\"forces\">
   <v>0.1 0.2 0.3</v>
   <v>0.0 0.0 0.0</v>
  </varray>
  <varray name=\"stress\">
   <v>1.0 2.0 3.0</v>
   <v>2.0 4.0 5.0</v>
   <v>3.0 5.0 6.0</v>
  </varray>
 </calculation>
 <calculation>
  <structure>
   <crystal>
    <varray name=\"basis\">
     <v>0.0 2.715 2.715</v>
     <v>2.715 0.0 2.715</v>
     <v>2.715 2.715 0.0</v>
    </varray>
   </crystal>
   <varray name=\"positions\">
    <v>0.0 0.0 0.0</v>
    <v>0.26 0.25 0.24</v>
   </varray>
  </structure>
  <energy>
   <i name=\"e_fr_energy\">-10.9</i>
  </energy>
  <varray name=\"forces\">
   <v>0.0 0.0 0.0</v>
   <v>0.0 0.0 0.0</v>
  </varray>
 </calculation>
 <structure name=\"finalpos\">
  <crystal>
   <varray name=\"basis\">
    <v>0.0 2.715 2.715</v>
    <v>2.715 0.0 2.715</v>
    <v>2.715 2.715 0.0</v>
   </varray>
  </crystal>
  <varray name=\"positions\">
   <v>0.0 0.0 0.0</v>
   <v>0.26 0.25 0.24</v>
  </varray>
 </structure>
</modeling>
";

#[test]
fn outcar_reads_steps_energies_forces_and_stress() {
    let trajectory = read_outcar(OUTCAR).unwrap();
    assert_eq!(trajectory.elements, vec!["Si", "Si"]);
    assert_eq!(trajectory.frames.len(), 2);
    assert!(trajectory.lattice.is_some());
    let first = &trajectory.frames[0];
    assert!((first.energy.unwrap() + 10.8).abs() < 1e-9);
    let forces = first.forces.as_ref().unwrap();
    assert!((forces[0][0] - 0.1).abs() < 1e-9);
    assert!((forces[0][2] - 0.3).abs() < 1e-9);
    let stress = first.stress.unwrap();
    assert!((stress[0][0] - 0.1).abs() < 1e-9);
    assert!((stress[1][1] - 0.2).abs() < 1e-9);
    assert!((stress[2][2] - 0.3).abs() < 1e-9);
    assert!((first.frac[1][0] - 0.25).abs() < 1e-9);
    assert!((trajectory.frames[1].energy.unwrap() + 10.9).abs() < 1e-9);
    assert!((trajectory.frames[1].stress.unwrap()[0][0] - 0.2).abs() < 1e-9);
}

#[test]
fn outcar_structure_is_the_last_step() {
    let structure = read_outcar_structure(OUTCAR).unwrap();
    assert_eq!(structure.formula(), "Si2");
    assert_eq!(structure.sites().len(), 2);
    assert!((structure.sites()[1].frac()[0] - 0.25).abs() < 1e-9);
    let expected = 2.0 * 2.715_f64.powi(3);
    assert!((structure.lattice().volume() - expected).abs() < 1e-9);
}

#[test]
fn vasprun_reads_calculations_forces_stress_and_final_structure() {
    let trajectory = read_vasprun(VASPRUN).unwrap();
    assert_eq!(trajectory.elements, vec!["Si", "Si"]);
    assert_eq!(trajectory.frames.len(), 2);
    assert!(trajectory.lattice.is_some());
    assert!(trajectory.frames[0].lattice.is_none());
    assert!((trajectory.frames[0].energy.unwrap() + 10.8).abs() < 1e-9);
    let forces = trajectory.frames[0].forces.as_ref().unwrap();
    assert!((forces[0][1] - 0.2).abs() < 1e-9);
    let stress = trajectory.frames[0].stress.unwrap();
    assert!((stress[0][1] - 0.2).abs() < 1e-9);
    assert!((stress[2][2] - 0.6).abs() < 1e-9);
    assert_eq!(
        trajectory.frames[0].info.get("program").map(String::as_str),
        Some("vasp")
    );
    assert!((trajectory.frames[1].frac[1][0] - 0.26).abs() < 1e-9);
    let structure = read_vasprun_structure(VASPRUN).unwrap();
    assert_eq!(structure.formula(), "Si2");
    assert!((structure.sites()[1].frac()[0] - 0.26).abs() < 1e-9);
}

#[test]
fn read_trajectory_detects_both_vasp_outputs() {
    assert_eq!(read_trajectory(OUTCAR).unwrap().frames.len(), 2);
    assert_eq!(read_trajectory(VASPRUN).unwrap().frames.len(), 2);
}

#[test]
fn vasp_outputs_report_malformed_files() {
    assert!(read_outcar("").is_err());
    assert!(read_outcar("TITEL = PAW_PBE Si\nions per type = 1\n").is_err());
    let truncated = VASPRUN.replace("</modeling>", "");
    assert!(read_vasprun(&truncated).is_err());
    let no_atoms = VASPRUN.replace("<c>Si</c><c>1</c>", "");
    assert!(read_vasprun_structure(&no_atoms).is_err());
}

#[test]
fn vasp_formats_are_registered_read_only() {
    let outcar = find(".outcar").unwrap();
    assert_eq!(outcar.name, "outcar");
    assert!(outcar.read.is_some());
    assert!(outcar.write.is_none());
    assert!(outcar.capabilities.frames);
    let vasprun = find(".vasprun.xml").unwrap();
    assert_eq!(vasprun.name, "vasprun");
    assert!(vasprun.read.is_some());
    assert!(vasprun.write.is_none());
    assert_eq!(sniff(OUTCAR).unwrap().name, "outcar");
    assert_eq!(sniff(VASPRUN).unwrap().name, "vasprun");
}
