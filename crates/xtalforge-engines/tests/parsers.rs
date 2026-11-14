//! Integration tests for the engine output parsers (Task 18.2).
//!
//! Fixtures under `tests/fixtures/engines/` are hand-written representative
//! files; their provenance and license are recorded in the README there.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_engines::Error;
use xtalforge_engines::parse::{
    parse_cp2k, parse_gulp, parse_lammps, parse_lammps_dump, parse_lammps_dump_with_types,
    parse_lammps_log, parse_qe, parse_qe_bands, parse_qe_projwfc, parse_qe_xml, parse_vasp_doscar,
    parse_vasp_eigenval, parse_vasp_oszicar, parse_vasp_outcar, parse_vasp_procar,
    parse_vasp_vasprun,
};
use xtalforge_engines::result::{
    CalculationResult, bader_charges, bader_table, ddec_charges, ddec_table,
};

const VASPRUN: &str = include_str!("fixtures/engines/vasp/vasprun.xml");
const OUTCAR: &str = include_str!("fixtures/engines/vasp/OUTCAR");
const OSZICAR: &str = include_str!("fixtures/engines/vasp/OSZICAR");
const EIGENVAL: &str = include_str!("fixtures/engines/vasp/EIGENVAL");
const DOSCAR: &str = include_str!("fixtures/engines/vasp/DOSCAR");
const PROCAR: &str = include_str!("fixtures/engines/vasp/PROCAR");
const QE_PW: &str = include_str!("fixtures/engines/qe/pw.out");
const QE_XML: &str = include_str!("fixtures/engines/qe/data-file-schema.xml");
const QE_BANDS: &str = include_str!("fixtures/engines/qe/bands.dat");
const QE_PROJWFC: &str = include_str!("fixtures/engines/qe/projwfc.out");
const CP2K: &str = include_str!("fixtures/engines/cp2k/geo_opt.out");
const LAMMPS_LOG: &str = include_str!("fixtures/engines/lammps/minimize.log");
const LAMMPS_DUMP: &str = include_str!("fixtures/engines/lammps/final.dump");
const GULP: &str = include_str!("fixtures/engines/gulp/optimise.out");
const BADER: &str = include_str!("fixtures/engines/bader/ACF.dat");
const DDEC: &str = include_str!("fixtures/engines/bader/ddec6_charges.xyz");
const DUMP_TYPES: &str = "\
ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS pp pp pp
0.0 4.21
0.0 4.21
0.0 4.21
ITEM: ATOMS id type x y z
1 1 0.0 0.0 0.0
2 2 2.105 2.105 2.105
";

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 1e-6,
        "expected {expected}, got {actual}"
    );
}

#[test]
fn vasprun_reads_trajectory_bands_and_dos() {
    let result = parse_vasp_vasprun(VASPRUN).unwrap();
    assert_eq!(result.engine.as_deref(), Some("vasp"));
    assert_eq!(result.trajectory.len(), 2);
    assert_eq!(result.energies.len(), 2);
    close(result.energies[1], -10.7);
    let structure = result.structure.as_ref().unwrap();
    assert_eq!(structure.sites().len(), 2);
    close(structure.sites()[1].frac()[0], 0.49);
    let forces = result.forces.as_ref().unwrap();
    close(forces[1][2], -0.02);
    let stress = result.stress.unwrap();
    close(stress[0][0], 0.15);
    close(stress[1][1], 0.25);
    close(stress[2][2], 0.35);
    assert_eq!(result.magnetic_moments, vec![0.02, -0.02]);
    close(result.fermi.unwrap(), 5.1234);
    assert!(result.convergence.scf);
    assert!(result.convergence.ionic);
    assert!(result.convergence.ionic_reported);
    assert_eq!(result.steps.len(), 2);
    close(result.steps[1].max_force.unwrap(), 0.02);
    let bands = result.bands.as_ref().unwrap();
    assert_eq!(bands.kpoints, vec![[0.0, 0.0, 0.0]]);
    assert_eq!(bands.energies, vec![vec![-1.0, 2.0]]);
    assert_eq!(bands.occupations, vec![vec![1.0, 0.0]]);
    let dos = result.dos.as_ref().unwrap();
    assert_eq!(dos.total, vec![0.0, 1.0, 2.0]);
    assert_eq!(dos.integrated, vec![0.0, 0.5, 2.0]);
    assert_eq!(result.timings.len(), 2);
    assert_eq!(result.warnings.len(), 1);
    assert_eq!(result.warnings[0].code, "vasp_scf_convergence_inferred");
}

#[test]
fn outcar_reads_forces_moments_and_timings() {
    let result = parse_vasp_outcar(OUTCAR).unwrap();
    assert_eq!(result.trajectory.len(), 1);
    close(result.energies[0], -10.5);
    let structure = result.structure.as_ref().unwrap();
    assert_eq!(structure.sites().len(), 2);
    close(result.forces.as_ref().unwrap()[0][2], 0.01);
    let stress = result.stress.unwrap();
    close(stress[0][0], 0.1);
    close(stress[1][1], 0.2);
    assert_eq!(result.magnetic_moments, vec![0.01, -0.01]);
    close(result.fermi.unwrap(), 5.1234);
    close(result.scalars["total_charge"], 20.0);
    assert!(result.convergence.scf && result.convergence.ionic);
    assert_eq!(result.timings[0].stage, "cpu");
    close(result.timings[0].seconds, 12.345);
    assert!(result.warnings.is_empty());
}

#[test]
fn oszicar_reads_electronic_and_ionic_steps() {
    let result = parse_vasp_oszicar(OSZICAR).unwrap();
    assert_eq!(result.energies, vec![-10.5, -10.7]);
    assert_eq!(result.steps.len(), 3);
    assert!(result.convergence.ionic);
    assert!(result.convergence.scf);
    close(result.scalars["electronic_steps"], 3.0);
    close(result.scalars["ionic_steps"], 2.0);
    assert!(result.structure.is_none());
}

#[test]
fn eigenval_reads_bands() {
    let result = parse_vasp_eigenval(EIGENVAL).unwrap();
    let bands = result.bands.as_ref().unwrap();
    assert_eq!(bands.kpoints.len(), 2);
    assert_eq!(bands.weights, vec![1.0, 1.0]);
    assert_eq!(bands.energies[0], vec![-10.0, 2.0, 5.0, 8.0]);
    assert_eq!(bands.energies[1], vec![-9.5, 3.0, 6.0, 9.0]);
    assert_eq!(bands.occupations[0], vec![1.0, 0.0, 0.0, 0.0]);
    close(result.scalars["nbands"], 4.0);
    assert_eq!(result.warnings[0].code, "eigenval_has_no_convergence");
    assert!(result.structure.is_none());
}

#[test]
fn doscar_reads_total_and_integrated_dos() {
    let result = parse_vasp_doscar(DOSCAR).unwrap();
    let dos = result.dos.as_ref().unwrap();
    assert_eq!(dos.energies, vec![-5.0, 0.0, 5.0]);
    assert_eq!(dos.total, vec![0.5, 1.5, 0.2]);
    assert_eq!(dos.integrated, vec![0.1, 0.6, 1.0]);
    close(result.fermi.unwrap(), 0.0);
    close(result.scalars["natoms"], 2.0);
    assert_eq!(result.warnings.len(), 1);
    assert_eq!(result.warnings[0].code, "doscar_has_no_convergence");
}

#[test]
fn procar_reads_bands_and_orbital_weights() {
    let result = parse_vasp_procar(PROCAR).unwrap();
    let bands = result.bands.as_ref().unwrap();
    assert_eq!(bands.kpoints, vec![[0.0, 0.0, 0.0]]);
    assert_eq!(bands.energies, vec![vec![-10.0, 2.0]]);
    let projections = result.projections.as_ref().unwrap();
    assert_eq!(projections.bands.len(), 2);
    let first = &projections.bands[0];
    assert_eq!(first.band, 0);
    close(first.energy, -10.0);
    assert_eq!(first.weights.len(), 6);
    assert_eq!(first.weights[0].atom, 0);
    assert_eq!(first.weights[0].orbital, "s");
    close(first.weights[0].weight, 0.1);
    assert_eq!(projections.bands[1].band, 1);
    close(projections.bands[1].weights[1].weight, 0.15);
    assert_eq!(result.warnings[0].code, "procar_has_no_convergence");
}

#[test]
fn qe_text_reads_structure_energy_forces_and_stress() {
    let result = parse_qe(QE_PW).unwrap();
    assert_eq!(result.engine.as_deref(), Some("quantum_espresso"));
    let structure = result.structure.as_ref().unwrap();
    assert_eq!(structure.sites().len(), 2);
    assert_eq!(result.energies.len(), 1);
    close(result.energies[0], -22.846 * 13.605_693_122_994);
    assert_eq!(result.steps.len(), 2);
    let forces = result.forces.as_ref().unwrap();
    close(forces[0][0], 0.001 * 25.711_043_095_4);
    close(forces[1][2], -0.003 * 25.711_043_095_4);
    let stress = result.stress.unwrap();
    close(stress[0][0], 0.00001 * 14_710.0);
    close(stress[2][2], 0.00003 * 14_710.0);
    close(result.fermi.unwrap(), 5.1234);
    assert!(result.convergence.scf && result.convergence.scf_reported);
    assert_eq!(result.timings[0].stage, "PWSCF");
    close(result.timings[0].seconds, 1.23);
    assert_eq!(result.trajectory.len(), 1);
    assert!(result.warnings.is_empty());
}

#[test]
fn qe_xml_reads_structure_bands_and_stress() {
    let result = parse_qe_xml(QE_XML).unwrap();
    let structure = result.structure.as_ref().unwrap();
    assert_eq!(structure.sites().len(), 2);
    assert!(result.convergence.scf && result.convergence.scf_reported);
    close(result.energies[0], -11.54321 * 27.211_386_245_988);
    close(result.fermi.unwrap(), 5.1234);
    let bands = result.bands.as_ref().unwrap();
    assert_eq!(bands.kpoints.len(), 2);
    assert_eq!(bands.weights, vec![1.0, 1.0]);
    close(bands.energies[0][0], -0.4 * 27.211_386_245_988);
    close(bands.energies[1][1], 0.3 * 27.211_386_245_988);
    assert_eq!(bands.occupations[0], vec![1.0, 0.0]);
    let forces = result.forces.as_ref().unwrap();
    close(forces[0][0], 0.001 * 51.422_067_476_3);
    let stress = result.stress.unwrap();
    close(stress[0][0], 0.1);
    close(stress[2][2], 0.3);
    close(result.timings[0].seconds, 1.25);
    assert_eq!(result.trajectory.len(), 1);
    assert!(result.warnings.is_empty());
}

#[test]
fn qe_bands_reads_kpoints_and_bands() {
    let result = parse_qe_bands(QE_BANDS).unwrap();
    let bands = result.bands.as_ref().unwrap();
    assert_eq!(bands.kpoints.len(), 2);
    assert_eq!(bands.energies, vec![vec![-10.0, 2.0], vec![-9.5, 3.0]]);
    close(result.scalars["nbnd"], 2.0);
    assert_eq!(result.warnings[0].code, "qe_bands_no_fermi");
}

#[test]
fn qe_projwfc_reads_bands_and_lowdin_charges() {
    let result = parse_qe_projwfc(QE_PROJWFC).unwrap();
    let bands = result.bands.as_ref().unwrap();
    assert_eq!(bands.kpoints.len(), 2);
    assert_eq!(bands.energies, vec![vec![-10.0, 2.0], vec![-9.5, 3.0]]);
    assert_eq!(bands.occupations[0], vec![1.0, 0.0]);
    assert_eq!(result.charges.len(), 2);
    close(result.charges[0], 10.1234);
    close(result.charges[1], 9.8766);
    assert_eq!(result.warnings[0].code, "qe_projwfc_no_fermi");
}

#[test]
fn cp2k_reads_energy_forces_stress_and_convergence() {
    let result = parse_cp2k(CP2K).unwrap();
    let structure = result.structure.as_ref().unwrap();
    assert_eq!(structure.sites().len(), 2);
    close(structure.lattice().parameters()[0], 4.21);
    close(result.energies[0], -17.123456789 * 27.211_386_245_988);
    let forces = result.forces.as_ref().unwrap();
    close(forces[0][2], 0.001 * 51.422_067_476_3);
    close(forces[1][2], -0.001 * 51.422_067_476_3);
    let stress = result.stress.unwrap();
    close(stress[0][0], 0.1);
    close(stress[2][2], 0.3);
    assert!(result.convergence.scf && result.convergence.ionic);
    close(result.timings[0].seconds, 1.5);
    assert_eq!(result.trajectory.len(), 1);
    assert!(result.warnings.is_empty());
}

#[test]
fn lammps_log_reads_thermo_and_timing() {
    let result = parse_lammps_log(LAMMPS_LOG).unwrap();
    assert_eq!(result.energies, vec![-99.5, -99.9]);
    close(result.scalars["pressure_gpa"], 0.09);
    close(result.scalars["temperature"], 310.0);
    assert!(result.convergence.ionic && result.convergence.ionic_reported);
    assert_eq!(result.timings[0].stage, "loop");
    close(result.timings[0].seconds, 0.125);
    assert_eq!(result.warnings[0].code, "lammps_log_has_no_structure");
    assert!(result.structure.is_none());
}

#[test]
fn lammps_dump_reads_final_structure() {
    let result = parse_lammps_dump(LAMMPS_DUMP).unwrap();
    let structure = result.structure.as_ref().unwrap();
    assert_eq!(structure.sites().len(), 2);
    close(structure.lattice().parameters()[0], 4.21);
    close(structure.sites()[1].frac()[0], 0.5);
    assert!(result.energies.is_empty());
    assert_eq!(result.warnings[0].code, "lammps_dump_has_no_thermo");
    let detected = parse_lammps(LAMMPS_DUMP).unwrap();
    assert_eq!(detected.trajectory.len(), 1);
    let typed = parse_lammps_dump_with_types(DUMP_TYPES, &["Mg", "O"]).unwrap();
    assert_eq!(typed.structure.unwrap().sites()[1].element(), "O");
    let error = parse_lammps_dump_with_types(DUMP_TYPES, &["Mg"]).unwrap_err();
    assert_eq!(error.code(), "engine_parse");
}

#[test]
fn gulp_reads_structure_energy_forces_and_stress() {
    let result = parse_gulp(GULP).unwrap();
    let structure = result.structure.as_ref().unwrap();
    assert_eq!(structure.sites().len(), 2);
    close(structure.lattice().parameters()[0], 4.21);
    close(result.energies[0], -123.456_789);
    let forces = result.forces.as_ref().unwrap();
    close(forces[0][2], 0.01);
    close(forces[1][2], -0.01);
    let stress = result.stress.unwrap();
    close(stress[0][0], 0.1);
    close(stress[1][1], 0.2);
    close(stress[2][2], 0.3);
    assert!(result.convergence.scf && result.convergence.ionic);
    assert_eq!(result.timings[0].stage, "cpu");
    close(result.timings[0].seconds, 0.25);
    assert_eq!(result.trajectory.len(), 1);
    assert!(result.warnings.is_empty());
}

#[test]
fn bader_and_ddec_charges_map_onto_site_order() {
    let charges = bader_charges(BADER).unwrap();
    assert_eq!(charges.len(), 2);
    close(charges[0], 10.1234);
    close(charges[1], 9.8766);
    let table = bader_table(BADER).unwrap();
    close(table[0].min_distance.unwrap(), 0.5234);
    close(table[1].volume.unwrap(), 13.1);
    let ddec = ddec_charges(DDEC).unwrap();
    assert_eq!(ddec.len(), 2);
    close(ddec[0], 10.123456);
    close(ddec[1], 9.876544);
    let entries = ddec_table(DDEC).unwrap();
    assert_eq!(entries[0].element.as_deref(), Some("Mg"));
    assert_eq!(entries[1].site, 2);
    close(entries[0].spin.unwrap(), 0.0);
}

#[test]
fn result_helpers_close_the_trajectory() {
    let parsed = parse_vasp_outcar(OUTCAR).unwrap();
    let mut result = CalculationResult::new("test");
    assert_eq!(result.engine.as_deref(), Some("test"));
    assert_eq!(result.final_energy(), None);
    result.energies.push(-1.0);
    assert_eq!(result.final_energy(), Some(-1.0));
    result.push_timing("ignored", -1.0);
    assert!(result.timings.is_empty());
    result.push_timing("total", 1.5);
    assert_eq!(result.timings.len(), 1);
    result.push_warning("code", "message");
    assert_eq!(result.warnings[0].code, "code");
    result.trajectory = parsed.trajectory;
    result.close_trajectory();
    assert_eq!(
        result.structure.as_ref().unwrap().sites().len(),
        parsed.structure.unwrap().sites().len()
    );
    assert!(result.forces.is_some());
}

#[test]
fn truncated_outputs_return_typed_parse_errors() {
    let cases: &[(&str, &str)] = &[
        ("vasprun", "<modeling><calculation>"),
        (
            "outcar",
            " vasp.6.4.2\n TITEL = PAW_PBE Mg\n ions per type = 1\n",
        ),
        ("oszicar", "DAV:   1    -0.105E+02\n"),
        (
            "eigenval",
            "   2   4   1   1\n  0.1\n  CAR\n  Mg O\n   1   1\n",
        ),
        (
            "doscar",
            "   2   1   1   1\n  0.1 0.0 0.0 0.0 1.0\n  CAR\n  Mg O\n   1   1\n  5.0 0.0 3.0 0.0 1.0\n",
        ),
        (
            "procar",
            " PRO/CAR\n# of k-points: 1 # of bands: 1 # of ions: 1\n",
        ),
        (
            "qe",
            "     Program PWSCF v.7.2 starts\n     lattice parameter (alat) = 10.26 a.u.\n",
        ),
        ("qe_xml", "<qes:espresso><output><total_energy>"),
        (
            "qe_bands",
            " &plot nbnd=   2, nks=   2 /\n 0.0 0.0 0.0\n -10.0 2.0\n",
        ),
        ("qe_projwfc", "     Program PROJWFC v.7.2 starts\n"),
        (
            "cp2k",
            " CP2K version 2024.1\n *** SCF run converged in 10 steps ***\n",
        ),
        ("lammps_log", " LAMMPS (2 Aug 2023)\n Step Temp E_pair\n"),
        (
            "lammps_dump",
            "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n2\n",
        ),
        ("gulp", " GULP 6.0.0\n"),
    ];
    for (label, text) in cases {
        let error = parse_case(label, text).unwrap_err();
        assert_eq!(error.code(), "engine_parse", "case {label}: {error}");
    }
    assert_eq!(
        bader_charges("no rows here").unwrap_err().code(),
        "engine_parse"
    );
    assert_eq!(
        ddec_charges("no charge table").unwrap_err().code(),
        "engine_parse"
    );
}

fn parse_case(label: &str, text: &str) -> Result<CalculationResult, Error> {
    match label {
        "vasprun" => parse_vasp_vasprun(text),
        "outcar" => parse_vasp_outcar(text),
        "oszicar" => parse_vasp_oszicar(text),
        "eigenval" => parse_vasp_eigenval(text),
        "doscar" => parse_vasp_doscar(text),
        "procar" => parse_vasp_procar(text),
        "qe" => parse_qe(text),
        "qe_xml" => parse_qe_xml(text),
        "qe_bands" => parse_qe_bands(text),
        "qe_projwfc" => parse_qe_projwfc(text),
        "cp2k" => parse_cp2k(text),
        "lammps_log" => parse_lammps_log(text),
        "lammps_dump" => parse_lammps_dump(text),
        "gulp" => parse_gulp(text),
        other => panic!("unknown case {other}"),
    }
}
