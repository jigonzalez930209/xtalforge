//! Integration tests for the UFF/UFF4MOF calculator (Tasks 4.1 and 4.2).
//!
//! Expected energies come from the Crystal-Builder reference expressions:
//! ethylene is checked term by term against `xtal/ff/uff/terms.py`, diamond
//! against the same expressions summed over the periodic images, and the
//! simple molecules are placed at their UFF equilibrium geometry where the
//! answer is zero by construction.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::{Lattice, Site, Structure, StructureBuilder};
use xtalforge_ff::uff::params;
use xtalforge_ff::{
    Calculator, CalculatorOptions, CalculatorRegistry, ParameterSet, RunRecord, UffCalculator,
    UffModel,
};

fn molecule(name: &str, atoms: &[(&str, [f64; 3])]) -> Structure {
    let points: Vec<(String, [f64; 3])> = atoms
        .iter()
        .map(|(element, point)| ((*element).to_owned(), *point))
        .collect();
    Structure::molecule(name, &points, 8.0).expect("molecule")
}

fn tetrahedral_methane() -> Structure {
    let length = xtalforge_ff::uff::natural_bond_length("C_3", "H_", 1.0).unwrap();
    let a = length / 3.0_f64.sqrt();
    molecule(
        "methane",
        &[
            ("C", [0.0, 0.0, 0.0]),
            ("H", [a, a, a]),
            ("H", [a, -a, -a]),
            ("H", [-a, a, -a]),
            ("H", [-a, -a, a]),
        ],
    )
}

fn water() -> Structure {
    let length = xtalforge_ff::uff::natural_bond_length("O_3", "H_", 1.0).unwrap();
    let angle = params::get("O_3").unwrap().theta0.to_radians();
    molecule(
        "water",
        &[
            ("O", [0.0, 0.0, 0.0]),
            ("H", [length, 0.0, 0.0]),
            ("H", [length * angle.cos(), length * angle.sin(), 0.0]),
        ],
    )
}

fn diamond() -> Structure {
    let lattice = Lattice::cubic(3.566_986).unwrap();
    let mut builder = StructureBuilder::new("diamond", lattice);
    for (index, point) in [
        [0.0, 0.0, 0.0],
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.5],
        [0.5, 0.5, 0.0],
    ]
    .iter()
    .enumerate()
    {
        builder = builder.site(Site::new(format!("C{index}"), "C", *point).unwrap());
    }
    for (index, point) in [
        [0.25, 0.25, 0.25],
        [0.25, 0.75, 0.75],
        [0.75, 0.25, 0.75],
        [0.75, 0.75, 0.25],
    ]
    .iter()
    .enumerate()
    {
        builder = builder.site(Site::new(format!("D{index}"), "C", *point).unwrap());
    }
    builder.build().unwrap()
}

fn mof_fragment() -> Structure {
    let dirs = [
        [1.0, 1.0, 1.0],
        [1.0, -1.0, -1.0],
        [-1.0, 1.0, -1.0],
        [-1.0, -1.0, 1.0],
    ];
    let mut atoms: Vec<(&str, [f64; 3])> = vec![("Zn", [0.0, 0.0, 0.0])];
    for direction in dirs {
        let norm = 3.0_f64.sqrt();
        let nitrogen = [
            direction[0] / norm * 2.0,
            direction[1] / norm * 2.0,
            direction[2] / norm * 2.0,
        ];
        atoms.push(("N", nitrogen));
        let mut perpendicular = [direction[1], -direction[0], 0.0];
        let length = (perpendicular[0] * perpendicular[0]
            + perpendicular[1] * perpendicular[1]
            + perpendicular[2] * perpendicular[2])
            .sqrt();
        perpendicular = [
            perpendicular[0] / length,
            perpendicular[1] / length,
            perpendicular[2] / length,
        ];
        atoms.push((
            "C",
            [
                nitrogen[0] + perpendicular[0] * 1.35,
                nitrogen[1] + perpendicular[1] * 1.35,
                nitrogen[2] + perpendicular[2] * 1.35,
            ],
        ));
    }
    molecule("mof-fragment", &atoms)
}

fn ethene() -> Structure {
    molecule(
        "ethene",
        &[
            ("C", [0.05, 0.02, 0.03]),
            ("C", [1.38, -0.04, 0.01]),
            ("H", [-0.51, 0.98, -0.02]),
            ("H", [-0.56, -0.94, 0.06]),
            ("H", [1.95, 0.93, 0.04]),
            ("H", [1.91, -0.99, -0.03]),
        ],
    )
}

fn options(cutoff: f64) -> CalculatorOptions {
    CalculatorOptions {
        vdw_cutoff: cutoff,
        ..CalculatorOptions::default()
    }
}

#[test]
fn parameter_table_has_the_reference_rows() {
    assert_eq!(params::all().unwrap().len(), 218);
    assert_eq!(params::schema_version().unwrap(), 1);
    let carbon = params::get("C_3").unwrap();
    assert!((carbon.r1 - 0.757).abs() < 1e-12);
    assert!((carbon.theta0 - 109.47).abs() < 1e-12);
    assert!((carbon.x1 - 3.851).abs() < 1e-12);
    assert_eq!(params::get("Zn3f2").unwrap().element, "Zn");
    assert!(params::get("Zn3f2").unwrap().fitted);
    assert_eq!(params::get("Lw6+3").unwrap().element, "Lr");
    assert_eq!(params::types_for("Lr").unwrap().len(), 1);
    assert!(!params::types_for("H").unwrap().is_empty());
    assert!(!params::data_source().unwrap().is_empty());
}

#[test]
fn methane_is_at_its_uff_bond_length() {
    let model = UffModel::build(&tetrahedral_methane(), options(8.0)).unwrap();
    assert_eq!(model.types(), vec!["C_3", "H_", "H_", "H_", "H_"]);
    assert_eq!(model.n_bonds(), 4);
    assert_eq!(model.n_angles(), 6);
    assert!((model.energy().unwrap()).abs() < 1e-6);
    let length = xtalforge_ff::uff::natural_bond_length("C_3", "H_", 1.0).unwrap();
    assert!((length - 1.109_400_794_877_744).abs() < 1e-12);
}

#[test]
fn water_angle_is_the_uff_equilibrium() {
    let model = UffModel::build(&water(), options(8.0)).unwrap();
    assert_eq!(model.types(), vec!["O_3", "H_", "H_"]);
    assert_eq!(model.n_angles(), 1);
    assert_eq!(model.n_vdw_pairs(), 0);
    assert!(model.energy().unwrap().abs() < 1e-9);
}

#[test]
fn noble_gas_dimer_is_one_lennard_jones_well() {
    let distance = params::get("Ar4+4").unwrap().x1;
    let dimer = molecule(
        "argon-dimer",
        &[("Ar", [0.0, 0.0, 0.0]), ("Ar", [distance, 0.0, 0.0])],
    );
    let model = UffModel::build(&dimer, options(12.0)).unwrap();
    assert_eq!(model.types(), vec!["Ar4+4", "Ar4+4"]);
    assert_eq!(model.n_bonds(), 0);
    assert_eq!(model.n_vdw_pairs(), 1);
    // At r = x the unshifted well is -D; the cutoff shift removes the
    // attractive tail at 12 A, giving -0.0080043731 eV.
    assert!((model.energy().unwrap() - -0.008_004_373_1).abs() < 1e-7);
    let (_, forces) = model.energy_forces().unwrap();
    assert!(forces[0][0].abs() < 1e-6);
    assert!(forces[1][0].abs() < 1e-6);
}

#[test]
fn diamond_is_periodic_and_symmetric() {
    let model = UffModel::build(&diamond(), options(12.0)).unwrap();
    assert_eq!(model.n_atoms(), 8);
    assert_eq!(model.types(), vec!["C_3"; 8]);
    assert_eq!(model.n_bonds(), 16);
    assert_eq!(model.n_angles(), 48);
    assert_eq!(model.n_torsions(), 144);
    assert_eq!(model.n_inversions(), 0);
    assert!((model.energy().unwrap() - 1.720_045_716_9).abs() < 1e-6);
    let (_, forces) = model.energy_forces().unwrap();
    let maximum = forces
        .iter()
        .map(|force| force[0].abs().max(force[1].abs()).max(force[2].abs()))
        .fold(0.0_f64, f64::max);
    assert!(maximum < 1e-10, "diamond force {maximum:e}");
}

#[test]
fn mof_fragment_gets_uff4mof_parameters() {
    let fragment = mof_fragment();
    let mof = UffModel::build(&fragment, options(8.0)).unwrap();
    assert_eq!(mof.n_atoms(), 9);
    assert_eq!(mof.types().first().copied(), Some("Zn3f2"));
    let energy = mof.energy().unwrap();
    assert!(energy.is_finite() && energy > 0.0);

    let plain = UffModel::build(
        &fragment,
        CalculatorOptions {
            parameter_set: ParameterSet::Uff,
            ..options(8.0)
        },
    )
    .unwrap();
    assert_eq!(plain.types().first().copied(), Some("Zn3+2"));
    assert!((plain.energy().unwrap() - energy).abs() > 1e-3);
}

#[test]
fn dummy_atoms_never_contribute() {
    let clean = tetrahedral_methane();
    let mut points = vec![("C", [0.0, 0.0, 0.0])];
    let length = xtalforge_ff::uff::natural_bond_length("C_3", "H_", 1.0).unwrap();
    let a = length / 3.0_f64.sqrt();
    for point in [[a, a, a], [a, -a, -a], [-a, a, -a], [-a, -a, a]] {
        points.push(("H", point));
    }
    points.push(("X", [0.0, 0.0, 8.0]));
    let with_dummy = molecule("methane-with-marker", &points);

    let clean_model = UffModel::build(&clean, options(8.0)).unwrap();
    let dummy_model = UffModel::build(&with_dummy, options(8.0)).unwrap();
    assert_eq!(dummy_model.n_atoms(), 5);
    assert_eq!(dummy_model.n_p1_atoms(), 6);
    assert_eq!(dummy_model.types().len(), 5);
    assert!((dummy_model.energy().unwrap() - clean_model.energy().unwrap()).abs() < 1e-12);
    let (_, forces) = dummy_model.energy_forces().unwrap();
    assert_eq!(forces.len(), 6);
    assert_eq!(forces[5], [0.0, 0.0, 0.0]);
}

#[test]
fn analytic_forces_match_finite_differences() {
    let model = UffModel::build(&ethene(), options(8.0)).unwrap();
    // This geometry is the one checked against terms.py: all five terms are
    // present and the total is 0.1113230794 eV.
    assert!((model.energy().unwrap() - 0.111_323_079_4).abs() < 1e-9);
    assert_eq!(model.n_bonds(), 5);
    assert_eq!(model.n_angles(), 6);
    assert_eq!(model.n_torsions(), 4);
    assert_eq!(model.n_inversions(), 6);
    let (_, forces) = model.energy_forces().unwrap();
    let step = 1e-4;
    let mut maximum = 0.0_f64;
    for atom in 0..model.n_atoms() {
        for axis in 0..3 {
            let mut displaced = model.positions().to_vec();
            displaced[atom][axis] += step;
            let up = model.energy_at(&displaced).unwrap();
            displaced[atom][axis] -= 2.0 * step;
            let down = model.energy_at(&displaced).unwrap();
            let finite_difference = -(up - down) / (2.0 * step);
            maximum = maximum.max((finite_difference - forces[atom][axis]).abs());
        }
    }
    assert!(maximum < 1e-4, "maximum force error {maximum:e} eV/A");
}

#[test]
fn benzene_is_resonant() {
    let mut atoms: Vec<(&str, [f64; 3])> = Vec::new();
    for index in 0..6 {
        let angle = std::f64::consts::TAU * (index as f64) / 6.0;
        atoms.push(("C", [1.39 * angle.cos(), 1.39 * angle.sin(), 0.0]));
    }
    for index in 0..6 {
        let angle = std::f64::consts::TAU * (index as f64) / 6.0;
        atoms.push(("H", [2.47 * angle.cos(), 2.47 * angle.sin(), 0.0]));
    }
    let model = UffModel::build(&molecule("benzene", &atoms), options(8.0)).unwrap();
    assert_eq!(model.typing().aromatic_rings.len(), 1);
    assert_eq!(model.typing().aromatic_rings[0].len(), 6);
    for atom_type in model.types().iter().take(6) {
        assert_eq!(*atom_type, "C_R");
    }
    assert!(model.energy().unwrap().is_finite());
}

#[test]
fn registry_and_run_record_round_trip() {
    let registry = CalculatorRegistry::builtin();
    assert_eq!(registry.names(), vec!["uff"]);
    let info = registry.info("uff").unwrap();
    assert!(info.available());
    assert!(info.capabilities.periodic && info.capabilities.forces);
    assert!(!info.description.is_empty());

    let calculator = registry.build("uff", &options(8.0)).unwrap();
    assert_eq!(calculator.name(), "uff");
    assert!(calculator.provides_forces());
    assert!(calculator.provides_stress());
    let energy = calculator.energy(&water()).unwrap();
    let (energy_forces, forces) = calculator.energy_forces(&water()).unwrap();
    assert!((energy - energy_forces).abs() < 1e-12);
    assert_eq!(forces.len(), 3);
    let stress = calculator.stress(&water()).unwrap().unwrap();
    assert!(stress[0][0].is_finite() && stress[1][1].is_finite());

    let model = UffModel::build(&water(), options(8.0)).unwrap();
    let record = RunRecord::new("uff", energy, model.n_atoms())
        .with_note("water at the UFF equilibrium angle")
        .with_steps(3, false);
    let text = record.to_json().unwrap();
    assert!(text.contains("\"calculator\":\"uff\""));
    let back = RunRecord::from_json(&text).unwrap();
    assert_eq!(back, record);
    assert_eq!(back.n_atoms, 3);
    assert_eq!(back.steps, 3);
    assert!(!back.converged);

    let unknown = match registry.build("not-an-engine", &options(8.0)) {
        Ok(_) => panic!("an unregistered engine was built"),
        Err(error) => error,
    };
    assert_eq!(unknown.code(), "unknown_calculator");
}

#[test]
fn calculator_trait_and_model_agree() {
    let calculator = UffCalculator::new(options(12.0)).unwrap();
    let model = calculator.model(&diamond()).unwrap();
    let (energy, forces) = calculator.energy_forces(&diamond()).unwrap();
    assert!((energy - model.energy().unwrap()).abs() < 1e-12);
    assert_eq!(forces.len(), 8);
    assert!(calculator.energy(&diamond()).unwrap().is_finite());
}
