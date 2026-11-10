//! Force-field verbs: hydrogens, energy and optimization.

use serde_json::Value;
use xtalforge_build::{HydrogenOptions, add_hydrogens};
use xtalforge_commands::{Origin, ReplaceStructure};
use xtalforge_ff::{
    CalculatorOptions, Minimizer, OptimizeOptions, UffCalculator, UffModel,
    apply_optimized_geometry, minimize,
};

use crate::answers::{Answer, OptimizationSummary};
use crate::ops::structure_answer;
use crate::run::{document_argument, integer_or, number_or, optional_string};
use crate::{Error, Session};

fn add_hydrogens_command(
    structure: &xtalforge_core::Structure,
    distance: f64,
) -> Result<xtalforge_core::Structure, Error> {
    let options = HydrogenOptions {
        distance,
        ..Default::default()
    };
    Ok(add_hydrogens(structure, &options)?)
}

pub(crate) fn hydrogens(session: &mut Session, arguments: &Value) -> Result<Answer, Error> {
    let name = session.resolve(document_argument(arguments)?)?;
    let distance = number_or(arguments, "distance", 1.09)?;
    let document = session.document_mut(Some(&name))?;
    let replacement = add_hydrogens_command(document.structure(), distance)?;
    document.apply(
        Box::new(ReplaceStructure::with_description(
            replacement,
            format!("Add hydrogens at {distance:.2} A"),
        )),
        Origin::Agent,
    )?;
    structure_answer(session, &name)
}

fn max_force(forces: &[[f64; 3]]) -> f64 {
    forces
        .iter()
        .map(|force| (force[0] * force[0] + force[1] * force[1] + force[2] * force[2]).sqrt())
        .fold(0.0_f64, f64::max)
}

pub(crate) fn energy(session: &mut Session, arguments: &Value) -> Result<Answer, Error> {
    let name = session.resolve(document_argument(arguments)?)?;
    let document = session.document(Some(&name))?;
    let structure = document.structure();
    let model = UffModel::build(structure, CalculatorOptions::default())?;
    let (energy, forces) = model.energy_forces()?;
    let p1 = structure.p1()?;
    let stress = model.stress().ok();
    Ok(Answer::Energy {
        document: name,
        energy_ev: energy,
        n_atoms: p1.atoms.len(),
        max_force_ev_per_angstrom: max_force(&forces),
        stress_gpa: stress,
    })
}

pub(crate) fn optimize(session: &mut Session, arguments: &Value) -> Result<Answer, Error> {
    let name = session.resolve(document_argument(arguments)?)?;
    let max_steps = integer_or(arguments, "max_steps", 200)?;
    let max_steps = usize::try_from(max_steps).map_err(|_| Error::Argument {
        name: "max_steps".to_owned(),
        reason: "does not fit in usize".to_owned(),
    })?;
    let algorithm = optional_string(arguments, "algorithm")?
        .unwrap_or("fire")
        .trim()
        .to_ascii_lowercase();
    let minimizer = match algorithm.as_str() {
        "fire" => Minimizer::Fire,
        "lbfgs" | "l-bfgs" => Minimizer::Lbfgs,
        other => {
            return Err(Error::Argument {
                name: "algorithm".to_owned(),
                reason: format!("unknown minimizer {other:?}; use fire or lbfgs"),
            });
        }
    };
    let calculator = UffCalculator::with_defaults()?;
    let document = session.document_mut(Some(&name))?;
    let options = OptimizeOptions {
        max_steps,
        ..Default::default()
    };
    let result = minimize(&calculator, document.structure(), &options, minimizer)?;
    let replacement = apply_optimized_geometry(&result, document.structure())?;
    let summary = OptimizationSummary {
        energy_ev: result.energy,
        converged: result.converged,
        steps: result.steps,
        max_force_ev_per_angstrom: result.max_force,
    };
    document.apply(
        Box::new(ReplaceStructure::with_description(
            replacement,
            format!("UFF optimization ({algorithm}, {} steps)", result.steps),
        )),
        Origin::Agent,
    )?;
    let new_structure = document.structure().summary()?;
    Ok(Answer::Optimization {
        document: name,
        optimization: summary,
        summary: crate::inspect::StructureSummaryView::from(new_structure),
    })
}
