//! The velocity-Verlet loop shared by the four ensembles.

use xtalforge_core::{Lattice, Structure};

use crate::Error;
use crate::calculator::Calculator;
use crate::math::{add, dot, scale, sub};
use crate::md::barostat::Barostat;
use crate::md::random::SeededRng;
use crate::md::sampling::{Recorder, SampleInput};
use crate::md::thermostat::Thermostat;
use crate::md::{AMU_TO_EV_FS2_PER_A2, KB_EV_PER_K, MdOptions, MdResult, degrees_of_freedom};
use crate::run::{ProgressEvent, RunControl};
use crate::structure_util::with_frac;

/// `Angstrom/fs^2` per `(eV/Angstrom)/amu`.
const ACCELERATION: f64 = 1.0 / AMU_TO_EV_FS2_PER_A2;

/// `GPa` per `eV/Angstrom^3`.
const GPA_PER_EV_PER_A3: f64 = 160.217_663_4;

/// One evaluated state of the system.
struct State {
    structure: Structure,
    positions: Vec<[f64; 3]>,
    forces: Vec<[f64; 3]>,
    energy: f64,
    lattice: Lattice,
}

/// Runs the integration loop.
pub(crate) fn run(
    calculator: &dyn Calculator,
    structure: &Structure,
    options: &MdOptions,
    control: &mut RunControl<'_>,
) -> Result<MdResult, Error> {
    options.validate()?;
    if !structure.is_p1() {
        return Err(Error::InvalidInput(
            "molecular dynamics needs a P1 structure; expand it with Structure::as_p1 first"
                .to_owned(),
        ));
    }
    let n_atoms = structure.sites().len();
    if n_atoms == 0 {
        return Err(Error::InvalidInput(
            "molecular dynamics needs at least one atom".to_owned(),
        ));
    }
    let masses = masses_of(structure)?;
    let lattice = *structure.lattice();
    if lattice.volume().abs() <= 0.0 {
        return Err(Error::InvalidInput(
            "molecular dynamics needs a cell with positive volume".to_owned(),
        ));
    }

    let mut rng = SeededRng::new(options.seed);
    let mut positions = structure.cart_sites();
    let mut velocities = initial_velocities(&masses, options.initial_temperature, &mut rng);
    remove_com_momentum(&masses, &mut velocities);
    let dof = degrees_of_freedom(n_atoms);
    let mut state = evaluate(calculator, structure, lattice, &positions)?;
    let mut kinetic = kinetic_energy(&velocities, &masses);
    let initial_energy = state.energy + kinetic;

    let mut thermostat = Thermostat::new(&options.ensemble, n_atoms);
    let barostat = Barostat::new(&options.ensemble);
    let mut recorder = Recorder::new(options);
    let dt = options.time_step_fs;
    let mut steps = 0_usize;
    let mut cancelled = false;

    if recorder.sample_due(0, options.steps) {
        record_step(&mut recorder, control, 0, options, &state, kinetic, dof);
    }

    for step in 1..=options.steps {
        if control.is_cancelled() {
            cancelled = true;
            break;
        }
        half_kick(&mut velocities, &state.forces, &masses, 0.5 * dt);
        thermostat.pre(&mut velocities, &masses, dt, &mut rng);
        for (position, velocity) in positions.iter_mut().zip(&velocities) {
            *position = add(*position, scale(*velocity, dt));
        }
        if options.wrap {
            wrap_positions(&mut positions, &state.lattice);
        }
        state = evaluate(calculator, structure, state.lattice, &positions)?;
        half_kick(&mut velocities, &state.forces, &masses, 0.5 * dt);
        kinetic = kinetic_energy(&velocities, &masses);
        thermostat.post(&mut velocities, dt, kinetic);

        if let Some(barostat) = &barostat
            && step % barostat.interval() == 0
        {
            let volume = state.lattice.volume().abs();
            let proposal = barostat.propose(&mut rng, volume, state.energy);
            let new_lattice = state.lattice.scaled_to_volume(proposal.new_volume)?;
            let factor = proposal.factor.cbrt();
            let candidate_positions: Vec<[f64; 3]> = state
                .positions
                .iter()
                .map(|position| scale(*position, factor))
                .collect();
            let candidate = evaluate(calculator, structure, new_lattice, &candidate_positions)?;
            if barostat.accepted(&mut rng, &proposal, candidate.energy, n_atoms) {
                state = candidate;
                positions = candidate_positions;
                for velocity in velocities.iter_mut() {
                    *velocity = scale(*velocity, factor);
                }
            }
        }

        kinetic = kinetic_energy(&velocities, &masses);
        steps = step;
        if recorder.sample_due(step, options.steps) {
            record_step(&mut recorder, control, step, options, &state, kinetic, dof);
        }
    }

    let (mean_temperature, mean_pressure_gpa) = recorder.means();
    let samples = recorder.samples().to_vec();
    let trajectory = recorder.into_trajectory();
    Ok(MdResult {
        structure: state.structure,
        ensemble: options.ensemble.clone(),
        samples,
        initial_energy,
        final_energy: state.energy + kinetic,
        mean_temperature,
        mean_pressure_gpa,
        steps,
        cancelled,
        trajectory,
    })
}

/// Records one sample, its frame and its progress event.
fn record_step(
    recorder: &mut Recorder,
    control: &mut RunControl<'_>,
    step: usize,
    options: &MdOptions,
    state: &State,
    kinetic: f64,
    dof: f64,
) {
    let volume = state.lattice.volume().abs();
    let pressure_gpa =
        pressure_ev_per_a3(&state.positions, &state.forces, kinetic, volume) * GPA_PER_EV_PER_A3;
    let temperature = 2.0 * kinetic / (dof * KB_EV_PER_K);
    let max_force = state
        .forces
        .iter()
        .map(|force| crate::math::norm(*force))
        .fold(0.0_f64, f64::max);
    let frac: Vec<[f64; 3]> = state
        .positions
        .iter()
        .map(|position| state.lattice.to_frac(*position))
        .collect();
    let sample = recorder.record(
        SampleInput {
            step,
            time_fs: step as f64 * options.time_step_fs,
            potential: state.energy,
            kinetic,
            temperature,
            volume,
            pressure_gpa,
            max_force,
            lattice: state.lattice,
        },
        &frac,
    );
    control.report(ProgressEvent::Dynamics {
        step,
        potential: sample.potential_energy,
        kinetic: sample.kinetic_energy,
        temperature: sample.temperature,
    });
}

/// Builds the structure at `positions` and evaluates energy and forces.
fn evaluate(
    calculator: &dyn Calculator,
    base: &Structure,
    lattice: Lattice,
    positions: &[[f64; 3]],
) -> Result<State, Error> {
    let frac: Vec<[f64; 3]> = positions
        .iter()
        .map(|position| lattice.to_frac(*position))
        .collect();
    let structure = with_frac(base, lattice, &frac)?;
    let (energy, forces) = calculator.energy_forces(&structure)?;
    if forces.len() != positions.len() {
        return Err(Error::InvalidInput(format!(
            "the calculator returned {} force rows for {} atoms",
            forces.len(),
            positions.len()
        )));
    }
    Ok(State {
        structure,
        positions: positions.to_vec(),
        forces,
        energy,
        lattice,
    })
}

/// Atomic masses of a P1 structure, in amu.
fn masses_of(structure: &Structure) -> Result<Vec<f64>, Error> {
    let mut masses = Vec::with_capacity(structure.sites().len());
    for site in structure.sites() {
        let mass = xtalforge_core::elements::mass_of(site.element()).map_err(|error| {
            Error::InvalidInput(format!(
                "molecular dynamics needs a mass for site {} ({}): {error}",
                site.label(),
                site.element()
            ))
        })?;
        if !mass.is_finite() || mass <= 0.0 {
            return Err(Error::InvalidInput(format!(
                "site {} has an unusable mass {mass}",
                site.label()
            )));
        }
        masses.push(mass);
    }
    Ok(masses)
}

/// Maxwell-Boltzmann velocities at `temperature`; zero at zero temperature.
fn initial_velocities(masses: &[f64], temperature: f64, rng: &mut SeededRng) -> Vec<[f64; 3]> {
    let mut velocities = vec![[0.0; 3]; masses.len()];
    if temperature <= 0.0 {
        return velocities;
    }
    for (velocity, mass) in velocities.iter_mut().zip(masses) {
        let sigma = (KB_EV_PER_K * temperature / (AMU_TO_EV_FS2_PER_A2 * mass)).sqrt();
        for value in velocity.iter_mut() {
            *value = sigma * rng.normal();
        }
    }
    velocities
}

/// Removes the centre-of-mass momentum.
fn remove_com_momentum(masses: &[f64], velocities: &mut [[f64; 3]]) {
    let total: f64 = masses.iter().sum();
    if total <= 0.0 {
        return;
    }
    let mut momentum = [0.0; 3];
    for (velocity, mass) in velocities.iter().zip(masses) {
        momentum = add(momentum, scale(*velocity, *mass));
    }
    for velocity in velocities.iter_mut() {
        *velocity = sub(*velocity, scale(momentum, 1.0 / total));
    }
}

/// Kinetic energy in eV.
fn kinetic_energy(velocities: &[[f64; 3]], masses: &[f64]) -> f64 {
    let mut total = 0.0;
    for (velocity, mass) in velocities.iter().zip(masses) {
        total += mass * dot(*velocity, *velocity);
    }
    0.5 * total * AMU_TO_EV_FS2_PER_A2
}

/// Half a velocity kick with the current forces.
fn half_kick(velocities: &mut [[f64; 3]], forces: &[[f64; 3]], masses: &[f64], half_dt: f64) {
    let factor = half_dt * ACCELERATION;
    for ((velocity, force), mass) in velocities.iter_mut().zip(forces).zip(masses) {
        if *mass <= 0.0 {
            continue;
        }
        *velocity = add(*velocity, scale(*force, factor / *mass));
    }
}

/// Wraps Cartesian positions into the cell.
fn wrap_positions(positions: &mut [[f64; 3]], lattice: &Lattice) {
    for position in positions.iter_mut() {
        let frac = Lattice::wrap(lattice.to_frac(*position));
        *position = lattice.to_cart(frac);
    }
}

/// Molecular virial pressure in eV/Angstrom^3.
fn pressure_ev_per_a3(
    positions: &[[f64; 3]],
    forces: &[[f64; 3]],
    kinetic: f64,
    volume: f64,
) -> f64 {
    if volume <= 0.0 {
        return 0.0;
    }
    let virial: f64 = positions
        .iter()
        .zip(forces)
        .map(|(position, force)| dot(*position, *force))
        .sum();
    (2.0 * kinetic + virial) / (3.0 * volume)
}
