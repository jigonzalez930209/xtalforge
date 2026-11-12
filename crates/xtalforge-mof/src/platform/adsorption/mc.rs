//! Widom insertion, GCMC, energy grids and RASPA export.

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Structure};

use super::{AVOGADRO, KB_EV, force_field, guest_by_name, guest_energy};

fn splitmix(seed: u64) -> impl FnMut() -> f64 {
    let mut state = seed;
    move || {
        state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        (z ^ (z >> 31)) as f64 / u64::MAX as f64
    }
}

/// Widom insertion options.
#[derive(Debug, Clone, PartialEq)]
pub struct WidomOptions {
    /// Temperature in K.
    pub temperature: f64,
    /// Number of random insertions.
    pub samples: usize,
    /// Deterministic seed.
    pub seed: u64,
}

impl Default for WidomOptions {
    fn default() -> Self {
        Self {
            temperature: 298.0,
            samples: 10_000,
            seed: 0,
        }
    }
}

/// Widom insertion result.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct WidomResult {
    /// Guest name.
    pub guest: String,
    /// Temperature in K.
    pub temperature: f64,
    /// Number of insertions attempted.
    pub samples: usize,
    /// Mean Boltzmann factor (dimensionless).
    pub henry: f64,
    /// Mean interaction energy weighted by the Boltzmann factor in eV.
    pub mean_energy: f64,
    /// Heat of adsorption in eV (positive for exothermic).
    pub heat_of_adsorption: f64,
    /// Lowest energy found in eV.
    pub minimum_energy: f64,
}

/// Widom insertion with a deterministic seed.
///
/// The Henry coefficient is reported as the dimensionless average Boltzmann
/// factor `<exp(-E/kT)>`; multiply by the cell free volume and the pressure
/// conversion for a dimensional isotherm.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for invalid options or a missing guest.
pub fn widom_insertion(
    structure: &Structure,
    guest_name: &str,
    options: &WidomOptions,
) -> Result<WidomResult, Error> {
    if !options.temperature.is_finite() || options.temperature <= 0.0 {
        return Err(Error::InvalidInput(
            "temperature must be positive and finite".to_owned(),
        ));
    }
    if options.samples == 0 {
        return Err(Error::InvalidInput(
            "Widom insertion needs at least one sample".to_owned(),
        ));
    }
    let force_field = force_field(structure)?;
    let guest = guest_by_name(&force_field, guest_name)?;
    let mut random = splitmix(options.seed);
    let beta = 1.0 / (KB_EV * options.temperature);
    let mut boltzmann = 0.0;
    let mut weighted_energy = 0.0;
    let mut minimum = f64::INFINITY;
    for _ in 0..options.samples {
        let frac = [random(), random(), random()];
        let angles = [
            random() * std::f64::consts::TAU,
            random() * std::f64::consts::TAU,
            random() * std::f64::consts::TAU,
        ];
        let energy = guest_energy(&force_field, guest, frac, angles);
        let factor = (-energy * beta).exp();
        boltzmann += factor;
        weighted_energy += factor * energy;
        minimum = minimum.min(energy);
    }
    let henry = boltzmann / options.samples as f64;
    let mean_energy = if boltzmann > 0.0 {
        weighted_energy / boltzmann
    } else {
        minimum
    };
    Ok(WidomResult {
        guest: guest.name.clone(),
        temperature: options.temperature,
        samples: options.samples,
        henry,
        mean_energy,
        heat_of_adsorption: -mean_energy,
        minimum_energy: minimum,
    })
}

/// GCMC options for a single component.
#[derive(Debug, Clone, PartialEq)]
pub struct GcmcOptions {
    /// Temperature in K.
    pub temperature: f64,
    /// Pressure points in bar.
    pub pressures: Vec<f64>,
    /// Monte Carlo steps per pressure.
    pub steps: usize,
    /// Deterministic seed.
    pub seed: u64,
}

/// One isotherm point.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct IsothermPoint {
    /// Pressure in bar.
    pub pressure: f64,
    /// Loading in molecules per cell.
    pub loading: f64,
    /// Acceptance ratio.
    pub acceptance: f64,
}

/// A single-component GCMC isotherm.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Isotherm {
    /// Guest name.
    pub guest: String,
    /// Temperature in K.
    pub temperature: f64,
    /// Isotherm points.
    pub points: Vec<IsothermPoint>,
}

/// Grand-canonical Monte Carlo with a rigid host and a single rigid guest.
///
/// Insertions and deletions use the host-guest energy only; guest-guest
/// interactions are not included. The chemical potential is the ideal-gas
/// `kT ln(P/P0)` with `P0 = 1 bar`, so the result is a Henry/Langmuir
/// isotherm and not a calibrated GCMC replacement.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for invalid options or a missing guest.
pub fn gcmc_isotherm(
    structure: &Structure,
    guest_name: &str,
    options: &GcmcOptions,
) -> Result<Isotherm, Error> {
    if !options.temperature.is_finite() || options.temperature <= 0.0 {
        return Err(Error::InvalidInput(
            "temperature must be positive and finite".to_owned(),
        ));
    }
    if options.pressures.is_empty() || options.steps == 0 {
        return Err(Error::InvalidInput(
            "GCMC needs pressures and steps".to_owned(),
        ));
    }
    let force_field = force_field(structure)?;
    let guest = guest_by_name(&force_field, guest_name)?;
    let beta = 1.0 / (KB_EV * options.temperature);
    let mut random = splitmix(options.seed);
    let mut points = Vec::new();
    for &pressure in &options.pressures {
        if !pressure.is_finite() || pressure <= 0.0 {
            return Err(Error::InvalidInput(
                "pressures must be positive and finite".to_owned(),
            ));
        }
        let mu = (pressure.ln()) / beta;
        let mut molecules = 0usize;
        let mut accepted = 0usize;
        let mut sum = 0.0;
        for _ in 0..options.steps {
            let frac = [random(), random(), random()];
            let angles = [
                random() * std::f64::consts::TAU,
                random() * std::f64::consts::TAU,
                random() * std::f64::consts::TAU,
            ];
            let energy = guest_energy(&force_field, guest, frac, angles);
            let activity = (mu - energy) * beta;
            if activity > 0.0 || random() < activity.exp().min(1.0) {
                molecules += 1;
                accepted += 1;
                sum += energy;
            } else if molecules > 0 && random() < 0.05 {
                molecules -= 1;
            }
        }
        points.push(IsothermPoint {
            pressure,
            loading: molecules as f64 / options.steps as f64 * AVOGADRO / 1e24,
            acceptance: accepted as f64 / options.steps as f64,
        });
        let _ = sum;
    }
    Ok(Isotherm {
        guest: guest.name.clone(),
        temperature: options.temperature,
        points,
    })
}

/// A precomputed guest energy grid.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct EnergyGrid {
    /// Grid samples in eV.
    pub grid: xtalforge_io::grid::Grid3,
    /// Guest name.
    pub guest: String,
    /// Minimum energy in eV.
    pub minimum: f64,
    /// Maximum energy in eV.
    pub maximum: f64,
}

/// Builds a guest energy grid for fast visualisation and insertion.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a bad spacing or missing guest, or the
/// errors of the force field and grid construction.
pub fn energy_grid(
    structure: &Structure,
    guest_name: &str,
    spacing: f64,
) -> Result<EnergyGrid, Error> {
    if !spacing.is_finite() || spacing <= 0.0 {
        return Err(Error::InvalidInput(
            "the grid spacing must be positive and finite".to_owned(),
        ));
    }
    let force_field = force_field(structure)?;
    let guest = guest_by_name(&force_field, guest_name)?;
    let lengths = [
        vector_length(force_field.matrix[0]),
        vector_length(force_field.matrix[1]),
        vector_length(force_field.matrix[2]),
    ];
    let shape = lengths.map(|length| (length / spacing).round().max(1.0) as usize);
    let step = [
        scale(force_field.matrix[0], 1.0 / shape[0] as f64),
        scale(force_field.matrix[1], 1.0 / shape[1] as f64),
        scale(force_field.matrix[2], 1.0 / shape[2] as f64),
    ];
    let mut data = Vec::with_capacity(shape[0] * shape[1] * shape[2]);
    let mut minimum = f64::INFINITY;
    let mut maximum = f64::NEG_INFINITY;
    for z in 0..shape[2] {
        for y in 0..shape[1] {
            for x in 0..shape[0] {
                let frac = [
                    x as f64 / shape[0] as f64,
                    y as f64 / shape[1] as f64,
                    z as f64 / shape[2] as f64,
                ];
                let energy = guest_energy(&force_field, guest, frac, [0.0, 0.0, 0.0]);
                minimum = minimum.min(energy);
                maximum = maximum.max(energy);
                data.push(energy);
            }
        }
    }
    let grid = xtalforge_io::grid::Grid3::new([0.0; 3], step, shape, data, true)?;
    Ok(EnergyGrid {
        grid,
        guest: guest.name.clone(),
        minimum,
        maximum,
    })
}

/// RASPA export files.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RaspaFiles {
    /// `simulation.input`.
    pub simulation_input: String,
    /// Framework CIF.
    pub framework_cif: String,
    /// Pseudo-atom force-field definitions.
    pub force_field: String,
    /// Mixing rules.
    pub mixing_rules: String,
}

/// Writes the RASPA input set for a Widom or GCMC run.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for invalid options, or the errors of the
/// force field and CIF writer.
pub fn raspa_input(
    structure: &Structure,
    guest_name: &str,
    temperature: f64,
    pressure: f64,
) -> Result<RaspaFiles, Error> {
    if !temperature.is_finite() || temperature <= 0.0 || !pressure.is_finite() || pressure <= 0.0 {
        return Err(Error::InvalidInput(
            "RASPA needs a positive temperature and pressure".to_owned(),
        ));
    }
    let force_field = force_field(structure)?;
    let guest = guest_by_name(&force_field, guest_name)?;
    let simulation_input = format!(
        "SimulationType                MonteCarlo\n\
NumberOfCycles               10000\n\
NumberOfInitializationCycles 1000\n\
PrintEvery                  1000\n\
Forcefield                  UFF-TraPPE\n\
FrameworkName               framework\n\
FrameworkNumberofUnitCells  1 1 1\n\
ChargeMethod                Equil\n\
Temperature                 {temperature}\n\
Pressure                    {pressure}\n\
Component 0 MoleculeName    {}\n\
           MoleculeDefinition Local\n\
           WidomInsertion    yes\n",
        guest.name
    );
    let framework_cif = xtalforge_io::write_structure(structure, "cif")?;
    let mut force_field_text = String::new();
    for host in &force_field.hosts {
        force_field_text.push_str(&format!(
            "{} lj {:.4} {:.4} charge {:.4}\n",
            host.element, host.sigma, host.epsilon_k, host.charge
        ));
    }
    for atom in &guest.atoms {
        force_field_text.push_str(&format!(
            "{} lj {:.4} {:.4} charge {:.4}\n",
            atom.element, atom.sigma, atom.epsilon_k, atom.charge
        ));
    }
    Ok(RaspaFiles {
        simulation_input,
        framework_cif,
        force_field: force_field_text,
        mixing_rules: "Lorentz-Berthelot\n".to_owned(),
    })
}

fn vector_length(vector: [f64; 3]) -> f64 {
    (vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]).sqrt()
}

fn scale(vector: [f64; 3], factor: f64) -> [f64; 3] {
    [vector[0] * factor, vector[1] * factor, vector[2] * factor]
}
