//! RASPA Widom-insertion input set (Task 18.1.6).
//!
//! This is a byte-compatible adapter for the existing
//! `xtalforge_mof::platform::raspa_input` export. It rebuilds the same four
//! files from the same UFF host parameters, EQeq charges and TraPPE guest
//! table; the default preset reproduces the Phase 10 output exactly.

use std::collections::BTreeMap;
use std::sync::OnceLock;

use serde::Deserialize;
use xtalforge_core::{P1Atom, Structure};

use crate::Error;
use crate::input_set::settings::Settings;
use crate::input_set::util;

const SIMULATION: &str = "simulation.input";
const FRAMEWORK: &str = "framework.cif";
const FORCE_FIELD: &str = "force_field";
const MIXING: &str = "mixing_rules";

/// Boltzmann constant in eV/K, the value used by the Phase 10 export.
const KB_EV: f64 = 8.617_333_262e-5;
/// Conversion from kcal/mol to eV, the value used by the Phase 10 export.
const KCAL_PER_MOL_TO_EV: f64 = 0.043_364_1;

const GUESTS: &str = include_str!("../../data/input-sets/raspa-guests.json");

#[derive(Debug, Deserialize)]
struct GuestFile {
    guests: Vec<Guest>,
}

#[derive(Debug, Clone, Deserialize)]
struct Guest {
    name: String,
    atoms: Vec<GuestAtom>,
}

#[derive(Debug, Clone, Deserialize)]
struct GuestAtom {
    element: String,
    sigma: f64,
    epsilon_k: f64,
    charge: f64,
}

#[derive(Debug, Clone)]
struct Host {
    element: String,
    sigma: f64,
    epsilon_k: f64,
    charge: f64,
}

static GUEST_TABLE: OnceLock<Result<BTreeMap<String, Guest>, String>> = OnceLock::new();

fn guests() -> Result<&'static BTreeMap<String, Guest>, Error> {
    match GUEST_TABLE.get_or_init(load_guests) {
        Ok(table) => Ok(table),
        Err(error) => Err(Error::Unavailable(format!("RASPA guest table: {error}"))),
    }
}

fn load_guests() -> Result<BTreeMap<String, Guest>, String> {
    let file: GuestFile =
        serde_json::from_str(GUESTS).map_err(|error| format!("invalid JSON: {error}"))?;
    Ok(file
        .guests
        .into_iter()
        .map(|guest| (guest.name.to_ascii_uppercase(), guest))
        .collect())
}

pub(crate) fn build(
    structure: &Structure,
    settings: &Settings,
) -> Result<Vec<super::InputFile>, Error> {
    let p1 = structure.p1()?;
    if p1.atoms.is_empty() {
        return Err(Error::Job(
            "RASPA needs a structure with at least one atom".to_owned(),
        ));
    }
    let guest_name = settings.text("guest")?;
    let table = guests()?;
    let guest = table
        .get(&guest_name.to_ascii_uppercase())
        .ok_or_else(|| Error::Job(format!("unknown RASPA guest {guest_name:?}")))?;
    let temperature = settings.float("temperature")?;
    let pressure = settings.float("pressure")?;
    if !temperature.is_finite() || temperature <= 0.0 || !pressure.is_finite() || pressure <= 0.0 {
        return Err(Error::Job(
            "RASPA needs a positive temperature and pressure".to_owned(),
        ));
    }
    let cycles = settings.integer("cycles")?;
    let initialization = settings.integer("initialization_cycles")?;
    let print_every = settings.integer("print_every")?;
    let forcefield = util::token(settings.text("forcefield")?, "force field")?;
    let simulation = format!(
        "SimulationType                MonteCarlo\n\
NumberOfCycles               {cycles}\n\
NumberOfInitializationCycles {initialization}\n\
PrintEvery                  {print_every}\n\
Forcefield                  {forcefield}\n\
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
    let framework = xtalforge_io::write_structure(structure, "cif")?;
    let hosts = host_parameters(structure)?;
    let mut force_field = String::new();
    for host in &hosts {
        force_field.push_str(&format!(
            "{} lj {:.4} {:.4} charge {:.4}\n",
            host.element, host.sigma, host.epsilon_k, host.charge
        ));
    }
    for atom in &guest.atoms {
        force_field.push_str(&format!(
            "{} lj {:.4} {:.4} charge {:.4}\n",
            atom.element, atom.sigma, atom.epsilon_k, atom.charge
        ));
    }
    Ok(vec![
        super::InputFile {
            name: SIMULATION.to_owned(),
            text: simulation,
        },
        super::InputFile {
            name: FRAMEWORK.to_owned(),
            text: framework,
        },
        super::InputFile {
            name: FORCE_FIELD.to_owned(),
            text: force_field,
        },
        super::InputFile {
            name: MIXING.to_owned(),
            text: "Lorentz-Berthelot\n".to_owned(),
        },
    ])
}

fn host_parameters(structure: &Structure) -> Result<Vec<Host>, Error> {
    let p1 = structure.p1()?;
    let snapshot = util::p1_snapshot(structure)?;
    let charges = xtalforge_ff::eqeq_charges(&snapshot, &xtalforge_ff::ChargeOptions::eqeq())
        .map(|charges| charges.per_atom)
        .unwrap_or_else(|_| vec![0.0; p1.atoms.len()]);
    let mut hosts = Vec::with_capacity(p1.atoms.len());
    for (index, atom) in p1.atoms.iter().enumerate() {
        let params = first_uff(atom)?;
        hosts.push(Host {
            element: atom.element.clone(),
            sigma: params.x1,
            epsilon_k: params.d1 * KCAL_PER_MOL_TO_EV / KB_EV,
            charge: charges.get(index).copied().unwrap_or(0.0),
        });
    }
    Ok(hosts)
}

fn first_uff(atom: &P1Atom) -> Result<&'static xtalforge_ff::uff::params::UffParams, Error> {
    xtalforge_ff::uff::params::types_for(&atom.element)
        .map_err(util::force_field_error)?
        .into_iter()
        .next()
        .ok_or_else(|| {
            Error::Job(format!(
                "no UFF parameters for host element {}",
                atom.element
            ))
        })
}
