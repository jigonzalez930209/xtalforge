//! Example plugin for xtalforge (Task 8.3.1).
//!
//! A separate crate can extend every open registry without editing the
//! library:
//!
//! - a file format ([`FormatRegistry`]);
//! - an energy calculator ([`CalculatorRegistry`]);
//! - a module ([`ModuleRegistry`]);
//! - a named scalar ([`ScalarRegistry`]);
//! - a view-style preset ([`StyleRegistry`]).
//!
//! The example is deliberately small: a two-line format, a zero-energy
//! calculator, a site-count module, a scalar that returns the first fractional
//! coordinate and a large-atom style.
//!
//! The draw-style enum itself is closed. A renderer must agree on the geometry
//! it draws; the open part is the settings preset, not the enum.

use std::sync::OnceLock;

use serde_json::{Value, json};
use xtalforge::commands::Document;
use xtalforge::core::{Error, Lattice, Site, Structure, StructureBuilder};
use xtalforge::ff::{Calculator, CalculatorInfo, CalculatorOptions, CalculatorRegistry};
use xtalforge::modules::{Module, ModuleInfo, ModuleRegistry};
use xtalforge::view::ViewSettings;

/// Name of the example format.
pub const FORMAT_NAME: &str = "two-line";

/// Reads `cell a` followed by `element x y z` lines in a cubic cell.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a malformed line or an unknown element.
pub fn read_two_line(text: &str) -> Result<Structure, Error> {
    let mut lines = text.lines().filter(|line| {
        let line = line.trim();
        !line.is_empty() && !line.starts_with('#')
    });
    let Some(header) = lines.next() else {
        return Err(Error::InvalidInput("two-line: empty input".to_owned()));
    };
    let Some(size) = header.trim().strip_prefix("cell ") else {
        return Err(Error::InvalidInput(
            "two-line: the first line must be `cell <a>`".to_owned(),
        ));
    };
    let size: f64 = size
        .parse()
        .map_err(|_| Error::InvalidInput("two-line: bad cell length".to_owned()))?;
    let lattice = Lattice::cubic(size)?;
    let mut builder = StructureBuilder::new("two-line", lattice);
    for line in lines {
        let parts: Vec<&str> = line.split_whitespace().collect();
        let [element, x, y, z] = parts.as_slice() else {
            return Err(Error::InvalidInput(
                "two-line: every atom line is `element x y z`".to_owned(),
            ));
        };
        let parse = |value: &str, axis: &str| {
            value
                .parse::<f64>()
                .map_err(|_| Error::InvalidInput(format!("two-line: bad {axis}")))
        };
        let point = [parse(x, "x")?, parse(y, "y")?, parse(z, "z")?];
        builder = builder.site(Site::new(*element, element, point)?);
    }
    builder.build()
}

/// Writes the two-line format.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a non-cubic cell.
pub fn write_two_line(structure: &Structure) -> Result<String, Error> {
    let parameters = structure.lattice().parameters();
    if (parameters[1] - parameters[0]).abs() > 1e-9
        || (parameters[2] - parameters[0]).abs() > 1e-9
        || parameters[3..]
            .iter()
            .any(|angle| (angle - 90.0).abs() > 1e-9)
    {
        return Err(Error::InvalidInput(
            "two-line: the format only stores a cubic cell".to_owned(),
        ));
    }
    let mut out = format!("cell {:.6}\n", parameters[0]);
    for site in structure.sites() {
        let [x, y, z] = site.frac();
        out.push_str(&format!("{} {x:.8} {y:.8} {z:.8}\n", site.element()));
    }
    Ok(out)
}

static ZERO_INFO: OnceLock<CalculatorInfo> = OnceLock::new();

fn zero_info() -> &'static CalculatorInfo {
    ZERO_INFO.get_or_init(|| CalculatorInfo {
        name: "zero".to_owned(),
        label: "Zero energy".to_owned(),
        description: "Example plugin: returns zero energy and zero forces.".to_owned(),
        availability: xtalforge::ff::Availability::Available,
        capabilities: xtalforge::ff::CalculatorCapabilities {
            periodic: true,
            forces: true,
            stress: false,
            charges: false,
            in_process: true,
            thread_safe: true,
        },
    })
}

/// A calculator that reports zero for everything.
#[derive(Debug, Default)]
pub struct ZeroCalculator;

impl Calculator for ZeroCalculator {
    fn name(&self) -> &'static str {
        "zero"
    }

    fn energy(&self, _structure: &Structure) -> Result<f64, xtalforge::ff::Error> {
        Ok(0.0)
    }

    fn energy_forces(
        &self,
        structure: &Structure,
    ) -> Result<(f64, Vec<[f64; 3]>), xtalforge::ff::Error> {
        let atoms = structure.p1()?.atoms.len();
        Ok((0.0, vec![[0.0; 3]; atoms]))
    }

    fn provides_forces(&self) -> bool {
        true
    }
}

fn build_zero(_options: CalculatorOptions) -> Result<Box<dyn Calculator>, xtalforge::ff::Error> {
    Ok(Box::new(ZeroCalculator))
}

/// Module that counts sites.
#[derive(Debug, Default)]
pub struct SiteCountModule;

impl Module for SiteCountModule {
    fn info(&self) -> ModuleInfo {
        ModuleInfo {
            name: "site-count".to_owned(),
            description: "Example plugin: counts the asymmetric sites.".to_owned(),
            version: env!("CARGO_PKG_VERSION").to_owned(),
            parameters: Vec::new(),
        }
    }

    fn run(&self, parameters: &Value) -> Result<Value, Error> {
        Ok(json!({
            "plugin": "site-count",
            "received": parameters,
        }))
    }
}

/// Scalar that returns the first fractional coordinate of every site.
///
/// # Errors
///
/// Never fails; the signature matches [`xtalforge::core::ScalarFn`].
pub fn first_fractional(structure: &Structure) -> Result<Vec<f64>, Error> {
    Ok(structure
        .sites()
        .iter()
        .map(|site| site.frac()[0])
        .collect())
}

/// A large-atom style.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an empty name or a duplicate.
pub fn large_atoms(
    registry: &mut xtalforge::view::StyleRegistry,
) -> Result<(), xtalforge::core::Error> {
    let settings = ViewSettings {
        atom_scale: 2.5,
        bond_radius: 0.2,
        ..ViewSettings::default()
    };
    registry.register("large-atoms", settings)
}

/// Registers every plugin in this crate.
///
/// # Errors
///
/// Returns the first registry error (duplicate name or invalid format).
pub fn register_all(
    formats: &mut xtalforge::io::registry::FormatRegistry,
    calculators: &mut CalculatorRegistry,
    modules: &mut ModuleRegistry,
    scalars: &mut xtalforge::core::ScalarRegistry,
    styles: &mut xtalforge::view::StyleRegistry,
) -> Result<(), Error> {
    formats.register(xtalforge::io::registry::FormatSpec {
        name: FORMAT_NAME,
        extensions: &[".2l"],
        read: Some(read_two_line),
        write: Some(write_two_line),
        capabilities: xtalforge::io::registry::Capabilities {
            cell: true,
            labels: true,
            ..xtalforge::io::registry::Capabilities::default()
        },
    })?;
    calculators
        .register(zero_info(), build_zero)
        .map_err(|error| Error::InvalidInput(error.to_string()))?;
    modules
        .register(Box::new(SiteCountModule))
        .map_err(|error| Error::InvalidInput(error.to_string()))?;
    scalars
        .register(
            xtalforge::core::ScalarInfo {
                name: "first_fractional".to_owned(),
                description: "First fractional coordinate of each site".to_owned(),
                unit: "".to_owned(),
            },
            first_fractional,
        )
        .map_err(|error| Error::InvalidInput(error.to_string()))?;
    large_atoms(styles)?;
    Ok(())
}

/// A document built from the example format, used by the test.
///
/// # Errors
///
/// Returns the reader, builder or command error.
pub fn example_document() -> Result<Document, Error> {
    let structure = read_two_line("cell 5.0\nNa 0 0 0\nCl 0.5 0.5 0.5\n")?;
    Ok(Document::new(structure))
}
