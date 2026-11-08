//! Built-in parameter descriptions shared by every surface (Task 3.8.2).
//!
//! The sets below are the single definition used by the CLI, the Python
//! bindings, the WASM bindings and the JSON schema the agent layer reads. A
//! front end renders them; it never duplicates the bounds or help text.

use serde_json::json;

use crate::params::{ParamSet, ParamSpec, ParamType};

/// Names of the built-in parameter sets.
pub const PRESET_NAMES: [&str; 3] = ["symmetry", "pxrd", "uff"];

fn spec(
    name: &str,
    kind: ParamType,
    default: serde_json::Value,
    help: &str,
    group: &str,
) -> ParamSpec {
    ParamSpec {
        name: name.to_owned(),
        kind,
        default,
        help: help.to_owned(),
        required: false,
        group: Some(group.to_owned()),
    }
}

/// Returns the named built-in parameter set.
///
/// # Errors
///
/// Returns [`crate::Error::InvalidInput`] for an unknown name or an invalid
/// built-in default.
pub fn preset(name: &str) -> Result<ParamSet, crate::Error> {
    let set = ParamSet::new();
    let set = match name {
        "symmetry" => set
            .with_param(spec(
                "symprec",
                ParamType::Float {
                    min: Some(1e-6),
                    max: Some(1.0),
                    unit: Some("Angstrom".to_owned()),
                },
                json!(1e-3),
                "Position tolerance for the space-group search",
                "Symmetry",
            ))?
            .with_param(spec(
                "angle_tolerance",
                ParamType::Float {
                    min: Some(0.1),
                    max: Some(10.0),
                    unit: Some("degree".to_owned()),
                },
                json!(5.0),
                "Angle tolerance for the space-group search",
                "Symmetry",
            ))?
            .with_param(spec(
                "hall_number",
                ParamType::Integer {
                    min: Some(1),
                    max: Some(530),
                },
                json!(1),
                "Preferred Hall setting; 1 means the spglib default",
                "Symmetry",
            ))?,
        "pxrd" => set
            .with_param(spec(
                "radiation",
                ParamType::Choice {
                    choices: [
                        "cu-ka", "cu-ka1", "mo-ka", "co-ka", "fe-ka", "cr-ka", "ag-ka",
                    ]
                    .map(str::to_owned)
                    .to_vec(),
                },
                json!("cu-ka"),
                "X-ray source",
                "Pattern",
            ))?
            .with_param(spec(
                "max_two_theta",
                ParamType::Float {
                    min: Some(5.0),
                    max: Some(180.0),
                    unit: Some("degree".to_owned()),
                },
                json!(80.0),
                "Highest two-theta in the pattern",
                "Pattern",
            ))?
            .with_param(spec(
                "step",
                ParamType::Float {
                    min: Some(0.001),
                    max: Some(1.0),
                    unit: Some("degree".to_owned()),
                },
                json!(0.02),
                "Pattern sampling step",
                "Pattern",
            ))?
            .with_param(spec(
                "fwhm",
                ParamType::Float {
                    min: Some(0.01),
                    max: Some(5.0),
                    unit: Some("degree".to_owned()),
                },
                json!(0.1),
                "Peak full width at half maximum",
                "Profile",
            ))?,
        "uff" => set
            .with_param(spec(
                "parameter_set",
                ParamType::Choice {
                    choices: ["uff", "uff4mof", "uff4mof-ii"].map(str::to_owned).to_vec(),
                },
                json!("uff4mof"),
                "Parameter table",
                "Force field",
            ))?
            .with_param(spec(
                "vdw",
                ParamType::Boolean,
                json!(true),
                "Whether the Lennard-Jones term is evaluated",
                "Force field",
            ))?
            .with_param(spec(
                "vdw_cutoff",
                ParamType::Float {
                    min: Some(1.0),
                    max: Some(50.0),
                    unit: Some("Angstrom".to_owned()),
                },
                json!(12.0),
                "van der Waals cutoff",
                "Force field",
            ))?
            .with_param(spec(
                "total_charge",
                ParamType::Float {
                    min: Some(-20.0),
                    max: Some(20.0),
                    unit: Some("e".to_owned()),
                },
                json!(0.0),
                "Total charge for the charge equilibration",
                "Charges",
            ))?,
        other => {
            return Err(crate::Error::InvalidInput(format!(
                "unknown parameter preset {other:?}; expected one of {PRESET_NAMES:?}"
            )));
        }
    };
    Ok(set)
}

/// Every built-in set with its name, in [`PRESET_NAMES`] order.
///
/// # Errors
///
/// Returns the errors of [`preset`].
pub fn presets() -> Result<Vec<(&'static str, ParamSet)>, crate::Error> {
    PRESET_NAMES
        .iter()
        .map(|name| preset(name).map(|set| (*name, set)))
        .collect()
}
