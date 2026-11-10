//! Verb dispatch and argument parsing.
//!
//! The catalogue in [`crate::verbs`] declares the names; this module turns a
//! JSON object into the arguments each implementation expects. Every missing
//! or mistyped argument is reported as [`crate::Error::Argument`], never as a
//! panic.

use serde_json::Value;

use crate::{Answer, Error, Session};

pub(crate) fn run(session: &mut Session, verb: &str, arguments: &Value) -> Result<Answer, Error> {
    if !arguments.is_object() && verb != "capabilities" && verb != "list_documents" {
        return Err(Error::Argument {
            name: "arguments".to_owned(),
            reason: "must be a JSON object".to_owned(),
        });
    }
    match verb {
        "capabilities" => Ok(Answer::Capabilities {
            capabilities: crate::capability_report(),
        }),
        "list_documents" => Ok(session.documents_answer()),
        "load" => crate::ops::load(session, arguments),
        "save" => crate::ops::save(session, arguments),
        "info" => crate::ops::info(session, arguments),
        "diagnostics" => crate::ops::diagnostics(session, arguments),
        "supercell" => crate::ops::supercell(session, arguments),
        "wrap" => crate::ops::wrap(session, arguments),
        "hydrogens" => crate::ops_ff::hydrogens(session, arguments),
        "symmetry" => crate::ops_analysis::symmetry(session, arguments),
        "bonds" => crate::ops_analysis::bonds(session, arguments),
        "rings" => crate::ops_analysis::rings(session, arguments),
        "pxrd" => crate::ops_analysis::pxrd(session, arguments),
        "porosity" => crate::ops_analysis::porosity(session, arguments),
        "uff_energy" => crate::ops_ff::energy(session, arguments),
        "optimize" => crate::ops_ff::optimize(session, arguments),
        "render" => crate::ops_picture::render(session, arguments),
        "undo" => crate::ops::undo(session, arguments, true),
        "redo" => crate::ops::undo(session, arguments, false),
        "history" => crate::ops::history(session, arguments),
        "changes" => crate::ops::changes(session, arguments),
        "close" => crate::ops::close(session, arguments),
        _ => Err(Error::UnknownVerb(verb.to_owned())),
    }
}

pub(crate) fn optional_string<'a>(
    arguments: &'a Value,
    name: &str,
) -> Result<Option<&'a str>, Error> {
    match arguments.get(name) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::String(value)) => Ok(Some(value)),
        Some(_) => Err(Error::Argument {
            name: name.to_owned(),
            reason: "must be a string".to_owned(),
        }),
    }
}

pub(crate) fn required_string<'a>(arguments: &'a Value, name: &str) -> Result<&'a str, Error> {
    optional_string(arguments, name)?.ok_or_else(|| Error::Argument {
        name: name.to_owned(),
        reason: "is required".to_owned(),
    })
}

pub(crate) fn document_argument(arguments: &Value) -> Result<Option<&str>, Error> {
    optional_string(arguments, "document")
}

pub(crate) fn optional_number(arguments: &Value, name: &str) -> Result<Option<f64>, Error> {
    match arguments.get(name) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(number)) => number.as_f64().map(Some).ok_or_else(|| Error::Argument {
            name: name.to_owned(),
            reason: "is not representable as a number".to_owned(),
        }),
        Some(_) => Err(Error::Argument {
            name: name.to_owned(),
            reason: "must be a number".to_owned(),
        }),
    }
}

pub(crate) fn optional_integer(arguments: &Value, name: &str) -> Result<Option<u64>, Error> {
    match arguments.get(name) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(number)) => number.as_u64().map(Some).ok_or_else(|| Error::Argument {
            name: name.to_owned(),
            reason: "must be a non-negative integer".to_owned(),
        }),
        Some(_) => Err(Error::Argument {
            name: name.to_owned(),
            reason: "must be a non-negative integer".to_owned(),
        }),
    }
}

pub(crate) fn number_or(arguments: &Value, name: &str, default: f64) -> Result<f64, Error> {
    Ok(optional_number(arguments, name)?.unwrap_or(default))
}

pub(crate) fn integer_or(arguments: &Value, name: &str, default: u64) -> Result<u64, Error> {
    Ok(optional_integer(arguments, name)?.unwrap_or(default))
}

pub(crate) fn repeats(arguments: &Value) -> Result<[u32; 3], Error> {
    let Some(value) = arguments.get("repeats") else {
        return Err(Error::Argument {
            name: "repeats".to_owned(),
            reason: "is required".to_owned(),
        });
    };
    let Some(items) = value.as_array() else {
        return Err(Error::Argument {
            name: "repeats".to_owned(),
            reason: "must be an array of three integers".to_owned(),
        });
    };
    if items.len() != 3 {
        return Err(Error::Argument {
            name: "repeats".to_owned(),
            reason: "must contain exactly three integers".to_owned(),
        });
    }
    let mut counts = [1_u32; 3];
    for (slot, item) in counts.iter_mut().zip(items.iter()) {
        let Some(count) = item.as_u64() else {
            return Err(Error::Argument {
                name: "repeats".to_owned(),
                reason: "every entry must be a positive integer".to_owned(),
            });
        };
        *slot = u32::try_from(count).map_err(|_| Error::Argument {
            name: "repeats".to_owned(),
            reason: "every entry must fit in 32 bits".to_owned(),
        })?;
        if *slot == 0 {
            return Err(Error::Argument {
                name: "repeats".to_owned(),
                reason: "every entry must be at least 1".to_owned(),
            });
        }
    }
    Ok(counts)
}
