//! ML verbs: prediction, selection and metadata append (Phase 21).

use serde_json::Value;
use xtalforge_commands::{Origin, ReplaceStructure};
use xtalforge_ml::{
    AcsfOptions, MbtrOptions, Predictor, PropertyModel, SelectionOptions, SelectionStrategy, acsf,
    mbtr, select_batch, soap_power_spectrum,
};

use crate::answers::Answer;
use crate::run::{document_argument, number_or, optional_string, required_string};
use crate::{Error, Session};

fn model(arguments: &Value) -> Result<PropertyModel, Error> {
    let text = required_string(arguments, "model")?;
    let json = if text.trim_start().starts_with('{') {
        text.to_owned()
    } else {
        std::fs::read_to_string(text)
            .map_err(|error| Error::Io(format!("could not read {text}: {error}")))?
    };
    PropertyModel::from_json(&json).map_err(Error::Core)
}

fn features(structure: &xtalforge_core::Structure, kind: &str) -> Result<Vec<f64>, Error> {
    match kind {
        "soap" => {
            let descriptor = soap_power_spectrum(structure, &Default::default())?;
            descriptor
                .rows
                .first()
                .cloned()
                .ok_or_else(|| Error::Argument {
                    name: "features".to_owned(),
                    reason: "SOAP produced no row".to_owned(),
                })
        }
        "acsf" => {
            let descriptor = acsf(
                structure,
                &AcsfOptions {
                    average: true,
                    ..AcsfOptions::default()
                },
            )?;
            descriptor
                .rows
                .first()
                .cloned()
                .ok_or_else(|| Error::Argument {
                    name: "features".to_owned(),
                    reason: "ACSF produced no row".to_owned(),
                })
        }
        "mbtr" => Ok(mbtr(structure, &MbtrOptions::default())?.features),
        other => Err(Error::Argument {
            name: "features".to_owned(),
            reason: format!("`{other}` is not soap, acsf or mbtr"),
        }),
    }
}

/// Predicts a scalar property for the active document.
pub(crate) fn predict(session: &mut Session, arguments: &Value) -> Result<Answer, Error> {
    let name = session.resolve(document_argument(arguments)?)?;
    let model = model(arguments)?;
    let kind = optional_string(arguments, "features")?.unwrap_or("soap");
    let structure = session.document(Some(&name))?.structure();
    let input = features(structure, kind)?;
    let value = model.predict(&input)?;
    Ok(Answer::MlPrediction {
        document: name,
        model: model.name.clone(),
        value,
    })
}

/// Selects a batch from a feature table.
pub(crate) fn select(_session: &mut Session, arguments: &Value) -> Result<Answer, Error> {
    let features: Vec<Vec<f64>> =
        serde_json::from_value(arguments.get("features").cloned().ok_or_else(|| {
            Error::Argument {
                name: "features".to_owned(),
                reason: "is required".to_owned(),
            }
        })?)
        .map_err(|error| Error::Argument {
            name: "features".to_owned(),
            reason: error.to_string(),
        })?;
    let uncertainties: Vec<f64> = match arguments.get("uncertainties") {
        None | Some(Value::Null) => vec![0.0; features.len()],
        Some(value) => serde_json::from_value(value.clone()).map_err(|error| Error::Argument {
            name: "uncertainties".to_owned(),
            reason: error.to_string(),
        })?,
    };
    let strategy = match optional_string(arguments, "strategy")?.unwrap_or("uncertainty") {
        "uncertainty" => SelectionStrategy::Uncertainty,
        "farthest_point" => SelectionStrategy::FarthestPoint,
        "diversity" => SelectionStrategy::Diversity,
        other => {
            return Err(Error::Argument {
                name: "strategy".to_owned(),
                reason: format!("`{other}` is not a selection strategy"),
            });
        }
    };
    let batch = crate::run::integer_or(arguments, "batch", 8)? as usize;
    let selected = select_batch(
        &features,
        &uncertainties,
        &SelectionOptions {
            batch,
            strategy,
            diversity_radius: 0.0,
        },
    )?;
    Ok(Answer::MlSelection { selected })
}

/// Appends a numeric metadata property to the active document.
pub(crate) fn append(session: &mut Session, arguments: &Value) -> Result<Answer, Error> {
    let name = session.resolve(document_argument(arguments)?)?;
    let key = required_string(arguments, "key")?.to_owned();
    let value = number_or(arguments, "value", 0.0)?;
    let document = session.document_mut(Some(&name))?;
    let structure = document
        .structure()
        .clone()
        .with_metadata(key.clone(), serde_json::json!(value));
    document.apply(
        Box::new(ReplaceStructure::with_description(
            structure,
            format!("Append metadata `{key}` = {value}"),
        )),
        Origin::Agent,
    )?;
    Ok(Answer::MlAppend {
        document: name,
        key,
        value,
    })
}
