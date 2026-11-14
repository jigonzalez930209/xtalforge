//! Generation verbs: checkpointed structure search (Phase 20).

use serde_json::Value;

use crate::answers::Answer;
use crate::run::{integer_or, optional_string};
use crate::{Error, Session};

/// Runs one chunk of a random or evolutionary search.
pub(crate) fn search(_session: &mut Session, arguments: &Value) -> Result<Answer, Error> {
    let kind = optional_string(arguments, "kind")?.unwrap_or("random");
    let kind = match kind {
        "random" => xtalforge_generation::SearchKind::Random,
        "evolutionary" => xtalforge_generation::SearchKind::Evolutionary,
        other => {
            return Err(Error::Argument {
                name: "kind".to_owned(),
                reason: format!("`{other}` is not random or evolutionary"),
            });
        }
    };
    let options = arguments
        .get("options")
        .cloned()
        .ok_or_else(|| Error::Argument {
            name: "options".to_owned(),
            reason: "is required".to_owned(),
        })?;
    let total = match kind {
        xtalforge_generation::SearchKind::Random => options
            .get("n_structures")
            .and_then(Value::as_u64)
            .unwrap_or(1) as usize,
        xtalforge_generation::SearchKind::Evolutionary => options
            .get("generations")
            .and_then(Value::as_u64)
            .unwrap_or(1) as usize,
    };
    let chunk = integer_or(arguments, "chunk", 1)?;
    let checkpoint = optional_string(arguments, "checkpoint")?;
    let result = if let Some(path) = checkpoint {
        xtalforge_generation::run_checkpointed(kind, options, total, path, chunk as usize)
            .map_err(Error::Core)?
    } else {
        let mut checkpoint = xtalforge_generation::SearchCheckpoint::new(kind, options, total);
        let result = match kind {
            xtalforge_generation::SearchKind::Random => {
                xtalforge_generation::run_random_chunk(&mut checkpoint, chunk as usize)
            }
            xtalforge_generation::SearchKind::Evolutionary => {
                xtalforge_generation::run_evolutionary_chunk(&mut checkpoint, chunk as usize)
            }
        };
        result.map_err(Error::Core)?
    };
    Ok(Answer::Search { chunk: result })
}
