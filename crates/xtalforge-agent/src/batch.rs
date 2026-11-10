//! Stateless batch runner (Task 8.1.5).
//!
//! `run_verbs` executes a JSON script against a fresh session. It is the
//! convenience form the Python and WebAssembly packages expose as
//! `agent_run`, and the same input the CLI `agent` command accepts.

use serde_json::Value;

use crate::{Answer, Error, Session};

/// Runs a JSON script and returns one answer per call.
///
/// The input is either an array of `{"verb": ..., "arguments": {...}}`
/// objects or an object with a `verbs` array. Every call shares one session,
/// so a later verb sees the documents and edits of an earlier one.
///
/// # Errors
///
/// Returns [`Error::Json`] for invalid input and the verb errors unchanged.
pub fn run_verbs(script: &str) -> Result<Vec<Answer>, Error> {
    let document: Value =
        serde_json::from_str(script).map_err(|error| Error::Json(error.to_string()))?;
    let calls = document
        .get("verbs")
        .and_then(Value::as_array)
        .or_else(|| document.as_array())
        .ok_or_else(|| Error::Argument {
            name: "script".to_owned(),
            reason: "must be an array of {verb, arguments} objects".to_owned(),
        })?;
    let mut session = Session::new();
    let mut answers = Vec::with_capacity(calls.len());
    for call in calls {
        let verb = call
            .get("verb")
            .and_then(Value::as_str)
            .ok_or_else(|| Error::Argument {
                name: "verb".to_owned(),
                reason: "every script entry needs a string `verb`".to_owned(),
            })?;
        let arguments = call
            .get("arguments")
            .cloned()
            .unwrap_or_else(|| Value::Object(serde_json::Map::new()));
        answers.push(session.run(verb, &arguments)?);
    }
    Ok(answers)
}
