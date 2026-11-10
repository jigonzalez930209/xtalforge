//! Document verbs: load, save, info, diagnostics, edits and history.

use std::path::Path;

use serde_json::Value;
use xtalforge_commands::{Origin, Supercell, WrapIntoCell};

use crate::answers::{Answer, HistoryEntry};
use crate::run::{document_argument, optional_string, required_string};
use crate::session::entry_summary;
use crate::{Error, Session, changes, diagnostics, inspect};

pub(crate) fn load(session: &mut Session, arguments: &Value) -> Result<Answer, Error> {
    let path = optional_string(arguments, "path")?;
    let text = optional_string(arguments, "text")?;
    let (content, source) = match (path, text) {
        (Some(path), None) => {
            let content = std::fs::read_to_string(path)
                .map_err(|error| Error::Io(format!("{path}: {error}")))?;
            (content, Some(path.to_owned()))
        }
        (None, Some(text)) => (text.to_owned(), None),
        (Some(_), Some(_)) => {
            return Err(Error::Argument {
                name: "path/text".to_owned(),
                reason: "give exactly one of them".to_owned(),
            });
        }
        (None, None) => {
            return Err(Error::Argument {
                name: "path".to_owned(),
                reason: "a file path or the text content is required".to_owned(),
            });
        }
    };
    let hint = optional_string(arguments, "format")?
        .map(str::to_owned)
        .or_else(|| {
            path.and_then(|path| {
                Path::new(path)
                    .extension()
                    .and_then(|extension| extension.to_str())
                    .map(str::to_owned)
            })
        })
        .or_else(|| xtalforge_io::registry::sniff(&content).map(|format| format.name.to_owned()))
        .ok_or_else(|| Error::Argument {
            name: "format".to_owned(),
            reason: "could not detect the format from the content; give one".to_owned(),
        })?;
    let structure = xtalforge_io::registry::read_with(&content, &hint)?;
    let name = match optional_string(arguments, "name")? {
        Some(name) => {
            if name.trim().is_empty() {
                return Err(Error::Argument {
                    name: "name".to_owned(),
                    reason: "must not be empty".to_owned(),
                });
            }
            session.insert(name, structure, source.clone());
            name.to_owned()
        }
        None => {
            let name = session.available_name(path, "structure");
            session.insert(name.clone(), structure, source.clone());
            name
        }
    };
    let document = session.document(Some(&name))?;
    let summary = entry_summary(&name, document, source.as_deref())?;
    Ok(Answer::Loaded { document: summary })
}

pub(crate) fn save(session: &mut Session, arguments: &Value) -> Result<Answer, Error> {
    let name = session.resolve(document_argument(arguments)?)?;
    let path = required_string(arguments, "path")?.to_owned();
    let hint = optional_string(arguments, "format")?
        .map(str::to_owned)
        .or_else(|| {
            Path::new(&path)
                .extension()
                .and_then(|extension| extension.to_str())
                .map(str::to_owned)
        })
        .ok_or_else(|| Error::Argument {
            name: "format".to_owned(),
            reason: "the path has no extension; give a format".to_owned(),
        })?;
    let format = xtalforge_io::registry::find(&hint).ok_or_else(|| Error::Argument {
        name: "format".to_owned(),
        reason: format!("unknown format {hint}"),
    })?;
    let document = session.document_mut(Some(&name))?;
    let text = xtalforge_io::registry::write_with(document.structure(), format.name)?;
    std::fs::write(&path, text.as_bytes())
        .map_err(|error| Error::Io(format!("{path}: {error}")))?;
    document.mark_clean();
    Ok(Answer::Saved {
        document: name,
        path,
        format: format.name.to_owned(),
        bytes: text.len(),
    })
}

pub(crate) fn info(session: &mut Session, arguments: &Value) -> Result<Answer, Error> {
    let name = session.resolve(document_argument(arguments)?)?;
    let document = session.document(Some(&name))?;
    let inspection = inspect::inspect_structure(document.structure())?;
    let report = diagnostics::diagnose(document.structure())?;
    Ok(Answer::Info {
        document: name,
        inspection,
        diagnostics: report,
    })
}

pub(crate) fn diagnostics(session: &mut Session, arguments: &Value) -> Result<Answer, Error> {
    let name = session.resolve(document_argument(arguments)?)?;
    let document = session.document(Some(&name))?;
    let report = diagnostics::diagnose(document.structure())?;
    Ok(Answer::Diagnostics {
        document: name,
        report,
    })
}

pub(crate) fn supercell(session: &mut Session, arguments: &Value) -> Result<Answer, Error> {
    let name = session.resolve(document_argument(arguments)?)?;
    let counts = crate::run::repeats(arguments)?;
    let document = session.document_mut(Some(&name))?;
    document.apply(
        Box::new(Supercell::new(counts[0], counts[1], counts[2])),
        Origin::Agent,
    )?;
    structure_answer(session, &name)
}

pub(crate) fn wrap(session: &mut Session, arguments: &Value) -> Result<Answer, Error> {
    let name = session.resolve(document_argument(arguments)?)?;
    let document = session.document_mut(Some(&name))?;
    document.apply(Box::new(WrapIntoCell), Origin::Agent)?;
    structure_answer(session, &name)
}

pub(crate) fn structure_answer(session: &Session, name: &str) -> Result<Answer, Error> {
    let document = session.document(Some(name))?;
    let summary = entry_summary(name, document, session.source(name))?;
    let description = document
        .history()
        .last()
        .map(|record| record.description.clone())
        .unwrap_or_default();
    Ok(Answer::Structure {
        document: name.to_owned(),
        summary,
        description,
    })
}

pub(crate) fn undo(session: &mut Session, arguments: &Value, undo: bool) -> Result<Answer, Error> {
    let name = session.resolve(document_argument(arguments)?)?;
    let document = session.document_mut(Some(&name))?;
    let applied = if undo {
        document.undo()
    } else {
        document.redo()
    };
    history_answer(session, &name, applied)
}

pub(crate) fn history(session: &mut Session, arguments: &Value) -> Result<Answer, Error> {
    let name = session.resolve(document_argument(arguments)?)?;
    history_answer(session, &name, false)
}

fn history_answer(session: &Session, name: &str, applied: bool) -> Result<Answer, Error> {
    let document = session.document(Some(name))?;
    let entries = document
        .history()
        .iter()
        .map(|record| HistoryEntry {
            description: record.description.clone(),
            origin: match record.origin {
                Origin::User => "user",
                Origin::Script => "script",
                Origin::Agent => "agent",
            }
            .to_owned(),
        })
        .collect();
    Ok(Answer::History {
        document: name.to_owned(),
        applied,
        entries,
        can_undo: document.can_undo(),
        can_redo: document.can_redo(),
    })
}

pub(crate) fn changes(session: &mut Session, arguments: &Value) -> Result<Answer, Error> {
    let name = session.resolve(document_argument(arguments)?)?;
    let document = session.document(Some(&name))?;
    let bursts = changes::changes_for(document.history());
    Ok(Answer::Changes {
        document: name,
        bursts,
    })
}

pub(crate) fn close(session: &mut Session, arguments: &Value) -> Result<Answer, Error> {
    let name = session.resolve(document_argument(arguments)?)?;
    session.remove(&name);
    Ok(session.documents_answer())
}
