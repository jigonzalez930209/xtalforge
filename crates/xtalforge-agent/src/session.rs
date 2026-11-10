//! Stateful session (Task 8.1.2).
//!
//! A session owns documents, runs verbs by name with JSON arguments and
//! returns typed answers. Every mutation goes through
//! [`xtalforge_commands::Document`] with [`xtalforge_commands::Origin::Agent`], so the history and
//! the change bursts are the same objects a user front end sees.

use std::collections::BTreeMap;
use std::path::Path;

use serde_json::Value;
use xtalforge_commands::Document;
use xtalforge_core::{Limits, Profile, Structure};

use crate::Error;
use crate::answers::{Answer, DocumentSummary};

/// Options a session starts with.
#[derive(Debug, Clone, PartialEq)]
pub struct SessionOptions {
    /// Interactive size limits applied to every verb.
    pub limits: Limits,
}

impl Default for SessionOptions {
    fn default() -> Self {
        Self {
            limits: Limits::new(Profile::Standard),
        }
    }
}

#[derive(Debug)]
pub(crate) struct Entry {
    document: Document,
    source: Option<String>,
}

/// A stateful agent session.
#[derive(Debug)]
pub struct Session {
    documents: BTreeMap<String, Entry>,
    active: Option<String>,
    options: SessionOptions,
}

impl Default for Session {
    fn default() -> Self {
        Self::new()
    }
}

impl Session {
    /// Creates an empty session with the interactive limits.
    #[must_use]
    pub fn new() -> Self {
        Self::with_options(SessionOptions::default())
    }

    /// Creates an empty session with explicit options.
    #[must_use]
    pub fn with_options(options: SessionOptions) -> Self {
        Self {
            documents: BTreeMap::new(),
            active: None,
            options,
        }
    }

    /// Options the session runs under.
    #[must_use]
    pub fn options(&self) -> &SessionOptions {
        &self.options
    }

    /// Name of the active document.
    #[must_use]
    pub fn active_name(&self) -> Option<&str> {
        self.active.as_deref()
    }

    /// The structure of one document, or of the active one.
    ///
    /// # Errors
    ///
    /// Returns [`Error::UnknownDocument`] or [`Error::NoDocument`].
    pub fn structure(&self, name: Option<&str>) -> Result<&Structure, Error> {
        Ok(self.document(name)?.structure())
    }

    /// The active document, when there is one.
    #[must_use]
    pub fn active(&self) -> Option<&Document> {
        self.active_name()
            .and_then(|name| self.documents.get(name))
            .map(|entry| &entry.document)
    }

    /// Whether a document is open.
    #[must_use]
    pub fn contains(&self, name: &str) -> bool {
        self.documents.contains_key(name)
    }

    /// Number of open documents.
    #[must_use]
    pub fn len(&self) -> usize {
        self.documents.len()
    }

    /// Whether no document is open.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.documents.is_empty()
    }

    /// Resolves a document name, falling back to the active document.
    ///
    /// # Errors
    ///
    /// Returns [`Error::UnknownDocument`] or [`Error::NoDocument`].
    pub(crate) fn resolve(&self, name: Option<&str>) -> Result<String, Error> {
        match name {
            Some(name) => {
                if self.documents.contains_key(name) {
                    Ok(name.to_owned())
                } else {
                    Err(Error::UnknownDocument(name.to_owned()))
                }
            }
            None => self
                .active_name()
                .map(str::to_owned)
                .ok_or(Error::NoDocument),
        }
    }

    pub(crate) fn document(&self, name: Option<&str>) -> Result<&Document, Error> {
        let name = self.resolve(name)?;
        self.documents
            .get(&name)
            .map(|entry| &entry.document)
            .ok_or(Error::UnknownDocument(name))
    }

    pub(crate) fn document_mut(&mut self, name: Option<&str>) -> Result<&mut Document, Error> {
        let name = self.resolve(name)?;
        self.documents
            .get_mut(&name)
            .map(|entry| &mut entry.document)
            .ok_or(Error::UnknownDocument(name))
    }

    pub(crate) fn source(&self, name: &str) -> Option<&str> {
        self.documents
            .get(name)
            .and_then(|entry| entry.source.as_deref())
    }

    /// Adds a document; an existing name is replaced and made active.
    pub(crate) fn insert(
        &mut self,
        name: impl Into<String>,
        structure: Structure,
        source: Option<String>,
    ) -> String {
        let name = name.into();
        self.documents.insert(
            name.clone(),
            Entry {
                document: Document::new(structure),
                source,
            },
        );
        self.active = Some(name.clone());
        name
    }

    /// Suggests a free document name from a path or a fallback.
    pub(crate) fn available_name(&self, path: Option<&str>, fallback: &str) -> String {
        let stem = path
            .and_then(|path| Path::new(path).file_stem())
            .and_then(|stem| stem.to_str())
            .filter(|stem| !stem.is_empty())
            .unwrap_or(fallback);
        if !self.documents.contains_key(stem) {
            return stem.to_owned();
        }
        let mut index = 2;
        loop {
            let candidate = format!("{stem}-{index}");
            if !self.documents.contains_key(&candidate) {
                return candidate;
            }
            index += 1;
        }
    }

    pub(crate) fn remove(&mut self, name: &str) -> bool {
        let removed = self.documents.remove(name).is_some();
        if removed && self.active.as_deref() == Some(name) {
            self.active = self.documents.keys().next().cloned();
        }
        removed
    }

    /// Summaries of the open documents, in name order.
    #[must_use]
    pub fn summaries(&self) -> Vec<DocumentSummary> {
        self.documents
            .iter()
            .filter_map(|(name, entry)| summary(name, entry).ok())
            .collect()
    }

    /// The `documents` answer.
    #[must_use]
    pub fn documents_answer(&self) -> Answer {
        Answer::Documents {
            documents: self.summaries(),
            active: self.active.clone(),
        }
    }

    /// Runs one verb.
    ///
    /// # Errors
    ///
    /// Returns [`Error::UnknownVerb`] for a name outside the catalogue, the
    /// argument errors of the verb and the errors of the operation it runs.
    pub fn run(&mut self, verb: &str, arguments: &Value) -> Result<Answer, Error> {
        if crate::verbs::lookup(verb).is_none() {
            return Err(Error::UnknownVerb(verb.to_owned()));
        }
        crate::run::run(self, verb, arguments)
    }
}

pub(crate) fn summary(name: &str, entry: &Entry) -> Result<DocumentSummary, Error> {
    let structure = entry.document.structure();
    let summary = structure.summary()?;
    Ok(DocumentSummary {
        name: name.to_owned(),
        title: summary.name,
        formula: summary.formula,
        n_sites: summary.n_sites,
        volume: summary.volume,
        density: summary.density,
        source: entry.source.clone(),
        dirty: entry.document.dirty(),
    })
}

pub(crate) fn entry_summary(
    name: &str,
    document: &Document,
    source: Option<&str>,
) -> Result<DocumentSummary, Error> {
    let summary = document.structure().summary()?;
    Ok(DocumentSummary {
        name: name.to_owned(),
        title: summary.name,
        formula: summary.formula,
        n_sites: summary.n_sites,
        volume: summary.volume,
        density: summary.density,
        source: source.map(str::to_owned),
        dirty: document.dirty(),
    })
}
