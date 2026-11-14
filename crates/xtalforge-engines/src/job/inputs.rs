//! Input files written into a compute job directory (Task 18.1/18.4).

use std::collections::BTreeSet;

use serde::{Deserialize, Serialize};

use crate::Error;
use crate::input_set::InputFile;

use super::util::{validate_relative_path, write_atomic};

/// The input files written into the job scratch directory.
///
/// The files come from [`crate::input_set::write_input_set`] through
/// [`super::ComputeStep::request`]. Writers that are not registered yet (or
/// surfaces that have already materialized files) can pass them directly with
/// [`super::ComputeStep::request_with_inputs`] or [`InputFiles::new`]; that is
/// the documented integration point for Task 18.1.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct InputFiles {
    /// Named files, in input-set order.
    pub files: Vec<InputFile>,
}

impl InputFiles {
    /// Wraps an explicit file list.
    #[must_use]
    pub fn new(files: Vec<InputFile>) -> Self {
        Self { files }
    }

    /// One file.
    #[must_use]
    pub fn one(name: impl Into<String>, text: impl Into<String>) -> Self {
        Self {
            files: vec![InputFile {
                name: name.into(),
                text: text.into(),
            }],
        }
    }

    /// File names in input-set order.
    #[must_use]
    pub fn names(&self) -> Vec<&str> {
        self.files.iter().map(|file| file.name.as_str()).collect()
    }

    /// Rejects duplicate, absolute and escaping names.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Job`] for an unsafe or duplicate file name.
    pub fn validate(&self) -> Result<(), Error> {
        let mut seen = BTreeSet::new();
        for file in &self.files {
            validate_relative_path(&file.name, "input file name")?;
            if !seen.insert(file.name.as_str()) {
                return Err(Error::Job(format!(
                    "input file {:?} appears more than once",
                    file.name
                )));
            }
        }
        Ok(())
    }

    /// Writes every file below `directory`, creating parents as needed.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Job`] for an unsafe name or a filesystem failure.
    pub fn write_to(&self, directory: &std::path::Path) -> Result<(), Error> {
        self.validate()?;
        for file in &self.files {
            write_atomic(&directory.join(&file.name), &file.text)?;
        }
        Ok(())
    }
}
