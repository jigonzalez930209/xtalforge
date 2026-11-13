//! Skill documents shipped with the CLI (Task 8.2.3).
//!
//! The files under `skill/` are compiled into the crate. `xtalforge skill
//! show` prints one, `xtalforge skill install` writes them into a directory an
//! assistant can read, and `xtalforge skill path` prints the default place.
//! The same material is published as the automation guide.

use std::path::{Path, PathBuf};

use crate::Error;

/// One embedded skill document.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SkillDocument {
    /// Relative path, for example `references/porosity.md`.
    pub name: &'static str,
    /// One-line title.
    pub title: &'static str,
    /// Markdown body.
    pub body: &'static str,
}

macro_rules! document {
    ($name:literal, $title:literal) => {
        SkillDocument {
            name: $name,
            title: $title,
            body: include_str!(concat!("../skill/", $name)),
        }
    };
}

/// Every embedded document, main skill first.
pub const DOCUMENTS: &[SkillDocument] = &[
    document!("SKILL.md", "How to drive xtalforge"),
    document!("mcp.md", "MCP client setup"),
    document!("prompts/build-mof.md", "Prompt: build a MOF"),
    document!("prompts/porosity-report.md", "Prompt: porosity report"),
    document!(
        "prompts/prepare-deposited-cif.md",
        "Prompt: prepare a deposited CIF"
    ),
    document!("references/api.md", "Reference: API and surfaces"),
    document!("references/calculations.md", "Reference: calculations"),
    document!("references/carbon.md", "Reference: carbon builders"),
    document!("references/diagnostics.md", "Reference: diagnostics"),
    document!("references/mof.md", "Reference: MOF builder"),
    document!("references/polymer.md", "Reference: polymers"),
    document!("references/porosity.md", "Reference: porosity"),
    document!("references/prepare.md", "Reference: preparation"),
    document!(
        "references/workspace.md",
        "Reference: documents and history"
    ),
    document!("references/workflows.md", "Reference: workflows"),
];

/// Looks up a document by path or by base name without extension.
#[must_use]
pub fn find(name: &str) -> Option<&'static SkillDocument> {
    let wanted = name.trim().trim_start_matches('/');
    DOCUMENTS.iter().find(|document| {
        document.name == wanted
            || document
                .name
                .strip_suffix(".md")
                .is_some_and(|stem| stem == wanted)
    })
}

/// Default install directory: `$XDG_DATA_HOME/xtalforge/skill`, or
/// `~/.local/share/xtalforge/skill` when the variable is unset.
#[must_use]
pub fn default_directory() -> PathBuf {
    if let Some(base) = std::env::var_os("XDG_DATA_HOME")
        && !base.is_empty()
    {
        return PathBuf::from(base).join("xtalforge").join("skill");
    }
    if let Some(home) = std::env::var_os("HOME")
        && !home.is_empty()
    {
        return PathBuf::from(home)
            .join(".local")
            .join("share")
            .join("xtalforge")
            .join("skill");
    }
    PathBuf::from("xtalforge-skill")
}

/// Writes every document under `directory`, creating parents.
///
/// # Errors
///
/// Returns [`Error::Io`] when a directory cannot be created or a file cannot
/// be written.
pub fn install(directory: &Path) -> Result<Vec<PathBuf>, Error> {
    let mut written = Vec::with_capacity(DOCUMENTS.len());
    for document in DOCUMENTS {
        let path = directory.join(document.name);
        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent)
                .map_err(|error| Error::Io(format!("{}: {error}", parent.display())))?;
        }
        std::fs::write(&path, document.body)
            .map_err(|error| Error::Io(format!("{}: {error}", path.display())))?;
        written.push(path);
    }
    Ok(written)
}
