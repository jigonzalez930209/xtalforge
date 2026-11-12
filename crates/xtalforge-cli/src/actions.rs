//! Nested CLI actions kept out of the main command tree.

use clap::Subcommand;

#[derive(Debug, Subcommand)]
pub(crate) enum ProjectAction {
    /// Pack a structure into a `.xtalproj` archive.
    Pack {
        /// Input structure path.
        path: String,
        /// Output archive path.
        out: String,
    },
    /// Print the parts of a project archive.
    Info {
        /// Archive path.
        path: String,
    },
    /// Extract a project into a directory.
    Extract {
        /// Archive path.
        path: String,
        /// Destination directory.
        dir: String,
    },
}

/// Model cache actions.
#[derive(Debug, Subcommand)]
pub(crate) enum ModelAction {
    /// List the catalogue and the cache state.
    List,
    /// Download a model into the cache and record its checksum.
    Fetch {
        /// Catalogue name: mace-mp-0, orb-v3 or mattersim.
        name: String,
        /// Artifact URL chosen by the user.
        #[arg(long)]
        url: String,
    },
    /// Recompute the checksum of a cached model.
    Verify {
        /// Catalogue name.
        name: String,
    },
    /// Remove a cached model and its metadata.
    Remove {
        /// Catalogue name.
        name: String,
    },
}

/// Skill document actions.
#[derive(Debug, Subcommand)]
pub(crate) enum SkillAction {
    /// Print the default install directory.
    Path,
    /// List the embedded documents.
    List,
    /// Print one document.
    Show {
        /// Path or base name, for example `SKILL` or `references/api.md`.
        name: String,
    },
    /// Write every document into a directory.
    Install {
        /// Destination; defaults to the user data directory.
        dir: Option<String>,
    },
}
