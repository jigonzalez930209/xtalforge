//! Input pattern expansion for the runner (Task 16.1.2).
//!
//! `--inputs` accepts a literal file, a directory or a shell-style pattern
//! with `*` (within one path component), `?` (one character) and `**` (any
//! number of components). Expansion is deterministic: results are sorted and
//! deduplicated, so the report order never depends on the filesystem.
//!
//! The matcher does not follow symbolic links and refuses to descend more
//! than 64 levels.

use std::collections::BTreeSet;
use std::path::{Component, Path, PathBuf};

use super::WorkflowError;

const MAX_DEPTH: usize = 64;

/// Expands every pattern into a sorted, deduplicated list of files.
///
/// A literal path must exist; a pattern that matches nothing contributes
/// nothing, so a caller can tell "no inputs" apart from a typo.
///
/// # Errors
///
/// Returns [`WorkflowError::Inputs`] for a missing literal path or a pattern
/// whose base directory does not exist, and [`WorkflowError::Io`] when a
/// directory cannot be read.
pub fn expand_inputs<I, S>(patterns: I) -> Result<Vec<PathBuf>, WorkflowError>
where
    I: IntoIterator<Item = S>,
    S: AsRef<str>,
{
    let mut found = BTreeSet::new();
    for pattern in patterns {
        let text = pattern.as_ref().trim();
        if text.is_empty() {
            continue;
        }
        expand_one(text, &mut found)?;
    }
    Ok(found.into_iter().collect())
}

fn expand_one(text: &str, found: &mut BTreeSet<PathBuf>) -> Result<(), WorkflowError> {
    let normalized = text.replace('\\', "/");
    if !has_wildcard(&normalized) {
        let path = PathBuf::from(&normalized);
        if path.is_file() {
            found.insert(path);
            return Ok(());
        }
        if path.is_dir() {
            for entry in read_dir(&path)? {
                if entry.file_type().is_ok_and(|kind| kind.is_file()) {
                    found.insert(entry.path());
                }
            }
            return Ok(());
        }
        return Err(WorkflowError::Inputs {
            pattern: text.to_owned(),
            reason: "does not exist".to_owned(),
        });
    }
    let absolute = normalized.starts_with('/');
    let parts: Vec<&str> = normalized
        .split('/')
        .filter(|part| !part.is_empty() && *part != ".")
        .collect();
    let first_wildcard = parts.iter().position(|part| has_wildcard(part));
    let Some(first_wildcard) = first_wildcard else {
        return Ok(());
    };
    let (base_parts, match_parts) = parts.split_at(first_wildcard);
    let base = base_path(absolute, base_parts);
    if !base.is_dir() {
        return Err(WorkflowError::Inputs {
            pattern: text.to_owned(),
            reason: format!("base directory {} does not exist", base.display()),
        });
    }
    let mut candidates = Vec::new();
    walk(&base, 0, &mut candidates)?;
    for candidate in candidates {
        let relative = candidate.strip_prefix(&base).unwrap_or(&candidate);
        let components: Vec<String> = relative
            .components()
            .filter_map(|component| match component {
                Component::Normal(value) => Some(value.to_string_lossy().into_owned()),
                _ => None,
            })
            .collect();
        let references: Vec<&str> = components.iter().map(String::as_str).collect();
        if matches(match_parts, &references) {
            found.insert(candidate);
        }
    }
    Ok(())
}

fn base_path(absolute: bool, parts: &[&str]) -> PathBuf {
    let joined = parts.join("/");
    if absolute {
        if joined.is_empty() {
            PathBuf::from("/")
        } else {
            PathBuf::from("/").join(joined)
        }
    } else if joined.is_empty() {
        PathBuf::from(".")
    } else {
        PathBuf::from(joined)
    }
}

fn read_dir(path: &Path) -> Result<Vec<std::fs::DirEntry>, WorkflowError> {
    let entries = std::fs::read_dir(path).map_err(|error| WorkflowError::Io {
        path: path.display().to_string(),
        reason: error.to_string(),
    })?;
    let mut collected = Vec::new();
    for entry in entries {
        collected.push(entry.map_err(|error| WorkflowError::Io {
            path: path.display().to_string(),
            reason: error.to_string(),
        })?);
    }
    Ok(collected)
}

fn walk(directory: &Path, depth: usize, files: &mut Vec<PathBuf>) -> Result<(), WorkflowError> {
    if depth > MAX_DEPTH {
        return Ok(());
    }
    for entry in read_dir(directory)? {
        let Ok(kind) = entry.file_type() else {
            continue;
        };
        if kind.is_file() {
            files.push(entry.path());
        } else if kind.is_dir() {
            walk(&entry.path(), depth + 1, files)?;
        }
    }
    Ok(())
}

fn has_wildcard(text: &str) -> bool {
    text.contains('*') || text.contains('?')
}

fn matches(pattern: &[&str], path: &[&str]) -> bool {
    let mut pattern = pattern;
    let mut path = path;
    loop {
        match (pattern.split_first(), path.split_first()) {
            (None, None) => return true,
            (Some((&"**", rest)), _) => {
                if matches(rest, path) {
                    return true;
                }
                let Some((_, tail)) = path.split_first() else {
                    return false;
                };
                path = tail;
            }
            (Some((part, rest)), Some((name, tail))) => {
                if !match_component(part, name) {
                    return false;
                }
                pattern = rest;
                path = tail;
            }
            _ => return false,
        }
    }
}

fn match_component(pattern: &str, text: &str) -> bool {
    let pattern: Vec<char> = pattern.chars().collect();
    let text: Vec<char> = text.chars().collect();
    let (mut pattern_index, mut text_index) = (0, 0);
    let mut star: Option<usize> = None;
    let mut mark = 0;
    while text_index < text.len() {
        match (pattern.get(pattern_index), text.get(text_index)) {
            (Some('*'), _) => {
                star = Some(pattern_index);
                pattern_index += 1;
                mark = text_index;
            }
            (Some('?'), _) => {
                pattern_index += 1;
                text_index += 1;
            }
            (Some(expected), Some(found)) if expected == found => {
                pattern_index += 1;
                text_index += 1;
            }
            _ => {
                let Some(star_index) = star else {
                    return false;
                };
                pattern_index = star_index + 1;
                mark += 1;
                text_index = mark;
            }
        }
    }
    while pattern.get(pattern_index) == Some(&'*') {
        pattern_index += 1;
    }
    pattern_index == pattern.len()
}
