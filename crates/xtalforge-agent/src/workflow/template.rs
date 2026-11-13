//! Path templates for `read` and `write` steps (Task 16.1.1).
//!
//! A template is plain text with `{{name}}` placeholders. Validation happens
//! when the workflow is parsed, so a run never fails because of a typo in a
//! path.

use std::path::Path;

use super::WorkflowError;

/// Placeholders accepted in a `read` or `write` path.
pub(super) const PLACEHOLDERS: &[&str] = &["input", "stem", "ext", "dir", "index"];

/// Values available to the placeholders of a `read` or `write` path.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TemplateContext {
    /// Input path as the caller wrote it, for example `in/NaCl.cif`.
    pub input: String,
    /// File stem, for example `NaCl`.
    pub stem: String,
    /// Extension without the dot, for example `cif`.
    pub ext: String,
    /// Parent directory, for example `in`.
    pub dir: String,
    /// 1-based position in the deterministic input order.
    pub index: usize,
}

impl TemplateContext {
    /// Builds the context for one input path.
    #[must_use]
    pub fn new(path: &Path, index: usize) -> Self {
        let input = path.to_string_lossy().into_owned();
        let stem = path
            .file_stem()
            .and_then(|stem| stem.to_str())
            .unwrap_or_default()
            .to_owned();
        let ext = path
            .extension()
            .and_then(|extension| extension.to_str())
            .unwrap_or_default()
            .to_owned();
        let dir = path
            .parent()
            .filter(|parent| !parent.as_os_str().is_empty())
            .map(|parent| parent.to_string_lossy().into_owned())
            .unwrap_or_default();
        Self {
            input,
            stem,
            ext,
            dir,
            index,
        }
    }

    fn value(&self, name: &str) -> Option<String> {
        match name {
            "input" => Some(self.input.clone()),
            "stem" => Some(self.stem.clone()),
            "ext" => Some(self.ext.clone()),
            "dir" => Some(self.dir.clone()),
            "index" => Some(self.index.to_string()),
            _ => None,
        }
    }
}

/// Checks a path template and every placeholder it uses.
pub(super) fn validate_template(template: &str, path: &str) -> Result<(), WorkflowError> {
    if template.trim().is_empty() {
        return Err(invalid(path, "must not be empty"));
    }
    let mut rest = template;
    while let Some(open) = rest.find("{{") {
        let (_head, tail) = rest.split_at(open);
        let Some(after) = tail.strip_prefix("{{") else {
            return Err(invalid(path, "unclosed `{{` placeholder"));
        };
        let Some(close) = after.find("}}") else {
            return Err(invalid(path, "unclosed `{{` placeholder"));
        };
        let (name, tail) = after.split_at(close);
        let Some(next) = tail.strip_prefix("}}") else {
            return Err(invalid(path, "unclosed `}}` placeholder"));
        };
        if !PLACEHOLDERS.contains(&name.trim()) {
            return Err(invalid(
                path,
                format!("unknown placeholder {{{{{}}}}}", name.trim()),
            ));
        }
        rest = next;
    }
    Ok(())
}

/// Expands a validated path template.
///
/// # Errors
///
/// Returns [`WorkflowError::Invalid`] for an unknown or unclosed placeholder;
/// [`super::Workflow::validate`] rejects those before a run starts.
pub fn render_template(template: &str, context: &TemplateContext) -> Result<String, WorkflowError> {
    let mut output = String::with_capacity(template.len());
    let mut rest = template;
    loop {
        let Some(open) = rest.find("{{") else {
            output.push_str(rest);
            return Ok(output);
        };
        let (head, tail) = rest.split_at(open);
        output.push_str(head);
        let Some(after) = tail.strip_prefix("{{") else {
            return Err(invalid("path", "unclosed `{{` placeholder"));
        };
        let Some(close) = after.find("}}") else {
            return Err(invalid("path", "unclosed `{{` placeholder"));
        };
        let (name, tail) = after.split_at(close);
        let Some(next) = tail.strip_prefix("}}") else {
            return Err(invalid("path", "unclosed `}}` placeholder"));
        };
        let name = name.trim();
        let value = context
            .value(name)
            .ok_or_else(|| invalid("path", format!("unknown placeholder {{{{{name}}}}}")))?;
        output.push_str(&value);
        rest = next;
    }
}

fn invalid(path: impl Into<String>, reason: impl Into<String>) -> WorkflowError {
    WorkflowError::Invalid {
        path: path.into(),
        reason: reason.into(),
    }
}
