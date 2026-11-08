//! Run reports (Tasks 5.10.4 and 5.10.6).
//!
//! A report is presentation, not results: every number in it was already
//! written to the run folder. [`Report::to_markdown`] is the one renderer
//! needed by the CLI and the docs; the structured form serializes to JSON for
//! the bindings.

use serde::{Deserialize, Serialize};

/// One headed section of a report.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ReportSection {
    /// Section heading.
    pub heading: String,
    /// Section body, in Markdown.
    pub body: String,
}

/// A human-readable report attached to a run.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Report {
    /// Stable id of the run that produced it.
    pub run_id: String,
    /// Report title.
    pub title: String,
    /// One-paragraph summary.
    pub summary: String,
    /// Body sections in presentation order.
    pub sections: Vec<ReportSection>,
    /// Citations of methods, programs and databases used by the run.
    pub references: Vec<String>,
}

impl Report {
    /// Renders the report as Markdown.
    ///
    /// Empty summaries and reference lists are omitted; sections keep their
    /// order. The output ends with exactly one newline.
    #[must_use]
    pub fn to_markdown(&self) -> String {
        let mut out = String::new();
        out.push_str(&format!("# {}\n\n", self.title));
        if !self.run_id.trim().is_empty() {
            out.push_str(&format!("Run: `{}`\n\n", self.run_id));
        }
        if !self.summary.trim().is_empty() {
            out.push_str(self.summary.trim());
            out.push_str("\n\n");
        }
        for section in &self.sections {
            out.push_str(&format!("## {}\n\n", section.heading));
            let body = section.body.trim();
            if !body.is_empty() {
                out.push_str(body);
                out.push_str("\n\n");
            }
        }
        if !self.references.is_empty() {
            out.push_str("## References\n\n");
            for reference in &self.references {
                out.push_str(&format!("- {reference}\n"));
            }
        }
        let mut trimmed = out.trim_end().to_owned();
        trimmed.push('\n');
        trimmed
    }
}
