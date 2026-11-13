//! Run reports (Tasks 5.10.4, 5.10.6 and 14.5.1).
//!
//! A report is presentation, not results: every number in it was already
//! written to the run folder. [`Report::to_markdown`] is the renderer used by
//! the CLI and the docs, [`Report::to_html`] produces a standalone page, and
//! [`Report::to_json`] with [`Report::from_json`] round-trips the structured
//! form used by the bindings.

use serde::{Deserialize, Serialize};
use xtalforge_core::Error;

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

    /// Renders the report as a standalone HTML 5 document.
    ///
    /// Every dynamic string is HTML-escaped. Section bodies are Markdown in
    /// the data model but are not converted: they appear verbatim inside a
    /// `<pre>` block, so the renderer never interprets untrusted HTML. Empty
    /// run ids, summaries and reference lists are omitted, matching
    /// [`Report::to_markdown`]. The output ends with exactly one newline.
    #[must_use]
    pub fn to_html(&self) -> String {
        let mut out = String::new();
        out.push_str("<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n");
        out.push_str("<meta charset=\"utf-8\">\n");
        out.push_str("<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n");
        out.push_str(&format!("<title>{}</title>\n", escape_html(&self.title)));
        out.push_str("<style>\n");
        out.push_str(
            "body{font-family:system-ui,sans-serif;max-width:48rem;margin:2rem auto;padding:0 1rem;line-height:1.5;color:#1f2328}\n",
        );
        out.push_str("h1,h2{line-height:1.25}\n");
        out.push_str(
            "pre{background:#f6f8fa;padding:.75rem;border-radius:.375rem;overflow-x:auto;white-space:pre-wrap}\n",
        );
        out.push_str("code{background:#f6f8fa;padding:.1rem .25rem;border-radius:.25rem}\n");
        out.push_str("</style>\n</head>\n<body>\n");
        out.push_str(&format!("<h1>{}</h1>\n", escape_html(&self.title)));
        if !self.run_id.trim().is_empty() {
            out.push_str(&format!(
                "<p class=\"run-id\">Run: <code>{}</code></p>\n",
                escape_html(&self.run_id)
            ));
        }
        if !self.summary.trim().is_empty() {
            out.push_str(&format!(
                "<p class=\"summary\">{}</p>\n",
                escape_html(self.summary.trim())
            ));
        }
        for section in &self.sections {
            out.push_str("<section>\n");
            out.push_str(&format!("<h2>{}</h2>\n", escape_html(&section.heading)));
            let body = section.body.trim();
            if !body.is_empty() {
                out.push_str(&format!("<pre>{}</pre>\n", escape_html(body)));
            }
            out.push_str("</section>\n");
        }
        if !self.references.is_empty() {
            out.push_str("<h2>References</h2>\n<ul>\n");
            for reference in &self.references {
                out.push_str(&format!("<li>{}</li>\n", escape_html(reference)));
            }
            out.push_str("</ul>\n");
        }
        out.push_str("</body>\n</html>\n");
        out
    }

    /// Serializes the report as pretty JSON.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when serialization fails.
    pub fn to_json(&self) -> Result<String, Error> {
        serde_json::to_string_pretty(self).map_err(|error| {
            Error::InvalidInput(format!("could not serialize the report: {error}"))
        })
    }

    /// Parses a report from JSON.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for invalid JSON or a missing field.
    pub fn from_json(text: &str) -> Result<Self, Error> {
        serde_json::from_str(text)
            .map_err(|error| Error::InvalidInput(format!("invalid report: {error}")))
    }
}

/// Escapes `&`, `<`, `>`, `"` and `'` for HTML text and attribute contexts.
fn escape_html(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    for character in text.chars() {
        match character {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&#39;"),
            _ => out.push(character),
        }
    }
    out
}
