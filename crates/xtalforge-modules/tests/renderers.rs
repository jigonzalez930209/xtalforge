//! HTML and JSON report renderers plus the module listing (Task 14.5.1).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use serde_json::Value;
use xtalforge_modules::{
    Module, ModuleInfo, ModuleListing, ModuleRegistry, Report, ReportSection, test_module,
};

fn sample() -> Report {
    Report {
        run_id: "run-7".to_owned(),
        title: "PXRD of Si".to_owned(),
        summary: "One phase, three reflections.".to_owned(),
        sections: vec![ReportSection {
            heading: "Peaks".to_owned(),
            body: "The three strongest peaks are 111, 220 and 311.".to_owned(),
        }],
        references: vec!["Crystal-Builder golden data".to_owned()],
    }
}

#[derive(Debug)]
struct Unavailable;

impl Module for Unavailable {
    fn info(&self) -> ModuleInfo {
        ModuleInfo {
            name: "zeta".to_owned(),
            description: "Needs an external program".to_owned(),
            version: "0.9".to_owned(),
            parameters: vec!["cutoff".to_owned()],
        }
    }

    fn available(&self) -> bool {
        false
    }
}

#[test]
fn html_has_the_document_structure() {
    let html = sample().to_html();
    assert!(html.starts_with("<!DOCTYPE html>\n<html lang=\"en\">\n"));
    assert!(html.contains("<meta charset=\"utf-8\">"));
    assert!(html.contains("<title>PXRD of Si</title>"));
    assert!(html.contains("<h1>PXRD of Si</h1>"));
    assert!(html.contains("<p class=\"run-id\">Run: <code>run-7</code></p>"));
    assert!(html.contains("<p class=\"summary\">One phase, three reflections.</p>"));
    assert!(html.contains("<section>\n<h2>Peaks</h2>"));
    assert!(html.contains("<pre>The three strongest peaks are 111, 220 and 311.</pre>"));
    assert!(html.contains("<h2>References</h2>"));
    assert!(html.contains("<li>Crystal-Builder golden data</li>"));
    assert!(html.ends_with("</body>\n</html>\n"));
    assert!(!html.ends_with("\n\n"));
}

#[test]
fn html_omits_empty_parts_like_markdown() {
    let report = Report {
        run_id: "  ".to_owned(),
        title: "Empty".to_owned(),
        summary: String::new(),
        sections: Vec::new(),
        references: Vec::new(),
    };
    let html = report.to_html();
    assert!(html.contains("<h1>Empty</h1>"));
    assert!(!html.contains("run-id"));
    assert!(!html.contains("summary"));
    assert!(!html.contains("<h2>"));
    assert!(!html.contains("<ul>"));
    assert!(report.to_markdown().contains("# Empty"));
}

#[test]
fn html_escapes_every_dynamic_string() {
    let report = Report {
        run_id: "run<&>".to_owned(),
        title: "A <script>alert(\"x\")</script> & 'quote'".to_owned(),
        summary: "1 < 2 & 3 > 2".to_owned(),
        sections: vec![ReportSection {
            heading: "H <b>".to_owned(),
            body: "<img src=x onerror=alert(1)> & `code`".to_owned(),
        }],
        references: vec!["R&D <https://example.org?a=1&b=2>".to_owned()],
    };
    let html = report.to_html();
    assert!(
        html.contains("A &lt;script&gt;alert(&quot;x&quot;)&lt;/script&gt; &amp; &#39;quote&#39;")
    );
    assert!(html.contains("<code>run&lt;&amp;&gt;</code>"));
    assert!(html.contains("<p class=\"summary\">1 &lt; 2 &amp; 3 &gt; 2</p>"));
    assert!(html.contains("<h2>H &lt;b&gt;</h2>"));
    assert!(html.contains("<pre>&lt;img src=x onerror=alert(1)&gt; &amp; `code`</pre>"));
    assert!(html.contains("<li>R&amp;D &lt;https://example.org?a=1&amp;b=2&gt;</li>"));
    assert!(!html.contains("<script>"));
    assert!(!html.contains("<img "));
}

#[test]
fn json_round_trips() {
    let report = sample();
    let text = report.to_json().unwrap();
    assert!(text.starts_with("{\n"));
    assert!(text.ends_with('}'));
    let value: Value = serde_json::from_str(&text).unwrap();
    assert_eq!(value["title"], "PXRD of Si");
    assert_eq!(value["sections"][0]["heading"], "Peaks");
    assert_eq!(value["references"][0], "Crystal-Builder golden data");
    let parsed = Report::from_json(&text).unwrap();
    assert_eq!(parsed, report);
}

#[test]
fn json_escapes_and_keeps_unicode() {
    let report = Report {
        run_id: "run-1".to_owned(),
        title: "A \"quoted\" \\ title — ñ".to_owned(),
        summary: "line one\nline two".to_owned(),
        sections: Vec::new(),
        references: Vec::new(),
    };
    let text = report.to_json().unwrap();
    let parsed = Report::from_json(&text).unwrap();
    assert_eq!(parsed, report);
    assert!(text.contains("\\\"quoted\\\""));
    assert!(text.contains("\\\\"));
    assert!(text.contains("ñ"));
}

#[test]
fn invalid_json_is_a_typed_error() {
    let error = Report::from_json("{\"title\": 1}").unwrap_err();
    assert_eq!(error.code(), "invalid_input");
    assert!(error.message().contains("invalid report"));
}

#[test]
fn the_listing_is_sorted_and_reports_availability() {
    let mut registry = ModuleRegistry::new();
    registry.register(Box::new(Unavailable)).unwrap();
    registry.register(test_module()).unwrap();
    let listing = registry.listing();
    assert_eq!(listing.len(), 2);
    assert_eq!(listing[0].name, "test");
    assert_eq!(listing[1].name, "zeta");
    assert!(listing[0].available);
    assert!(!listing[1].available);
    assert_eq!(listing[1].description, "Needs an external program");
    assert_eq!(listing[1].version, "0.9");
    assert_eq!(listing[1].parameters, vec!["cutoff".to_owned()]);
    assert_eq!(registry.list().len(), 2);
    let text = registry.listing_json().unwrap();
    let parsed: Vec<ModuleListing> = serde_json::from_str(&text).unwrap();
    assert_eq!(parsed, listing);
    assert!(text.starts_with("[\n"));
}

#[test]
fn an_empty_registry_lists_nothing() {
    let registry = ModuleRegistry::new();
    assert!(registry.listing().is_empty());
    let parsed: Vec<ModuleListing> =
        serde_json::from_str(&registry.listing_json().unwrap()).unwrap();
    assert!(parsed.is_empty());
}
