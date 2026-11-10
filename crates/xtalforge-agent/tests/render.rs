//! Pictures for the agent (Task 8.1.6).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_agent::render::{MAX_DIMENSION, RenderOptions, render_png, render_svg};
use xtalforge_io::registry::read_with;

const NACL: &str = include_str!("../../../tests/corpus/simple/NaCl.cif");

fn structure() -> xtalforge_core::Structure {
    read_with(NACL, "cif").expect("NaCl")
}

#[test]
fn svg_has_the_requested_size_and_draws_atoms() {
    let options = RenderOptions {
        width: 320,
        height: 240,
        ..Default::default()
    };
    let svg = render_svg(&structure(), &options).expect("svg");
    assert!(svg.starts_with("<svg"));
    assert!(svg.contains("width=\"320\""));
    assert!(svg.contains("height=\"240\""));
    assert!(svg.contains("<circle"));
    assert!(svg.trim_end().ends_with("</svg>"));
}

#[test]
fn png_has_a_valid_signature_and_dimensions() {
    let options = RenderOptions {
        width: 160,
        height: 120,
        ..Default::default()
    };
    let png = render_png(&structure(), &options).expect("png");
    assert_eq!(&png[..8], &[137, 80, 78, 71, 13, 10, 26, 10]);
    assert_eq!(&png[12..16], b"IHDR");
    let width = u32::from_be_bytes([png[16], png[17], png[18], png[19]]);
    let height = u32::from_be_bytes([png[20], png[21], png[22], png[23]]);
    assert_eq!(width, 160);
    assert_eq!(height, 120);
}

#[test]
fn an_empty_structure_is_refused() {
    let hydrogen = xtalforge_core::Structure::molecule(
        "H2",
        &[
            ("H".to_owned(), [0.0, 0.0, 0.0]),
            ("H".to_owned(), [0.74, 0.0, 0.0]),
        ],
        8.0,
    )
    .expect("H2");
    let mut options = RenderOptions::default();
    options.settings.hide_hydrogens = true;
    options.settings.show_cell = false;
    let error = render_svg(&hydrogen, &options).unwrap_err();
    assert_eq!(error.code(), "invalid_input");
}

#[test]
fn dimensions_are_checked() {
    let options = RenderOptions {
        width: MAX_DIMENSION + 1,
        ..Default::default()
    };
    let error = render_svg(&structure(), &options).unwrap_err();
    assert_eq!(error.code(), "invalid_argument");
}

#[test]
fn the_document_verb_writes_a_file() {
    use serde_json::json;
    use xtalforge_agent::Session;
    let mut session = Session::new();
    session
        .run(
            "load",
            &json!({"text": NACL, "format": "cif", "name": "nacl"}),
        )
        .expect("load");
    let mut path = std::env::temp_dir();
    path.push(format!("xtalforge-agent-render-{}.png", std::process::id()));
    let answer = session
        .run(
            "render",
            &json!({"format": "png", "path": path.to_string_lossy(), "width": 200, "height": 150}),
        )
        .expect("render");
    let bytes = std::fs::read(&path).expect("written image");
    assert_eq!(&bytes[..8], &[137, 80, 78, 71, 13, 10, 26, 10]);
    assert!(matches!(
        answer,
        xtalforge_agent::Answer::Picture { ref format, .. } if format == "png"
    ));
    let _ = std::fs::remove_file(&path);
}
