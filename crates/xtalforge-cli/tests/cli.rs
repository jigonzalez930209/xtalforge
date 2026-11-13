//! Integration tests for the `xtalforge` CLI hello-world surface.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::process::{Command, Output};

fn run(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_xtalforge"))
        .args(args)
        .output()
        .expect("the CLI binary must run")
}

#[test]
fn version_human_output() {
    let output = run(&["version"]);
    assert!(output.status.success());
    let text = String::from_utf8(output.stdout).expect("stdout is UTF-8");
    assert_eq!(
        text.trim(),
        format!("xtalforge {}", env!("CARGO_PKG_VERSION"))
    );
}

#[test]
fn version_json_output() {
    let output = run(&["--json", "version"]);
    assert!(output.status.success());
    let value: serde_json::Value =
        serde_json::from_slice(&output.stdout).expect("stdout must be JSON");
    assert_eq!(
        value["version"],
        serde_json::Value::String(env!("CARGO_PKG_VERSION").to_owned())
    );
}

#[test]
fn cell_volume_human_output() {
    let output = run(&["cell-volume", "1", "1", "1", "90", "90", "90"]);
    assert!(output.status.success());
    let volume: f64 = String::from_utf8(output.stdout)
        .expect("stdout is UTF-8")
        .trim()
        .parse()
        .expect("the volume must be a number");
    assert!((volume - 1.0).abs() < 1e-12);
}

#[test]
fn cell_volume_json_output() {
    let output = run(&["--json", "cell-volume", "2", "2", "3", "90", "90", "90"]);
    assert!(output.status.success());
    let value: serde_json::Value =
        serde_json::from_slice(&output.stdout).expect("stdout must be JSON");
    let volume = value["volume"].as_f64().expect("volume must be a number");
    assert!((volume - 12.0).abs() < 1e-12);
}

#[test]
fn invalid_cell_exits_with_code_three() {
    let output = run(&["cell-volume", "0", "1", "1", "90", "90", "90"]);
    assert_eq!(output.status.code(), Some(3));
    let text = String::from_utf8(output.stderr).expect("stderr is UTF-8");
    assert!(text.contains("invalid input"));
}

#[test]
fn invalid_cell_json_error_has_a_code() {
    let output = run(&["--json", "cell-volume", "0", "1", "1", "90", "90", "90"]);
    assert_eq!(output.status.code(), Some(3));
    let value: serde_json::Value =
        serde_json::from_slice(&output.stderr).expect("stderr must be JSON");
    assert_eq!(value["error"]["code"], "invalid_input");
}

#[test]
fn missing_arguments_use_the_usage_exit_code() {
    let output = run(&["cell-volume", "1"]);
    assert_eq!(output.status.code(), Some(2));
}

#[test]
fn shared_cell_volume_vectors() {
    let vectors: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/vectors/hello_world.json"))
            .expect("shared hello-world vectors must be valid JSON");
    let cases = vectors["cell_volume"]
        .as_array()
        .expect("cell_volume must be an array");
    for case in cases {
        let number = |key: &str| {
            case[key]
                .as_f64()
                .unwrap_or_else(|| panic!("missing numeric field {key}"))
        };
        let output = Command::new(env!("CARGO_BIN_EXE_xtalforge"))
            .arg("--json")
            .arg("cell-volume")
            .arg(number("a").to_string())
            .arg(number("b").to_string())
            .arg(number("c").to_string())
            .arg(number("alpha").to_string())
            .arg(number("beta").to_string())
            .arg(number("gamma").to_string())
            .output()
            .expect("the CLI binary must run");
        assert!(output.status.success());
        let value: serde_json::Value =
            serde_json::from_slice(&output.stdout).expect("stdout must be JSON");
        let got = value["volume"].as_f64().expect("volume must be a number");
        assert!((got - number("expected")).abs() <= number("tolerance"));
    }
}

#[test]
fn element_human_output() {
    let output = run(&["element", "Fe"]);
    assert!(output.status.success());
    let text = String::from_utf8(output.stdout).expect("stdout is UTF-8");
    assert!(text.contains("Fe"));
    assert!(text.contains("Iron"));
    assert!(text.contains("Z=26"));
}

#[test]
fn element_json_output() {
    let output = run(&["--json", "element", "26"]);
    assert!(output.status.success());
    let value: serde_json::Value =
        serde_json::from_slice(&output.stdout).expect("stdout must be JSON");
    assert_eq!(value["symbol"], "Fe");
    assert_eq!(value["z"], 26);
    assert_eq!(value["color"], serde_json::json!([224, 102, 51]));
}

#[test]
fn invalid_element_exits_with_code_three() {
    let output = run(&["element", "Zz"]);
    assert_eq!(output.status.code(), Some(3));
    let text = String::from_utf8(output.stderr).expect("stderr is UTF-8");
    assert!(text.contains("unknown element symbol"));
}

#[test]
fn shared_element_vectors() {
    let vectors: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/vectors/elements.json"))
            .expect("element vectors must be valid JSON");
    for case in vectors["elements"].as_array().expect("elements array") {
        let input = case["input"].as_str().expect("input");
        let output = Command::new(env!("CARGO_BIN_EXE_xtalforge"))
            .args(["--json", "element", input])
            .output()
            .expect("the CLI binary must run");
        assert!(output.status.success());
        let value: serde_json::Value =
            serde_json::from_slice(&output.stdout).expect("stdout must be JSON");
        assert_eq!(value["symbol"], case["symbol"]);
    }
}

#[test]
fn cell_human_output() {
    let output = run(&["cell", "--params", "5", "5", "5", "90", "90", "90"]);
    assert!(output.status.success());
    let text = String::from_utf8(output.stdout).expect("stdout is UTF-8");
    assert!(text.contains("volume=125.000000"));
    assert!(text.contains("a=5.000000"));
}

#[test]
fn cell_json_output() {
    let output = run(&[
        "--json", "cell", "--params", "3", "4", "5", "60", "70", "80",
    ]);
    assert!(output.status.success());
    let value: serde_json::Value =
        serde_json::from_slice(&output.stdout).expect("stdout must be JSON");
    let volume = value["volume"].as_f64().expect("volume");
    assert!((volume - 48.827_604_320_084_98).abs() < 1e-9);
    assert_eq!(value["parameters"][0], 3.0);
}

#[test]
fn invalid_cell_parameters_exit_with_code_three() {
    let output = run(&["cell", "--params", "0", "1", "1", "90", "90", "90"]);
    assert_eq!(output.status.code(), Some(3));
}

#[test]
fn shared_lattice_vectors() {
    let vectors: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/vectors/lattice.json"))
            .expect("lattice vectors must be valid JSON");
    for case in vectors["cases"].as_array().expect("cases array") {
        let number = |name: &str| case[name].as_f64().expect("numeric field").to_string();
        let output = Command::new(env!("CARGO_BIN_EXE_xtalforge"))
            .args([
                "--json",
                "cell",
                "--params",
                &number("a"),
                &number("b"),
                &number("c"),
                &number("alpha"),
                &number("beta"),
                &number("gamma"),
            ])
            .output()
            .expect("the CLI binary must run");
        assert!(output.status.success());
        let value: serde_json::Value =
            serde_json::from_slice(&output.stdout).expect("stdout must be JSON");
        let volume = value["volume"].as_f64().expect("volume");
        assert!((volume - case["volume"].as_f64().expect("expected volume")).abs() < 1e-9);
    }
}

#[test]
fn info_human_and_json() {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/vectors/structure.json");
    let output = run(&["info", &path.to_string_lossy()]);
    assert!(output.status.success());
    let text = String::from_utf8(output.stdout).expect("stdout is UTF-8");
    assert!(text.contains("NaCl"));
    assert!(text.contains("Na Cl"));

    let output = run(&["--json", "info", &path.to_string_lossy()]);
    assert!(output.status.success());
    let value: serde_json::Value =
        serde_json::from_slice(&output.stdout).expect("stdout must be JSON");
    assert_eq!(value["name"], "NaCl");
    assert_eq!(value["formula"], "Na Cl");
    assert!((value["density"].as_f64().unwrap() - 0.540_931_879_932_264_1).abs() < 1e-9);
}

#[test]
fn info_reads_a_cif_file() {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/corpus/simple/Si.cif");
    let output = run(&["--json", "info", &path.to_string_lossy()]);
    assert!(output.status.success());
    let value: serde_json::Value =
        serde_json::from_slice(&output.stdout).expect("stdout must be JSON");
    assert_eq!(value["formula"], "Si8");
    assert_eq!(value["n_sites"], 1);
}

#[test]
fn info_invalid_input_exits_with_code_three() {
    let output = run(&["info", "/definitely/not/a/file"]);
    assert_eq!(output.status.code(), Some(3));
}

#[test]
fn help_lists_the_hello_world_commands_and_flags() {
    let output = run(&["--help"]);
    assert!(output.status.success());
    let text = String::from_utf8(output.stdout).expect("stdout is UTF-8");
    assert!(text.contains("cell-volume"));
    assert!(text.contains("--json"));
    assert!(text.contains("--threads"));
    assert!(text.contains("--limits"));
}

#[test]
fn a_closed_output_pipe_is_not_a_panic() {
    use std::io::Read;
    use std::process::{Command, Stdio};
    let sample = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../docs/public/samples/Si.cif"
    );
    let mut child = Command::new(env!("CARGO_BIN_EXE_xtalforge"))
        .args(["--json", "pxrd", "--step", "0.001", sample])
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("spawn");
    // Read a few bytes, then close the pipe like `| head -c 16`.
    let mut head = [0u8; 16];
    child
        .stdout
        .take()
        .expect("stdout")
        .read_exact(&mut head)
        .expect("read");
    let output = child.wait_with_output().expect("wait");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(!stderr.contains("panicked"), "{stderr}");
}

#[test]
fn fetch_help_lists_the_remote_providers() {
    let output = run(&["fetch", "--help"]);
    assert!(output.status.success());
    let text = String::from_utf8(output.stdout).expect("stdout is UTF-8");
    for provider in ["optimade", "cod", "mp"] {
        assert!(text.contains(provider), "missing {provider} in fetch help");
    }
}
