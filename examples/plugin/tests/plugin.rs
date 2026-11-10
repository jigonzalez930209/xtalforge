//! The example plugin, exercised through the public registries (Task 8.3.1).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge::core::ScalarRegistry;
use xtalforge::ff::{CalculatorOptions, CalculatorRegistry};
use xtalforge::io::registry::FormatRegistry;
use xtalforge::modules::ModuleRegistry;
use xtalforge::view::StyleRegistry;
use xtalforge_plugin_example::{
    FORMAT_NAME, example_document, read_two_line, register_all, write_two_line,
};

#[test]
fn every_registry_accepts_the_plugin() {
    let mut formats = FormatRegistry::new();
    let mut calculators = CalculatorRegistry::default();
    let mut modules = ModuleRegistry::new();
    let mut scalars = ScalarRegistry::new();
    let mut styles = StyleRegistry::new();
    register_all(
        &mut formats,
        &mut calculators,
        &mut modules,
        &mut scalars,
        &mut styles,
    )
    .expect("register");

    assert!(formats.find(FORMAT_NAME).is_some());
    assert!(calculators.names().contains(&"zero"));
    assert!(modules.get("site-count").is_some());
    assert!(scalars.names().contains(&"first_fractional"));
    assert!(styles.get("large-atoms").is_some());
}

#[test]
fn the_custom_format_reads_and_writes() {
    let mut formats = FormatRegistry::new();
    formats
        .register(xtalforge::io::registry::FormatSpec {
            name: FORMAT_NAME,
            extensions: &[".2l"],
            read: Some(read_two_line),
            write: Some(write_two_line),
            capabilities: xtalforge::io::registry::Capabilities::default(),
        })
        .expect("format");
    let text = "cell 4.0\nNa 0 0 0\nCl 0.5 0.5 0.5\n";
    let structure = formats.read_with(text, ".2l").expect("read");
    assert_eq!(structure.sites().len(), 2);
    assert_eq!(structure.formula(), "Na Cl");
    let written = formats.write_with(&structure, FORMAT_NAME).expect("write");
    let again = formats.read_with(&written, ".2l").expect("re-read");
    assert_eq!(again.sites().len(), 2);
    assert!(formats.find("nope").is_none());
}

#[test]
fn the_custom_calculator_answers() {
    let mut calculators = CalculatorRegistry::default();
    register_all(
        &mut FormatRegistry::new(),
        &mut calculators,
        &mut ModuleRegistry::new(),
        &mut ScalarRegistry::new(),
        &mut StyleRegistry::new(),
    )
    .expect("register");
    let calculator = calculators
        .build("zero", &CalculatorOptions::default())
        .expect("build");
    let document = example_document().expect("document");
    assert_eq!(calculator.name(), "zero");
    assert!(calculator.energy(document.structure()).expect("energy") == 0.0);
    let (energy, forces) = calculator
        .energy_forces(document.structure())
        .expect("forces");
    assert_eq!(energy, 0.0);
    assert_eq!(forces.len(), 2);
}

#[test]
fn the_custom_module_scalar_and_style_work() {
    let mut modules = ModuleRegistry::new();
    let mut scalars = ScalarRegistry::new();
    let mut styles = StyleRegistry::new();
    register_all(
        &mut FormatRegistry::new(),
        &mut CalculatorRegistry::default(),
        &mut modules,
        &mut scalars,
        &mut styles,
    )
    .expect("register");
    let result = modules
        .get("site-count")
        .expect("module")
        .run(&serde_json::json!({"n": 2}))
        .expect("run");
    assert_eq!(result["plugin"], "site-count");
    let document = example_document().expect("document");
    let values = scalars
        .compute("first_fractional", document.structure())
        .expect("scalar");
    assert_eq!(values, vec![0.0, 0.5]);
    let settings = styles.get("large-atoms").expect("style");
    assert!(settings.atom_scale > 1.0);
}

#[test]
fn duplicate_registration_is_refused() {
    let mut formats = FormatRegistry::new();
    let mut calculators = CalculatorRegistry::default();
    let mut modules = ModuleRegistry::new();
    let mut scalars = ScalarRegistry::new();
    let mut styles = StyleRegistry::new();
    register_all(
        &mut formats,
        &mut calculators,
        &mut modules,
        &mut scalars,
        &mut styles,
    )
    .expect("first");
    let error = register_all(
        &mut formats,
        &mut calculators,
        &mut modules,
        &mut scalars,
        &mut styles,
    )
    .expect_err("second registration");
    assert_eq!(error.code(), "invalid_input");
}
