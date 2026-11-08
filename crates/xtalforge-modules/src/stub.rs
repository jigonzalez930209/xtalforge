//! Test stub module (Task 5.10.5).
//!
//! A module that does nothing, so the registry, the parameter forms of front
//! ends, cancellation and run folders can be exercised without a real
//! calculation.

use serde_json::json;
use xtalforge_core::Error;

use crate::registry::{Module, ModuleInfo};

/// The stub module implementation.
#[derive(Debug, Clone, Copy, Default)]
pub struct TestModule;

impl Module for TestModule {
    fn info(&self) -> ModuleInfo {
        ModuleInfo {
            name: "test".to_owned(),
            description: "Stub module that reports availability and does nothing.".to_owned(),
            version: env!("CARGO_PKG_VERSION").to_owned(),
            parameters: Vec::new(),
        }
    }

    fn available(&self) -> bool {
        true
    }

    fn run(&self, _parameters: &serde_json::Value) -> Result<serde_json::Value, Error> {
        Ok(json!({
            "ok": true,
            "module": "test",
            "message": "test module ran",
        }))
    }
}

/// Creates the stub module for registration and surface tests.
#[must_use]
pub fn test_module() -> Box<dyn Module> {
    Box::new(TestModule)
}
