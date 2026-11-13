//! Built-in modules available without a plugin (Task 14.5.1).

use serde_json::{Value, json};

use crate::registry::{Module, ModuleInfo, ModuleRegistry};
use xtalforge_core::Error;

/// Reports the xtalforge version and the module registry protocol.
#[derive(Debug, Clone, Copy, Default)]
pub struct VersionModule;

impl Module for VersionModule {
    fn info(&self) -> ModuleInfo {
        ModuleInfo {
            name: "version".to_owned(),
            description: "Reports the xtalforge version and module protocol.".to_owned(),
            version: env!("CARGO_PKG_VERSION").to_owned(),
            parameters: Vec::new(),
        }
    }

    fn run(&self, _parameters: &Value) -> Result<Value, Error> {
        Ok(json!({
            "version": xtalforge_core::version(),
            "module_protocol": 1,
        }))
    }
}

/// Registry with the built-in modules.
///
/// # Errors
///
/// Returns an error only if a built-in name collides, which is a programming
/// error and not reachable with the current list.
pub fn registry() -> Result<ModuleRegistry, Error> {
    let mut registry = ModuleRegistry::new();
    registry.register(Box::new(VersionModule))?;
    Ok(registry)
}
