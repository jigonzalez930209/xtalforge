//! Generic module registry (Task 3.13).

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use xtalforge_core::Error;

/// Description of one module.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ModuleInfo {
    /// Stable name.
    pub name: String,
    /// One-line description.
    pub description: String,
    /// Version of the module implementation.
    pub version: String,
    /// Parameter names it accepts.
    pub parameters: Vec<String>,
}

/// A module that can be registered.
pub trait Module: std::fmt::Debug + Send + Sync {
    /// Description.
    fn info(&self) -> ModuleInfo;
}

/// Registry of modules by name.
#[derive(Debug, Default)]
pub struct ModuleRegistry {
    modules: BTreeMap<String, Box<dyn Module>>,
}

impl ModuleRegistry {
    /// Empty registry.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Registers a module.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a duplicate or empty name.
    pub fn register(&mut self, module: Box<dyn Module>) -> Result<(), Error> {
        let info = module.info();
        if info.name.trim().is_empty() {
            return Err(Error::InvalidInput(
                "a module needs a non-empty name".to_owned(),
            ));
        }
        if self.modules.contains_key(&info.name) {
            return Err(Error::InvalidInput(format!(
                "module {} is already registered",
                info.name
            )));
        }
        self.modules.insert(info.name, module);
        Ok(())
    }

    /// Unregisters a module.
    pub fn unregister(&mut self, name: &str) -> bool {
        self.modules.remove(name).is_some()
    }

    /// Gets a module.
    #[must_use]
    pub fn get(&self, name: &str) -> Option<&dyn Module> {
        self.modules.get(name).map(Box::as_ref)
    }

    /// Lists all module descriptions.
    #[must_use]
    pub fn list(&self) -> Vec<ModuleInfo> {
        self.modules.values().map(|module| module.info()).collect()
    }
}
