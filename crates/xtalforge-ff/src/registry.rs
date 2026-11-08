//! Built-in calculator registration (Task 4.1.2).
//!
//! Ported from Crystal-Builder `xtal/ff/registry.py`. The registry is the one
//! list the CLI `--engine` flag, the Python and WASM bindings and the agent
//! tools read. UFF is compiled in, so it is always available; engines that
//! need a binary, a feature or model weights will register an
//! [`crate::Availability::Unavailable`] reason from their probe.

use std::sync::OnceLock;

use crate::Error;
use crate::calculator::{
    Availability, Calculator, CalculatorCapabilities, CalculatorInfo, CalculatorOptions,
    CalculatorRegistry,
};
use crate::uff::calculator::UffCalculator;

static UFF_INFO: OnceLock<CalculatorInfo> = OnceLock::new();

fn uff_info() -> &'static CalculatorInfo {
    UFF_INFO.get_or_init(|| CalculatorInfo {
        name: "uff".to_owned(),
        label: "UFF".to_owned(),
        description: "Universal Force Field (Rappe et al. 1992) with the \
                      UFF4MOF and UFF4MOF-II extensions -- covers the whole \
                      periodic table, native and in-process"
            .to_owned(),
        availability: Availability::Available,
        capabilities: CalculatorCapabilities {
            periodic: true,
            forces: true,
            stress: true,
            charges: false,
            in_process: true,
            thread_safe: true,
        },
    })
}

fn build_uff(options: CalculatorOptions) -> Result<Box<dyn Calculator>, Error> {
    Ok(Box::new(UffCalculator::new(options)?))
}

/// The registry with every built-in engine.
#[must_use]
pub fn builtin_registry() -> CalculatorRegistry {
    let mut registry = CalculatorRegistry::new();
    registry.add_builtin(uff_info(), build_uff);
    registry
}
