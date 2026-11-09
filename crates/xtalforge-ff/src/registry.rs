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
use crate::quantum::{QuantumMethod, adapter, capabilities};
use crate::uff::calculator::UffCalculator;

static UFF_INFO: OnceLock<CalculatorInfo> = OnceLock::new();
static GFN0_INFO: OnceLock<CalculatorInfo> = OnceLock::new();
static GFN1_INFO: OnceLock<CalculatorInfo> = OnceLock::new();
static GFN2_INFO: OnceLock<CalculatorInfo> = OnceLock::new();

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

fn quantum_info(
    method: QuantumMethod,
    slot: &'static OnceLock<CalculatorInfo>,
) -> &'static CalculatorInfo {
    slot.get_or_init(|| {
        let capabilities = capabilities(method);
        CalculatorInfo {
            name: method.name().to_owned(),
            label: format!("sci-form {}", method.name().to_uppercase()),
            description: format!(
                "In-process {} semi-empirical method from sci-form; molecular only",
                method.name().to_uppercase()
            ),
            availability: if cfg!(feature = "sci-form") {
                Availability::Available
            } else {
                Availability::Unavailable(
                    "sci-form support is not compiled; rebuild with --features sci-form".to_owned(),
                )
            },
            capabilities: CalculatorCapabilities {
                periodic: capabilities.periodic,
                forces: capabilities.forces,
                stress: capabilities.stress,
                charges: true,
                in_process: true,
                thread_safe: true,
            },
        }
    })
}

fn build_gfn0(_options: CalculatorOptions) -> Result<Box<dyn Calculator>, Error> {
    Ok(Box::new(adapter(QuantumMethod::Gfn0)?))
}

fn build_gfn1(_options: CalculatorOptions) -> Result<Box<dyn Calculator>, Error> {
    Ok(Box::new(adapter(QuantumMethod::Gfn1)?))
}

fn build_gfn2(_options: CalculatorOptions) -> Result<Box<dyn Calculator>, Error> {
    Ok(Box::new(adapter(QuantumMethod::Gfn2)?))
}

/// The registry with every built-in engine.
#[must_use]
pub fn builtin_registry() -> CalculatorRegistry {
    let mut registry = CalculatorRegistry::new();
    registry.add_builtin(uff_info(), build_uff);
    let gfn0 = quantum_info(QuantumMethod::Gfn0, &GFN0_INFO);
    let gfn1 = quantum_info(QuantumMethod::Gfn1, &GFN1_INFO);
    let gfn2 = quantum_info(QuantumMethod::Gfn2, &GFN2_INFO);
    registry.add_builtin(gfn0, build_gfn0);
    registry.add_builtin(gfn1, build_gfn1);
    registry.add_builtin(gfn2, build_gfn2);
    registry
}
