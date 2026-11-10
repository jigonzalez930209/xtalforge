//! Machine-learning potential bridge (Task 8.4.1).
//!
//! The decision recorded in Task 4.11 is the Python bridge: MACE, ORB-v3 and
//! MatterSim are torch models whose graph is not a fixed ONNX operator set, so
//! they run through ASE in Python rather than through a native Rust backend.
//!
//! [`MlCalculator`] wraps **any** ASE calculator. The caller builds the
//! backend object (`mace.calculators.mace_mp`, an ORB calculator, a MatterSim
//! calculator) and passes it in; this class only converts structure JSON to an
//! ASE `Atoms`, runs the calculator and maps energy, forces and stress back.
//!
//! Numbers are eV, eV/Angstrom and GPa at this boundary, like every other
//! force-field surface. ASE reports stress in eV/Angstrom^3; the conversion is
//! applied here and the mapping from its six components to a 3x3 tensor is
//! documented on [`MlCalculator::energy_forces`].

use pyo3::prelude::*;
use serde_json::{Value, json};

use crate::XtalforgeError;

/// `1 eV/Angstrom^3` in GPa.
const EV_PER_A3_TO_GPA: f64 = 160.217_662_08;

/// Reports which Python ML modules are importable.
#[pyfunction]
pub(crate) fn ml_backends(py: Python<'_>) -> String {
    let mut report = serde_json::Map::new();
    for module in ["torch", "ase", "mace", "orb_models", "mattersim"] {
        let version = py.import(module).ok().and_then(|module| {
            module
                .getattr("__version__")
                .ok()
                .and_then(|value| value.extract::<String>().ok())
        });
        report.insert(
            module.to_owned(),
            json!({
                "available": version.is_some() || py.import(module).is_ok(),
                "version": version,
            }),
        );
    }
    Value::Object(report).to_string()
}

fn backend_error(error: PyErr) -> PyErr {
    XtalforgeError::new_err(("ml_backend", error.to_string()))
}

fn parse_structure(text: &str) -> PyResult<xtalforge::core::Structure> {
    let value: Value = serde_json::from_str(text)
        .map_err(|error| XtalforgeError::new_err(("json", error.to_string())))?;
    let inner = value.get("structure").cloned().unwrap_or(value);
    serde_json::from_value(inner)
        .map_err(|error| XtalforgeError::new_err(("json", error.to_string())))
}

fn atoms<'py>(
    py: Python<'py>,
    structure: &xtalforge::core::Structure,
    periodic: bool,
) -> PyResult<Bound<'py, PyAny>> {
    let p1 = structure
        .p1()
        .map_err(|error| XtalforgeError::new_err((error.code(), error.to_string())))?;
    let symbols: Vec<String> = p1.atoms.iter().map(|atom| atom.element.clone()).collect();
    let positions: Vec<[f64; 3]> = p1
        .atoms
        .iter()
        .map(|atom| p1.lattice.to_cart(atom.frac))
        .collect();
    let cell = p1.lattice.matrix().to_vec();
    let ase = py.import("ase").map_err(backend_error)?;
    let arguments = pyo3::types::PyDict::new(py);
    arguments
        .set_item("symbols", symbols)
        .map_err(backend_error)?;
    arguments
        .set_item("positions", positions)
        .map_err(backend_error)?;
    arguments.set_item("cell", cell).map_err(backend_error)?;
    arguments.set_item("pbc", periodic).map_err(backend_error)?;
    ase.getattr("Atoms")
        .map_err(backend_error)?
        .call((), Some(&arguments))
        .map_err(backend_error)
}

/// An ASE calculator exposed through the xtalforge boundary.
#[pyclass]
pub struct MlCalculator {
    calculator: Py<PyAny>,
    periodic: bool,
    name: String,
}

#[pymethods]
impl MlCalculator {
    /// Wraps an ASE calculator.
    #[new]
    #[pyo3(signature = (calculator, periodic = true))]
    fn new(calculator: Py<PyAny>, periodic: bool) -> Self {
        let name = Python::attach(|py| {
            calculator
                .bind(py)
                .getattr("name")
                .ok()
                .and_then(|value| value.extract::<String>().ok())
                .unwrap_or_else(|| "ase".to_owned())
        });
        Self {
            calculator,
            periodic,
            name,
        }
    }

    /// Name reported by the ASE calculator.
    fn name(&self) -> String {
        self.name.clone()
    }

    /// Energy in eV for one structure JSON document.
    fn energy(&self, py: Python<'_>, structure_json: &str) -> PyResult<f64> {
        let structure = parse_structure(structure_json)?;
        let atoms = atoms(py, &structure, self.periodic)?;
        atoms
            .call_method1("set_calculator", (self.calculator.bind(py),))
            .map_err(backend_error)?;
        atoms
            .call_method0("get_potential_energy")
            .map_err(backend_error)?
            .extract::<f64>()
            .map_err(backend_error)
    }

    /// Energy in eV, forces in eV/Angstrom and stress in GPa as JSON.
    ///
    /// The stress tensor is assembled from ASE's six components
    /// `(xx, yy, zz, yz, xz, xy)` into `[[xx, xy, xz], [xy, yy, yz],
    /// [xz, yz, zz]]` and converted from `eV/Angstrom^3` to GPa. It is `null`
    /// when the calculator does not provide a stress.
    fn energy_forces(&self, py: Python<'_>, structure_json: &str) -> PyResult<String> {
        let structure = parse_structure(structure_json)?;
        let atoms = atoms(py, &structure, self.periodic)?;
        atoms
            .call_method1("set_calculator", (self.calculator.bind(py),))
            .map_err(backend_error)?;
        let energy = atoms
            .call_method0("get_potential_energy")
            .map_err(backend_error)?
            .extract::<f64>()
            .map_err(backend_error)?;
        let forces = atoms
            .call_method0("get_forces")
            .map_err(backend_error)?
            .extract::<Vec<Vec<f64>>>()
            .map_err(backend_error)?;
        let stress = atoms
            .call_method0("get_stress")
            .ok()
            .and_then(|value| value.extract::<Vec<f64>>().ok())
            .and_then(|values| {
                let [xx, yy, zz, yz, xz, xy] = values.as_slice() else {
                    return None;
                };
                let scale = EV_PER_A3_TO_GPA;
                Some([
                    [xx * scale, xy * scale, xz * scale],
                    [xy * scale, yy * scale, yz * scale],
                    [xz * scale, yz * scale, zz * scale],
                ])
            });
        let result = json!({
            "energy_ev": energy,
            "forces": forces,
            "stress_gpa": stress,
            "n_atoms": forces.len(),
        });
        serde_json::to_string(&result)
            .map_err(|error| XtalforgeError::new_err(("json", error.to_string())))
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(ml_backends, module)?)?;
    module.add_class::<MlCalculator>()?;
    Ok(())
}
