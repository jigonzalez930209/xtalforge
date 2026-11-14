//! Phase 19 bindings: magnetism, disorder, diffraction and spectroscopy.

use pyo3::prelude::*;

use xtalforge::analysis::{
    DerivativeOptions, MagneticOptions, MagneticOrderingOptions, NeutronOptions, SaedOptions,
    derivative_structures, disorder_composition, disorder_groups, magnetic_diffraction,
    magnetic_orderings, neutron, zone_axis_pattern,
};

fn json_value(py: Python<'_>, value: &impl serde::Serialize) -> PyResult<Py<PyAny>> {
    let text = serde_json::to_string(value)
        .map_err(|error| crate::invalid_input(format!("phase 19 JSON: {error}")))?;
    let module = py.import("json")?;
    module
        .call_method1("loads", (text,))
        .map(|value| value.unbind())
}

/// Enumerates symmetry-ranked collinear magnetic orderings.
#[pyfunction(name = "magnetic_orderings")]
#[pyo3(signature = (text, supercell = 1, max_orderings = 8, format_hint = None))]
pub(crate) fn magnetic_orderings_py(
    py: Python<'_>,
    text: &str,
    supercell: u32,
    max_orderings: usize,
    format_hint: Option<&str>,
) -> PyResult<Py<PyAny>> {
    let structure = crate::structures::load_structure(text, format_hint)?;
    let options = MagneticOrderingOptions {
        supercell: [supercell, supercell, supercell],
        max_orderings,
        ..MagneticOrderingOptions::default()
    };
    let result = magnetic_orderings(&structure, &options).map_err(crate::map_error)?;
    json_value(py, &result)
}

/// Reports disorder groups, composition and density.
#[pyfunction(name = "disorder_groups")]
#[pyo3(signature = (text, format_hint = None))]
pub(crate) fn disorder_groups_py(
    py: Python<'_>,
    text: &str,
    format_hint: Option<&str>,
) -> PyResult<Py<PyAny>> {
    let structure = crate::structures::load_structure(text, format_hint)?;
    let groups = disorder_groups(&structure).map_err(crate::map_error)?;
    let composition = disorder_composition(&structure).map_err(crate::map_error)?;
    let value = serde_json::json!({
        "formula": composition.formula(),
        "density": composition.density(structure.lattice().volume()).map_err(crate::map_error)?,
        "groups": groups,
    });
    json_value(py, &value)
}

/// Enumerates symmetry-distinct derivative structures.
#[pyfunction(name = "derivative_structures")]
#[pyo3(signature = (text, max_index = 2, max_structures = 8, format_hint = None))]
pub(crate) fn derivative_structures_py(
    py: Python<'_>,
    text: &str,
    max_index: u32,
    max_structures: usize,
    format_hint: Option<&str>,
) -> PyResult<Py<PyAny>> {
    let structure = crate::structures::load_structure(text, format_hint)?;
    let options = DerivativeOptions {
        max_index,
        max_structures,
        ..DerivativeOptions::default()
    };
    let result = derivative_structures(&structure, &options).map_err(crate::map_error)?;
    json_value(py, &result)
}

/// Calculates a nuclear neutron powder pattern.
#[pyfunction]
#[pyo3(signature = (text, wavelength = 1.54, step = 0.05, format_hint = None))]
pub(crate) fn neutron_powder(
    py: Python<'_>,
    text: &str,
    wavelength: f64,
    step: f64,
    format_hint: Option<&str>,
) -> PyResult<Py<PyAny>> {
    let structure = crate::structures::load_structure(text, format_hint)?;
    let options = NeutronOptions {
        geometry: xtalforge::analysis::NeutronGeometry::ConstantWavelength { wavelength },
        step,
        ..NeutronOptions::default()
    };
    let pattern = neutron::neutron_powder(&structure, &options).map_err(crate::map_error)?;
    json_value(py, &pattern)
}

/// Calculates the magnetic neutron powder pattern.
#[pyfunction]
#[pyo3(signature = (text, wavelength = 1.54, step = 0.05, allow_missing_moments = false, format_hint = None))]
pub(crate) fn magnetic_powder(
    py: Python<'_>,
    text: &str,
    wavelength: f64,
    step: f64,
    allow_missing_moments: bool,
    format_hint: Option<&str>,
) -> PyResult<Py<PyAny>> {
    let structure = crate::structures::load_structure(text, format_hint)?;
    let options = MagneticOptions {
        profile: xtalforge::analysis::InstrumentProfile {
            radiation: xtalforge::analysis::RadiationSource {
                kind: xtalforge::analysis::RadiationKind::Neutron,
                wavelength,
            },
            ..xtalforge::analysis::InstrumentProfile::default()
        },
        step,
        allow_missing_moments,
        ..MagneticOptions::default()
    };
    let pattern =
        magnetic_diffraction::magnetic_powder(&structure, &options).map_err(crate::map_error)?;
    json_value(py, &pattern)
}

/// Calculates a kinematic SAED pattern for a zone axis.
#[pyfunction]
#[pyo3(signature = (text, zone = (0, 0, 1), voltage = 200.0, size = 256, format_hint = None))]
pub(crate) fn electron_diffraction(
    py: Python<'_>,
    text: &str,
    zone: (i32, i32, i32),
    voltage: f64,
    size: usize,
    format_hint: Option<&str>,
) -> PyResult<Py<PyAny>> {
    let structure = crate::structures::load_structure(text, format_hint)?;
    let options = SaedOptions {
        acceleration_kv: voltage,
        zone_axis: [zone.0, zone.1, zone.2],
        image_size: size,
        ..SaedOptions::default()
    };
    let pattern = zone_axis_pattern(&structure, &options).map_err(crate::map_error)?;
    let value = serde_json::json!({
        "spots": pattern.spots,
        "max_intensity": pattern.max_intensity(),
    });
    json_value(py, &value)
}

/// Reads Raman tensors and returns the powder average per tensor.
#[pyfunction]
pub(crate) fn raman_powder_average(tensors: Vec<Vec<Vec<f64>>>) -> PyResult<Vec<f64>> {
    let mut values = Vec::with_capacity(tensors.len());
    for tensor in tensors {
        let matrix = tensor_to_matrix(&tensor)?;
        values.push(xtalforge::ff::powder_average(&matrix));
    }
    Ok(values)
}

/// Reads Born tensors and returns the isotropic diagonal average per atom.
#[pyfunction]
pub(crate) fn born_diagonal_average(tensors: Vec<Vec<Vec<f64>>>) -> PyResult<Vec<f64>> {
    let mut values = Vec::with_capacity(tensors.len());
    for tensor in tensors {
        let matrix = tensor_to_matrix(&tensor)?;
        values.push((matrix[0][0] + matrix[1][1] + matrix[2][2]) / 3.0);
    }
    Ok(values)
}

fn tensor_to_matrix(tensor: &[Vec<f64>]) -> PyResult<[[f64; 3]; 3]> {
    if tensor.len() != 3 || tensor.iter().any(|row| row.len() != 3) {
        return Err(crate::invalid_input("a tensor must be a 3x3 matrix"));
    }
    let mut matrix = [[0.0; 3]; 3];
    for (target, values) in matrix.iter_mut().zip(tensor.iter()) {
        for (slot, value) in target.iter_mut().zip(values.iter()) {
            *slot = *value;
        }
    }
    Ok(matrix)
}

/// Registers the Phase 19 functions on the extension module.
pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(magnetic_orderings_py, module)?)?;
    module.add_function(wrap_pyfunction!(disorder_groups_py, module)?)?;
    module.add_function(wrap_pyfunction!(derivative_structures_py, module)?)?;
    module.add_function(wrap_pyfunction!(neutron_powder, module)?)?;
    module.add_function(wrap_pyfunction!(magnetic_powder, module)?)?;
    module.add_function(wrap_pyfunction!(electron_diffraction, module)?)?;
    module.add_function(wrap_pyfunction!(raman_powder_average, module)?)?;
    module.add_function(wrap_pyfunction!(born_diagonal_average, module)?)?;
    Ok(())
}
