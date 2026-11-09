//! Force-field bindings (Task 4.12.2).

use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

use xtalforge::ff::{
    Calculator, CalculatorOptions, ChargeOptions, HessianOptions, IrOptions, Minimizer,
    OptimizeOptions, ThermoOptions, UffCalculator, UffModel, eqeq_charges, hessian, ir_spectrum,
    minimize, qeq_charges, thermochemistry, vibrational_modes,
};

use crate::structures::{load_structure, structure_payload};

fn calculator() -> PyResult<UffCalculator> {
    UffCalculator::with_defaults().map_err(crate::map_ff_error)
}

/// Returns the UFF energy in eV.
#[pyfunction]
#[pyo3(signature = (text, format_hint = None))]
pub(crate) fn ff_energy(text: &str, format_hint: Option<&str>) -> PyResult<f64> {
    let structure = load_structure(text, format_hint)?;
    calculator()?
        .energy(&structure)
        .map_err(crate::map_ff_error)
}

/// Relaxes a structure and returns its payload plus the run summary.
#[pyfunction]
#[pyo3(signature = (text, minimizer = "fire".to_string(), max_steps = 500, force_tolerance = 0.05, cell = false, format_hint = None))]
#[allow(clippy::too_many_arguments)] // one Python argument per option
pub(crate) fn ff_optimize<'py>(
    py: Python<'py>,
    text: &str,
    minimizer: String,
    max_steps: usize,
    force_tolerance: f64,
    cell: bool,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyDict>> {
    let structure = load_structure(text, format_hint)?;
    let minimizer = match minimizer.as_str() {
        "fire" => Minimizer::Fire,
        "lbfgs" => Minimizer::Lbfgs,
        other => {
            return Err(crate::invalid_input(format!(
                "unknown minimizer {other:?}; use fire or lbfgs"
            )));
        }
    };
    let options = OptimizeOptions {
        max_steps,
        force_tolerance,
        cell,
        ..OptimizeOptions::default()
    };
    let result =
        minimize(&calculator()?, &structure, &options, minimizer).map_err(crate::map_ff_error)?;
    let dict = structure_payload(py, &result.structure)?;
    dict.set_item("energy", result.energy)?;
    dict.set_item("max_force", result.max_force)?;
    dict.set_item("steps", result.steps)?;
    dict.set_item("converged", result.converged)?;
    Ok(dict)
}

/// Returns QEq or EQeq charges.
#[pyfunction]
#[pyo3(signature = (text, method = "qeq".to_string(), total_charge = 0.0, format_hint = None))]
pub(crate) fn ff_charges<'py>(
    py: Python<'py>,
    text: &str,
    method: String,
    total_charge: f64,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyDict>> {
    let structure = load_structure(text, format_hint)?;
    let options = ChargeOptions {
        total_charge,
        ..ChargeOptions::default()
    };
    let charges = match method.as_str() {
        "qeq" => qeq_charges(&structure, &options),
        "eqeq" => eqeq_charges(&structure, &options),
        other => {
            return Err(crate::invalid_input(format!(
                "unknown charge method {other:?}; use qeq or eqeq"
            )));
        }
    }
    .map_err(crate::map_ff_error)?;
    let dict = PyDict::new(py);
    dict.set_item("per_atom", charges.per_atom)?;
    dict.set_item("total", charges.total)?;
    dict.set_item("iterations", charges.iterations)?;
    Ok(dict)
}

/// Returns the Gamma-point frequencies, zero modes and imaginary flags.
#[pyfunction]
#[pyo3(signature = (text, step = 0.005, format_hint = None))]
pub(crate) fn ff_vibrations<'py>(
    py: Python<'py>,
    text: &str,
    step: f64,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyDict>> {
    let structure = load_structure(text, format_hint)?;
    let p1 = xtalforge::symmetry::to_p1(&structure).map_err(crate::map_error)?;
    let options = HessianOptions {
        step,
        ..HessianOptions::default()
    };
    let matrix = hessian(&calculator()?, &p1, &options).map_err(crate::map_ff_error)?;
    let modes = vibrational_modes(&matrix, &p1, &options).map_err(crate::map_ff_error)?;
    let dict = PyDict::new(py);
    dict.set_item("frequencies", modes.frequencies)?;
    dict.set_item("zero_modes", modes.zero_modes)?;
    dict.set_item("imaginary", modes.imaginary)?;
    let vectors = PyList::empty(py);
    for mode in &modes.modes {
        vectors.append(mode.to_vec())?;
    }
    dict.set_item("modes", vectors)?;
    Ok(dict)
}

/// Returns the broadened IR spectrum.
#[pyfunction]
#[pyo3(signature = (text, fwhm = 10.0, format_hint = None))]
pub(crate) fn ff_ir<'py>(
    py: Python<'py>,
    text: &str,
    fwhm: f64,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyDict>> {
    let structure = load_structure(text, format_hint)?;
    let p1 = xtalforge::symmetry::to_p1(&structure).map_err(crate::map_error)?;
    let options = HessianOptions::default();
    let matrix = hessian(&calculator()?, &p1, &options).map_err(crate::map_ff_error)?;
    let modes = vibrational_modes(&matrix, &p1, &options).map_err(crate::map_ff_error)?;
    let spectrum = ir_spectrum(
        &p1,
        &modes,
        &IrOptions {
            fwhm,
            ..IrOptions::default()
        },
    )
    .map_err(crate::map_ff_error)?;
    let dict = PyDict::new(py);
    dict.set_item("points", spectrum.points)?;
    let peaks = PyList::empty(py);
    for peak in &spectrum.peaks {
        let item = PyDict::new(py);
        item.set_item("mode", peak.mode)?;
        item.set_item("frequency", peak.frequency)?;
        item.set_item("intensity", peak.intensity)?;
        peaks.append(item)?;
    }
    dict.set_item("peaks", peaks)?;
    Ok(dict)
}

/// Returns harmonic thermochemistry.
#[pyfunction]
#[pyo3(signature = (text, temperature = 298.15, pressure = 1.0, format_hint = None))]
pub(crate) fn ff_thermo<'py>(
    py: Python<'py>,
    text: &str,
    temperature: f64,
    pressure: f64,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyDict>> {
    let frequencies = ff_vibrations(py, text, 0.005, format_hint)?
        .get_item("frequencies")?
        .ok_or_else(|| crate::invalid_input("no frequencies"))?;
    let frequencies: Vec<f64> = frequencies.extract()?;
    let result = thermochemistry(
        &frequencies,
        &ThermoOptions {
            temperature,
            pressure,
        },
    )
    .map_err(crate::map_ff_error)?;
    let dict = PyDict::new(py);
    dict.set_item("zpe", result.zpe)?;
    dict.set_item("thermal_u", result.thermal_u)?;
    dict.set_item("thermal_h", result.thermal_h)?;
    dict.set_item("entropy", result.entropy)?;
    dict.set_item("gibbs", result.gibbs)?;
    Ok(dict)
}

/// Returns the UFF atom types with their reasons.
#[pyfunction]
#[pyo3(signature = (text, format_hint = None))]
pub(crate) fn ff_types<'py>(
    py: Python<'py>,
    text: &str,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyList>> {
    let structure = load_structure(text, format_hint)?;
    let model =
        UffModel::build(&structure, CalculatorOptions::default()).map_err(crate::map_ff_error)?;
    let items = PyList::empty(py);
    for (atom, typing) in model.typing().types.iter().enumerate() {
        let dict = PyDict::new(py);
        dict.set_item("atom", atom)?;
        dict.set_item("type", &typing.name)?;
        dict.set_item("reason", &typing.reason)?;
        dict.set_item("confidence", format!("{:?}", typing.confidence))?;
        items.append(dict)?;
    }
    Ok(items)
}

/// Runs the preparation pipeline and returns its report.
#[pyfunction]
#[pyo3(signature = (text, minimize = true, charges = true, ewald = false, vibrations = false, format_hint = None, merge_duplicates = false, order_disorder = false, remove_solvent = false, check_charge_balance = false))]
#[allow(clippy::too_many_arguments)]
pub(crate) fn ff_prepare<'py>(
    py: Python<'py>,
    text: &str,
    minimize: bool,
    charges: bool,
    ewald: bool,
    vibrations: bool,
    format_hint: Option<&str>,
    merge_duplicates: bool,
    order_disorder: bool,
    remove_solvent: bool,
    check_charge_balance: bool,
) -> PyResult<Bound<'py, PyDict>> {
    let structure = load_structure(text, format_hint)?;
    let options = xtalforge::ff::PrepareOptions {
        minimize,
        charges,
        ewald,
        vibrations,
        merge_duplicates,
        order_disorder,
        remove_solvent,
        check_charge_balance,
        ..Default::default()
    };
    let report = xtalforge::ff::prepare(&calculator()?, &structure, &options)
        .map_err(crate::map_ff_error)?;
    let dict = PyDict::new(py);
    dict.set_item("energy", report.energy)?;
    dict.set_item("total_charge", report.total_charge)?;
    dict.set_item("ewald_energy", report.ewald_energy)?;
    dict.set_item("frequencies", report.frequencies)?;
    dict.set_item("removed_sites", report.removed_sites)?;
    if let Some(balance) = &report.charge_balance {
        let item = PyDict::new(py);
        item.set_item("total_charge", balance.total_charge)?;
        item.set_item("neutral", balance.neutral)?;
        dict.set_item("charge_balance", item)?;
    } else {
        dict.set_item("charge_balance", py.None())?;
    }
    dict.set_item(
        "structure",
        serde_json::to_string(&report.structure)
            .map_err(|error| crate::invalid_input(error.to_string()))?,
    )?;
    let steps = PyList::empty(py);
    for step in &report.steps {
        let item = PyDict::new(py);
        item.set_item("name", &step.name)?;
        item.set_item("done", step.done)?;
        item.set_item("detail", &step.detail)?;
        steps.append(item)?;
    }
    dict.set_item("steps", steps)?;
    Ok(dict)
}
