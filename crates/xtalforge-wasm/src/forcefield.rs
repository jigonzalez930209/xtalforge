//! Force-field exports (Task 4.12.2).

use wasm_bindgen::prelude::*;

use xtalforge::ff::Calculator;

use crate::common::{parse_any, set, summary_object, to_js_error};

fn calculator() -> Result<xtalforge::ff::UffCalculator, JsValue> {
    xtalforge::ff::UffCalculator::with_defaults().map_err(ff_error)
}

/// Returns the UFF energy in eV.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or typing.
#[wasm_bindgen]
pub fn ff_energy(text: &str) -> Result<f64, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    calculator()?.energy(&structure).map_err(ff_error)
}

/// Relaxes a structure and returns summary fields.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or an unknown minimizer.
#[wasm_bindgen]
pub fn ff_optimize_info(
    text: &str,
    minimizer: &str,
    max_steps: usize,
    force_tolerance: f64,
    cell: bool,
) -> Result<JsValue, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let minimizer = match minimizer {
        "fire" => xtalforge::ff::Minimizer::Fire,
        "lbfgs" => xtalforge::ff::Minimizer::Lbfgs,
        other => {
            return Err(JsValue::from_str(&format!(
                "invalid_input: unknown minimizer {other:?}"
            )));
        }
    };
    let options = xtalforge::ff::OptimizeOptions {
        max_steps,
        force_tolerance,
        cell,
        ..Default::default()
    };
    let result = xtalforge::ff::minimize(&calculator()?, &structure, &options, minimizer)
        .map_err(ff_error)?;
    let object = summary_object(&result.structure)?;
    let object: js_sys::Object = object.unchecked_into();
    set(&object, "energy", JsValue::from_f64(result.energy));
    set(&object, "max_force", JsValue::from_f64(result.max_force));
    set(&object, "steps", JsValue::from_f64(result.steps as f64));
    set(&object, "converged", JsValue::from_bool(result.converged));
    Ok(object.into())
}

/// Returns QEq or EQeq charges.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or an unknown method.
#[wasm_bindgen]
pub fn ff_charges_info(text: &str, method: &str, total_charge: f64) -> Result<JsValue, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let options = xtalforge::ff::ChargeOptions {
        total_charge,
        ..Default::default()
    };
    let charges = match method {
        "qeq" => xtalforge::ff::qeq_charges(&structure, &options),
        "eqeq" => xtalforge::ff::eqeq_charges(&structure, &options),
        other => {
            return Err(JsValue::from_str(&format!(
                "invalid_input: unknown charge method {other:?}"
            )));
        }
    }
    .map_err(ff_error)?;
    let array = js_sys::Array::new();
    for charge in &charges.per_atom {
        array.push(&JsValue::from_f64(*charge));
    }
    let object = js_sys::Object::new();
    set(&object, "per_atom", array.into());
    set(&object, "total", JsValue::from_f64(charges.total));
    set(
        &object,
        "iterations",
        JsValue::from_f64(charges.iterations as f64),
    );
    Ok(object.into())
}

fn p1_modes(
    text: &str,
    step: f64,
) -> Result<(xtalforge::core::Structure, xtalforge::ff::VibrationalModes), JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let p1 = xtalforge::symmetry::to_p1(&structure).map_err(to_js_error)?;
    let options = xtalforge::ff::HessianOptions {
        step,
        ..Default::default()
    };
    let matrix = xtalforge::ff::hessian(&calculator()?, &p1, &options).map_err(ff_error)?;
    let modes = xtalforge::ff::vibrational_modes(&matrix, &p1, &options).map_err(ff_error)?;
    Ok((p1, modes))
}

/// Returns the Gamma-point frequencies and flags.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or a failed Hessian.
#[wasm_bindgen]
pub fn ff_vibrations_info(text: &str, step: f64) -> Result<JsValue, JsValue> {
    let (_, modes) = p1_modes(text, step)?;
    let frequencies = js_sys::Array::new();
    for frequency in &modes.frequencies {
        frequencies.push(&JsValue::from_f64(*frequency));
    }
    let imaginary = js_sys::Array::new();
    for mode in &modes.imaginary {
        imaginary.push(&JsValue::from_f64(*mode as f64));
    }
    let zero_modes = js_sys::Array::new();
    for mode in &modes.zero_modes {
        zero_modes.push(&JsValue::from_f64(*mode as f64));
    }
    let object = js_sys::Object::new();
    set(&object, "frequencies", frequencies.into());
    set(&object, "imaginary", imaginary.into());
    set(&object, "zero_modes", zero_modes.into());
    Ok(object.into())
}

/// Returns the broadened IR spectrum.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or a failed Hessian.
#[wasm_bindgen]
pub fn ff_ir_info(text: &str, fwhm: f64) -> Result<JsValue, JsValue> {
    let (p1, modes) = p1_modes(text, 0.005)?;
    let spectrum = xtalforge::ff::ir_spectrum(
        &p1,
        &modes,
        &xtalforge::ff::IrOptions {
            fwhm,
            ..Default::default()
        },
    )
    .map_err(ff_error)?;
    let points = js_sys::Array::new();
    for point in &spectrum.points {
        let pair = js_sys::Array::new();
        pair.push(&JsValue::from_f64(point[0]));
        pair.push(&JsValue::from_f64(point[1]));
        points.push(&pair.into());
    }
    let peaks = js_sys::Array::new();
    for peak in &spectrum.peaks {
        let item = js_sys::Object::new();
        set(&item, "mode", JsValue::from_f64(peak.mode as f64));
        set(&item, "frequency", JsValue::from_f64(peak.frequency));
        set(&item, "intensity", JsValue::from_f64(peak.intensity));
        peaks.push(&item.into());
    }
    let object = js_sys::Object::new();
    set(&object, "points", points.into());
    set(&object, "peaks", peaks.into());
    Ok(object.into())
}

/// Returns harmonic thermochemistry.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or a failed Hessian.
#[wasm_bindgen]
pub fn ff_thermo_info(text: &str, temperature: f64, pressure: f64) -> Result<JsValue, JsValue> {
    let (_, modes) = p1_modes(text, 0.005)?;
    let result = xtalforge::ff::thermochemistry(
        &modes.frequencies,
        &xtalforge::ff::ThermoOptions {
            temperature,
            pressure,
        },
    )
    .map_err(ff_error)?;
    let object = js_sys::Object::new();
    set(&object, "zpe", JsValue::from_f64(result.zpe));
    set(&object, "thermal_u", JsValue::from_f64(result.thermal_u));
    set(&object, "thermal_h", JsValue::from_f64(result.thermal_h));
    set(&object, "entropy", JsValue::from_f64(result.entropy));
    set(&object, "gibbs", JsValue::from_f64(result.gibbs));
    Ok(object.into())
}

/// Returns the UFF atom types with reasons.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or typing.
#[wasm_bindgen]
pub fn ff_types_info(text: &str) -> Result<js_sys::Array, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let model = xtalforge::ff::UffModel::build(&structure, Default::default()).map_err(ff_error)?;
    let array = js_sys::Array::new();
    for (atom, typing) in model.typing().types.iter().enumerate() {
        let item = js_sys::Object::new();
        set(&item, "atom", JsValue::from_f64(atom as f64));
        set(&item, "type", JsValue::from_str(&typing.name));
        set(&item, "reason", JsValue::from_str(&typing.reason));
        array.push(&item.into());
    }
    Ok(array)
}

fn ff_error(error: xtalforge::ff::Error) -> JsValue {
    let value = js_sys::Object::new();
    let _ = js_sys::Reflect::set(
        &value,
        &JsValue::from_str("code"),
        &JsValue::from_str(error.code()),
    );
    let _ = js_sys::Reflect::set(
        &value,
        &JsValue::from_str("message"),
        &JsValue::from_str(&error.message()),
    );
    value.into()
}

/// Runs the preparation pipeline and returns its report.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or a failed step.
#[wasm_bindgen]
pub fn ff_prepare_info(
    text: &str,
    minimize: bool,
    charges: bool,
    ewald: bool,
    vibrations: bool,
    merge_duplicates: Option<bool>,
    order_disorder: Option<bool>,
    remove_solvent: Option<bool>,
    check_charge_balance: Option<bool>,
) -> Result<JsValue, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let options = xtalforge::ff::PrepareOptions {
        minimize,
        charges,
        ewald,
        vibrations,
        merge_duplicates: merge_duplicates.unwrap_or(false),
        order_disorder: order_disorder.unwrap_or(false),
        remove_solvent: remove_solvent.unwrap_or(false),
        check_charge_balance: check_charge_balance.unwrap_or(false),
        ..Default::default()
    };
    let report = xtalforge::ff::prepare(&calculator()?, &structure, &options).map_err(ff_error)?;
    let object = js_sys::Object::new();
    set(&object, "energy", JsValue::from_f64(report.energy));
    set(
        &object,
        "total_charge",
        JsValue::from_f64(report.total_charge),
    );
    set(
        &object,
        "ewald_energy",
        report.ewald_energy.map_or(JsValue::NULL, JsValue::from_f64),
    );
    let frequencies = js_sys::Array::new();
    for frequency in &report.frequencies {
        frequencies.push(&JsValue::from_f64(*frequency));
    }
    set(&object, "frequencies", frequencies.into());
    set(
        &object,
        "removed_sites",
        JsValue::from_f64(report.removed_sites as f64),
    );
    if let Some(balance) = &report.charge_balance {
        let item = js_sys::Object::new();
        set(
            &item,
            "total_charge",
            JsValue::from_f64(balance.total_charge),
        );
        set(&item, "neutral", JsValue::from_bool(balance.neutral));
        set(&object, "charge_balance", item.into());
    } else {
        set(&object, "charge_balance", JsValue::NULL);
    }
    let structure = serde_json::to_string(&report.structure)
        .map_err(|error| to_js_error(xtalforge::ff::Error::InvalidInput(error.to_string())))?;
    set(&object, "structure", JsValue::from_str(&structure));
    let steps = js_sys::Array::new();
    for step in &report.steps {
        let item = js_sys::Object::new();
        set(&item, "name", JsValue::from_str(&step.name));
        set(&item, "done", JsValue::from_bool(step.done));
        set(&item, "detail", JsValue::from_str(&step.detail));
        steps.push(&item.into());
    }
    set(&object, "steps", steps.into());
    Ok(object.into())
}
