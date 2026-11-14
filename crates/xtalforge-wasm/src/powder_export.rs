//! Powder refinement exports (Tasks 7.1-7.11).

use wasm_bindgen::prelude::*;

use crate::common::{array_i32, parse_any, set, to_js_error};

fn pattern(text: &str) -> Result<Vec<xtalforge::io::PatternPoint>, JsValue> {
    xtalforge::io::read_xy(text).map_err(to_js_error)
}

/// Runs the automatic peak/index/Le Bail/Rietveld pipeline.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text, pattern or options.
#[wasm_bindgen]
pub fn powder_automatic_info(text: &str, pattern_text: &str) -> Result<JsValue, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let observed = pattern(pattern_text)?;
    let result = xtalforge::analysis::run_automatic(
        &structure,
        &observed,
        &xtalforge::analysis::AutomaticOptions::default(),
    )
    .map_err(to_js_error)?;
    let steps = js_sys::Array::new();
    for step in &result.steps {
        steps.push(&JsValue::from_str(step));
    }
    let object = js_sys::Object::new();
    set(
        &object,
        "peaks",
        JsValue::from_f64(result.peaks.len() as f64),
    );
    set(&object, "steps", steps.into());
    set(
        &object,
        "cell",
        result
            .index
            .map_or(JsValue::NULL, |index| JsValue::from_f64(index.cell.a)),
    );
    // `cell` is only the edge of the quick cubic index (kept for the 1.0
    // signature); `best_cell` is the Pawley-refined cell of the leading
    // (cell, extinction class) candidate, `[a, b, c, alpha, beta, gamma]`.
    let best = js_sys::Array::new();
    if let Some(pawley) = &result.pawley {
        for value in pawley.cell {
            best.push(&JsValue::from_f64(value));
        }
    } else if let Some(row) = result.leading_cells.first() {
        let cell = &row.index.cell;
        for value in [cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma] {
            best.push(&JsValue::from_f64(value));
        }
    }
    set(
        &object,
        "best_cell",
        if best.length() == 0 {
            JsValue::NULL
        } else {
            best.into()
        },
    );
    set(
        &object,
        "rwp",
        result
            .phase
            .map_or(JsValue::NULL, |phase| JsValue::from_f64(phase.rwp)),
    );
    set(
        &object,
        "pawley_rwp",
        result
            .pawley
            .map_or(JsValue::NULL, |pawley| JsValue::from_f64(pawley.rwp)),
    );
    let leading = js_sys::Array::new();
    for row in &result.leading_cells {
        let item = js_sys::Object::new();
        set(&item, "system", JsValue::from_str(row.index.system.name()));
        set(&item, "cell_rank", JsValue::from_f64(row.cell_rank as f64));
        set(
            &item,
            "class_rank",
            JsValue::from_f64(row.class_rank as f64),
        );
        set(&item, "class_symbol", JsValue::from_str(&row.class_symbol));
        set(&item, "delta_bic", JsValue::from_f64(row.delta_bic));
        for (key, value) in [
            ("a", row.index.cell.a),
            ("b", row.index.cell.b),
            ("c", row.index.cell.c),
            ("alpha", row.index.cell.alpha),
            ("beta", row.index.cell.beta),
            ("gamma", row.index.cell.gamma),
            ("m20", row.index.m20),
            ("f_n", row.index.f_n),
            ("coverage", row.index.coverage),
        ] {
            set(&item, key, JsValue::from_f64(value));
        }
        set(
            &item,
            "pawley_rwp",
            row.pawley_rwp.map_or(JsValue::NULL, JsValue::from_f64),
        );
        let esds = js_sys::Array::new();
        if let Some(values) = row.pawley_cell_esds {
            for value in values {
                esds.push(&JsValue::from_f64(value));
            }
        }
        set(&item, "pawley_cell_esds", esds.into());
        let groups = js_sys::Array::new();
        for group in &row.space_groups {
            let entry = js_sys::Object::new();
            set(&entry, "number", JsValue::from_f64(f64::from(group.number)));
            set(&entry, "symbol", JsValue::from_str(&group.symbol));
            set(
                &entry,
                "centering",
                JsValue::from_str(&group.centering.to_string()),
            );
            set(&entry, "point_group", JsValue::from_str(&group.point_group));
            set(&entry, "order", JsValue::from_f64(group.order as f64));
            set(
                &entry,
                "unexplained_peaks",
                JsValue::from_f64(group.unexplained_peaks as f64),
            );
            groups.push(&entry.into());
        }
        set(&item, "space_groups", groups.into());
        leading.push(&item.into());
    }
    set(&object, "leading_cells", leading.into());
    Ok(object.into())
}

/// Refines a phase with the McCusker stages.
///
/// The March-Dollase axis is the `c` axis. Returns `rwp`, `converged`,
/// `iterations`, `stages`, `parameters` and `parameter_esds`.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or a failed fit.
#[wasm_bindgen]
#[allow(clippy::too_many_arguments)] // one boolean per documented stage switch
pub fn powder_refine_info(
    text: &str,
    pattern_text: &str,
    refine_positions: bool,
    refine_adps: bool,
    refine_occupancies: bool,
    refine_zero: bool,
    refine_cell: bool,
    refine_width: bool,
    refine_texture: bool,
    doublet: bool,
) -> Result<JsValue, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let observed = pattern(pattern_text)?;
    let result = xtalforge::analysis::refine_phase(
        &structure,
        &observed,
        &xtalforge::analysis::RietveldOptions::default(),
        &xtalforge::analysis::RefineOptions {
            refine_positions,
            refine_adps,
            refine_occupancies,
            refine_zero,
            refine_cell,
            refine_width,
            refine_texture,
            doublet,
            ..Default::default()
        },
    )
    .map_err(to_js_error)?;
    let object = js_sys::Object::new();
    set(&object, "rwp", JsValue::from_f64(result.rwp));
    set(&object, "converged", JsValue::from_bool(result.converged));
    set(
        &object,
        "iterations",
        JsValue::from_f64(result.iterations as f64),
    );
    let stages = js_sys::Array::new();
    for stage in &result.stages {
        stages.push(&JsValue::from_str(stage));
    }
    set(&object, "stages", stages.into());
    let parameters = js_sys::Object::new();
    for (name, value) in &result.parameters {
        set(&parameters, name, JsValue::from_f64(*value));
    }
    set(&object, "parameters", parameters.into());
    let esds = js_sys::Object::new();
    for (name, value) in &result.parameter_esds {
        set(&esds, name, JsValue::from_f64(*value));
    }
    set(&object, "parameter_esds", esds.into());
    Ok(object.into())
}

/// Returns the free cell parameters and extinction symbol.
///
/// # Errors
///
/// Returns a JavaScript `Error` for an invalid space-group number.
#[wasm_bindgen]
pub fn powder_cells_info(number: u32) -> Result<JsValue, JsValue> {
    let parameters = xtalforge::analysis::free_cell_parameters(number).map_err(to_js_error)?;
    let extinction = xtalforge::analysis::extinction_symbol(number).map_err(to_js_error)?;
    let array = js_sys::Array::new();
    for parameter in parameters {
        array.push(&JsValue::from_str(parameter));
    }
    let object = js_sys::Object::new();
    set(&object, "free_parameters", array.into());
    set(&object, "extinction", JsValue::from_str(extinction));
    Ok(object.into())
}

/// Returns the rietx bridge inventory.
#[wasm_bindgen]
pub fn rietx_inventory_info() -> Result<JsValue, JsValue> {
    let bridge = xtalforge::analysis::rietx_inventory();
    let object = js_sys::Object::new();
    set(&object, "name", JsValue::from_str(bridge.name));
    set(&object, "license", JsValue::from_str(bridge.license));
    set(&object, "source", JsValue::from_str(bridge.source));
    let capabilities = js_sys::Array::new();
    for capability in bridge.capabilities {
        capabilities.push(&JsValue::from_str(capability));
    }
    set(&object, "capabilities", capabilities.into());
    Ok(object.into())
}

/// Compares the cubic and general indexing engines on a measured XY pattern.
///
/// Returns `{ agree, cubic_a, general_a, delta_a, relative_delta,
/// relative_tolerance, cubic_completeness, general_completeness, n_peaks,
/// message }`.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or a failed search.
#[wasm_bindgen]
pub fn powder_index_agreement_info(
    pattern_text: &str,
    wavelength: f64,
    max_index: i32,
    tolerance: f64,
    relative_tolerance: f64,
) -> Result<JsValue, JsValue> {
    let observed = pattern(pattern_text)?;
    let peaks =
        xtalforge::analysis::find_peaks(&observed, &xtalforge::analysis::PeakOptions::default())
            .map_err(to_js_error)?;
    let result = xtalforge::analysis::indexing_agreement(
        &peaks,
        wavelength,
        &xtalforge::analysis::IndexingOptions {
            max_index,
            tolerance,
            ..xtalforge::analysis::IndexingOptions::default()
        },
        relative_tolerance,
    )
    .map_err(to_js_error)?;
    let object = js_sys::Object::new();
    set(&object, "agree", JsValue::from_bool(result.agree));
    for (key, value) in [
        ("cubic_a", result.cubic_a),
        ("general_a", result.general_a),
        ("delta_a", result.delta_a),
        ("relative_delta", result.relative_delta),
    ] {
        set(&object, key, value.map_or(JsValue::NULL, JsValue::from_f64));
    }
    for (key, value) in [
        ("relative_tolerance", result.relative_tolerance),
        ("cubic_completeness", result.cubic_completeness),
        ("general_completeness", result.general_completeness),
    ] {
        set(&object, key, JsValue::from_f64(value));
    }
    set(&object, "n_peaks", JsValue::from_f64(result.n_peaks as f64));
    set(&object, "message", JsValue::from_str(&result.message));
    Ok(object.into())
}

/// Finds peaks and indexes the metric families of a measured XY pattern.
///
/// Returns `{ peaks: [{ two_theta, intensity, fwhm }], candidates: [{ system,
/// a, b, c, alpha, beta, gamma, m20, f_n, coverage, completeness }] }`.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or a failed search.
#[wasm_bindgen]
pub fn powder_index_info(
    pattern_text: &str,
    wavelength: f64,
    max_index: i32,
    tolerance: f64,
) -> Result<JsValue, JsValue> {
    let observed = pattern(pattern_text)?;
    let peaks =
        xtalforge::analysis::find_peaks(&observed, &xtalforge::analysis::PeakOptions::default())
            .map_err(to_js_error)?;
    let result = xtalforge::analysis::index_pattern(
        &peaks,
        wavelength,
        &xtalforge::analysis::IndexingOptions {
            max_index,
            tolerance,
            ..xtalforge::analysis::IndexingOptions::default()
        },
    )
    .map_err(to_js_error)?;
    let peak_array = js_sys::Array::new();
    for peak in &peaks {
        let item = js_sys::Object::new();
        set(&item, "two_theta", JsValue::from_f64(peak.two_theta));
        set(&item, "intensity", JsValue::from_f64(peak.intensity));
        set(&item, "fwhm", JsValue::from_f64(peak.fwhm));
        peak_array.push(&item.into());
    }
    let candidate_array = js_sys::Array::new();
    for candidate in &result.candidates {
        let item = js_sys::Object::new();
        set(&item, "system", JsValue::from_str(candidate.system.name()));
        for (key, value) in [
            ("a", candidate.cell.a),
            ("b", candidate.cell.b),
            ("c", candidate.cell.c),
            ("alpha", candidate.cell.alpha),
            ("beta", candidate.cell.beta),
            ("gamma", candidate.cell.gamma),
            ("m20", candidate.m20),
            ("f_n", candidate.f_n),
            ("coverage", candidate.coverage),
            ("completeness", candidate.completeness),
        ] {
            set(&item, key, JsValue::from_f64(value));
        }
        candidate_array.push(&item.into());
    }
    let object = js_sys::Object::new();
    set(&object, "peaks", peak_array.into());
    set(&object, "candidates", candidate_array.into());
    Ok(object.into())
}

/// Fits a Pawley model: free intensities, cell, zero shift and width model.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or a failed fit.
#[wasm_bindgen]
pub fn powder_pawley_info(
    text: &str,
    pattern_text: &str,
    profile: &str,
) -> Result<JsValue, JsValue> {
    pawley_info(text, pattern_text, profile, false, false)
}

/// Pawley fit that can also refine the size (constant) and strain
/// (`tan theta`) width terms, like Python `powder_pawley(refine_size=...,
/// refine_strain=...)`.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or a failed fit.
#[wasm_bindgen]
pub fn powder_pawley_width_info(
    text: &str,
    pattern_text: &str,
    profile: &str,
    refine_size: bool,
    refine_strain: bool,
) -> Result<JsValue, JsValue> {
    pawley_info(text, pattern_text, profile, refine_size, refine_strain)
}

fn pawley_info(
    text: &str,
    pattern_text: &str,
    profile: &str,
    refine_size: bool,
    refine_strain: bool,
) -> Result<JsValue, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let observed = pattern(pattern_text)?;
    let kind = match profile.to_ascii_lowercase().as_str() {
        "gaussian" => xtalforge::analysis::ProfileKind::Gaussian,
        "lorentzian" => xtalforge::analysis::ProfileKind::Lorentzian,
        "pseudo_voigt" | "pseudovoigt" | "pseudo-voigt" => {
            xtalforge::analysis::ProfileKind::PseudoVoigt
        }
        other => {
            return Err(crate::common::invalid_input(&format!(
                "unknown profile {other:?}; use gaussian, lorentzian or pseudo_voigt"
            )));
        }
    };
    let result = xtalforge::analysis::pawley_refine(
        &structure,
        &observed,
        &xtalforge::analysis::PawleyOptions {
            kind,
            refine_size,
            refine_strain,
            ..xtalforge::analysis::PawleyOptions::default()
        },
    )
    .map_err(to_js_error)?;
    let object = js_sys::Object::new();
    let reflections = js_sys::Array::new();
    for reflection in &result.reflections {
        let item = js_sys::Object::new();
        set(&item, "hkl", array_i32(reflection.hkl).into());
        set(&item, "d", JsValue::from_f64(reflection.d));
        set(&item, "two_theta", JsValue::from_f64(reflection.two_theta));
        set(&item, "intensity", JsValue::from_f64(reflection.intensity));
        set(&item, "esd", JsValue::from_f64(reflection.esd));
        reflections.push(&item.into());
    }
    set(&object, "reflections", reflections.into());
    let cell = js_sys::Array::new();
    for value in result.cell {
        cell.push(&JsValue::from_f64(value));
    }
    set(&object, "cell", cell.into());
    for (key, value) in [
        ("zero_shift", result.zero_shift),
        ("zero_shift_esd", result.zero_shift_esd),
        ("rp", result.rp),
        ("rwp", result.rwp),
        ("rexp", result.rexp),
        ("chi_squared", result.chi_squared),
        ("size", result.size),
        ("strain", result.strain),
    ] {
        set(&object, key, JsValue::from_f64(value));
    }
    set(&object, "converged", JsValue::from_bool(result.converged));
    // The fitted profile at the observed points, for an observed/calculated
    // overlay.
    let calculated = js_sys::Array::new();
    for point in &result.pattern {
        calculated.push(&JsValue::from_f64(point.intensity));
    }
    set(&object, "calculated", calculated.into());
    Ok(object.into())
}
