//! Symmetry analysis exports (Task 2.8-2.12).

use wasm_bindgen::prelude::*;

use crate::common::{array3, parse_any, set, summary_object, to_js_error};

/// Standardizes a structure to the conventional cell.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or a failed search.
#[wasm_bindgen]
pub fn standardize_info(text: &str, symprec: f64) -> Result<JsValue, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let result = xtalforge::symmetry::standardize(
        &structure,
        xtalforge::symmetry::SymmetryOptions {
            symprec,
            ..Default::default()
        },
    )
    .map_err(to_js_error)?;
    summary_object(&result.structure)
}

/// Builds the standardized primitive cell.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or a failed search.
#[wasm_bindgen]
pub fn primitive_info(text: &str, symprec: f64) -> Result<JsValue, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let result = xtalforge::symmetry::primitive_cell(
        &structure,
        xtalforge::symmetry::SymmetryOptions {
            symprec,
            ..Default::default()
        },
    )
    .map_err(to_js_error)?;
    summary_object(&result.structure)
}

/// Reduces a structure with Niggli, Delaunay or Minkowski reduction.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or an unknown algorithm.
#[wasm_bindgen]
pub fn reduce_info(text: &str, algorithm: &str) -> Result<JsValue, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let reduction = match algorithm {
        "niggli" => xtalforge::symmetry::Reduction::Niggli,
        "delaunay" => xtalforge::symmetry::Reduction::Delaunay,
        "minkowski" => xtalforge::symmetry::Reduction::Minkowski,
        other => {
            return Err(JsValue::from_str(&format!(
                "invalid_input: unknown reduction {other:?}"
            )));
        }
    };
    let result =
        xtalforge::symmetry::reduce_structure(&structure, reduction).map_err(to_js_error)?;
    let object = summary_object(&result.structure)?;
    let object: js_sys::Object = object.unchecked_into();
    set(
        &object,
        "determinant",
        JsValue::from_f64(result.determinant),
    );
    Ok(object.into())
}

/// Returns Wyckoff assignments.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or a failed search.
#[wasm_bindgen]
pub fn wyckoff_info(text: &str, symprec: f64) -> Result<js_sys::Array, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let assignments = xtalforge::symmetry::wyckoff_assignments(
        &structure,
        xtalforge::symmetry::SymmetryOptions {
            symprec,
            ..Default::default()
        },
    )
    .map_err(to_js_error)?;
    let array = js_sys::Array::new();
    for assignment in assignments {
        let object = js_sys::Object::new();
        set(&object, "atom", JsValue::from_f64(assignment.atom as f64));
        set(&object, "element", JsValue::from_str(&assignment.element));
        set(&object, "label", JsValue::from_str(&assignment.label));
        set(
            &object,
            "letter",
            JsValue::from_str(&assignment.letter.to_string()),
        );
        set(
            &object,
            "multiplicity",
            JsValue::from_f64(assignment.multiplicity as f64),
        );
        set(
            &object,
            "site_symmetry",
            JsValue::from_str(&assignment.site_symmetry),
        );
        set(
            &object,
            "coordinates",
            JsValue::from_str(&assignment.coordinates),
        );
        array.push(&object.into());
    }
    Ok(array)
}

/// Returns the Wigner-Seitz Brillouin zone.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or a degenerate lattice.
#[wasm_bindgen]
pub fn brillouin_zone(text: &str) -> Result<JsValue, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let zone = xtalforge::symmetry::brillouin_zone(structure.lattice()).map_err(to_js_error)?;
    let object = js_sys::Object::new();
    let vertices = js_sys::Array::new();
    for vertex in &zone.vertices {
        vertices.push(&array3(*vertex).into());
    }
    set(&object, "vertices", vertices.into());
    let faces = js_sys::Array::new();
    for face in &zone.faces {
        let indices = js_sys::Array::new();
        for index in face {
            indices.push(&JsValue::from_f64(*index as f64));
        }
        faces.push(&indices.into());
    }
    set(&object, "faces", faces.into());
    let edges = js_sys::Array::new();
    for (first, second) in &zone.edges {
        let edge = js_sys::Array::new();
        edge.push(&JsValue::from_f64(*first as f64));
        edge.push(&JsValue::from_f64(*second as f64));
        edges.push(&edge.into());
    }
    set(&object, "edges", edges.into());
    Ok(object.into())
}

/// Detects across the default tolerance ladder.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or a failed search.
#[wasm_bindgen]
pub fn tolerance_ladder_info(text: &str) -> Result<js_sys::Array, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let ladder =
        xtalforge::symmetry::tolerance_ladder(&structure, &xtalforge::symmetry::DEFAULT_LADDER)
            .map_err(to_js_error)?;
    let array = js_sys::Array::new();
    for range in ladder {
        let object = js_sys::Object::new();
        set(&object, "min", JsValue::from_f64(range.min));
        set(&object, "max", JsValue::from_f64(range.max));
        set(
            &object,
            "number",
            JsValue::from_f64(f64::from(range.number)),
        );
        set(
            &object,
            "hall_number",
            JsValue::from_f64(f64::from(range.hall_number)),
        );
        set(&object, "hm_symbol", JsValue::from_str(&range.hm_symbol));
        array.push(&object.into());
    }
    Ok(array)
}

/// Reports per-atom deviations from a space group.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or an unknown symbol.
#[wasm_bindgen]
pub fn diagnose_info(text: &str, symbol: &str) -> Result<JsValue, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let group = xtalforge::symmetry::SpaceGroup::parse(symbol).map_err(to_js_error)?;
    let diagnostics = xtalforge::symmetry::diagnose(&structure, &group).map_err(to_js_error)?;
    let object = js_sys::Object::new();
    set(
        &object,
        "max_displacement",
        JsValue::from_f64(diagnostics.max_displacement),
    );
    set(
        &object,
        "mean_displacement",
        JsValue::from_f64(diagnostics.mean_displacement),
    );
    let atoms = js_sys::Array::new();
    for atom in diagnostics.atoms {
        let item = js_sys::Object::new();
        set(&item, "atom", JsValue::from_f64(atom.atom as f64));
        set(&item, "element", JsValue::from_str(&atom.element));
        set(&item, "displacement", JsValue::from_f64(atom.displacement));
        atoms.push(&item.into());
    }
    set(&object, "atoms", atoms.into());
    Ok(object.into())
}

/// Detects the layer group of a slab.
///
/// # Errors
///
/// Returns a JavaScript `Error` when the layer contract is violated or the
/// search fails.
#[wasm_bindgen]
pub fn layer_group_info(text: &str, symprec: f64) -> Result<JsValue, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let dataset =
        xtalforge::symmetry::detect_layer_symmetry(&structure, symprec).map_err(to_js_error)?;
    let object = js_sys::Object::new();
    set(
        &object,
        "number",
        JsValue::from_f64(f64::from(dataset.number)),
    );
    set(
        &object,
        "hall_number",
        JsValue::from_f64(f64::from(dataset.hall_number)),
    );
    set(&object, "hm_symbol", JsValue::from_str(&dataset.hm_symbol));
    set(
        &object,
        "pearson_symbol",
        JsValue::from_str(&dataset.pearson_symbol),
    );
    let operations = js_sys::Array::new();
    for operation in &dataset.operations {
        operations.push(&JsValue::from_str(&operation.triplet()));
    }
    set(&object, "operations", operations.into());
    let wyckoffs = js_sys::Array::new();
    for wyckoff in &dataset.wyckoffs {
        wyckoffs.push(&JsValue::from_str(&wyckoff.to_string()));
    }
    set(&object, "wyckoffs", wyckoffs.into());
    Ok(object.into())
}

/// Detects the collinear magnetic space group (experimental).
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or a failed search.
#[wasm_bindgen]
pub fn magnetic_group_info(text: &str, symprec: f64) -> Result<JsValue, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let dataset =
        xtalforge::symmetry::detect_magnetic_symmetry(&structure, symprec).map_err(to_js_error)?;
    let object = js_sys::Object::new();
    set(
        &object,
        "uni_number",
        JsValue::from_f64(f64::from(dataset.uni_number)),
    );
    set(
        &object,
        "operations",
        JsValue::from_f64(dataset.operations.len() as f64),
    );
    Ok(object.into())
}

/// Changes the setting of a structure to a Hall number.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or an unknown Hall number.
#[wasm_bindgen]
pub fn change_setting_info(text: &str, hall_number: u32) -> Result<JsValue, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let changed =
        xtalforge::symmetry::change_setting(&structure, hall_number).map_err(to_js_error)?;
    summary_object(&changed.structure)
}

/// Returns the Bravais type and band path of a structure.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or a k-path failure.
#[wasm_bindgen]
pub fn kpath_info(text: &str) -> Result<JsValue, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let path = xtalforge::symmetry::band_path_for_structure(
        &structure,
        xtalforge::symmetry::SymmetryOptions {
            symprec: 1e-3,
            ..Default::default()
        },
    )
    .map_err(to_js_error)?;
    let vertices = path.vertices(structure.lattice()).map_err(to_js_error)?;
    let object = js_sys::Object::new();
    set(&object, "bravais", JsValue::from_str(&path.bravais));
    let points = js_sys::Object::new();
    for (label, point) in &path.special_points {
        set(&points, label, array3(*point).into());
    }
    set(&object, "special_points", points.into());
    let segments = js_sys::Array::new();
    for (start, end) in &path.path {
        let pair = js_sys::Array::new();
        pair.push(&JsValue::from_str(start));
        pair.push(&JsValue::from_str(end));
        segments.push(&pair.into());
    }
    set(&object, "path", segments.into());
    let list = js_sys::Array::new();
    for (label, point, distance) in vertices {
        let item = js_sys::Object::new();
        set(&item, "label", JsValue::from_str(&label));
        set(&item, "cartesian", array3(point).into());
        set(&item, "distance", JsValue::from_f64(distance));
        list.push(&item.into());
    }
    set(&object, "vertices", list.into());
    Ok(object.into())
}

/// Detects the molecular point group of a finite structure.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or an unsupported group.
#[wasm_bindgen]
pub fn pointgroup_info(text: &str, tolerance: f64) -> Result<JsValue, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let p1 = structure.p1().map_err(to_js_error)?;
    let points: Vec<(String, [f64; 3])> = p1
        .atoms
        .iter()
        .map(|atom| (atom.element.clone(), p1.lattice.to_cart(atom.frac)))
        .collect();
    let group =
        xtalforge::symmetry::molecular_point_group(&points, tolerance).map_err(to_js_error)?;
    let object = js_sys::Object::new();
    set(
        &object,
        "schoenflies",
        JsValue::from_str(&group.schoenflies),
    );
    set(&object, "order", JsValue::from_f64(group.order as f64));
    set(
        &object,
        "proper_rotations",
        JsValue::from_f64(group.proper_rotations as f64),
    );
    set(
        &object,
        "improper_operations",
        JsValue::from_f64(group.improper_operations as f64),
    );
    set(
        &object,
        "has_inversion",
        JsValue::from_bool(group.has_inversion),
    );
    set(&object, "mirrors", JsValue::from_f64(group.mirrors as f64));
    Ok(object.into())
}
