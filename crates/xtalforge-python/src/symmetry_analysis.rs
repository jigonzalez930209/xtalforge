//! Symmetry analysis and workflow bindings (Task 2.5-2.12).

use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

use xtalforge::symmetry::SymmetryOptions;

use crate::structures::load_structure;

/// Standardizes a structure to the conventional cell.
#[pyfunction]
#[pyo3(signature = (text, symprec = 1e-3, format_hint = None))]
pub(crate) fn standardize<'py>(
    py: Python<'py>,
    text: &str,
    symprec: f64,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyDict>> {
    let structure = load_structure(text, format_hint)?;
    let result = xtalforge::symmetry::standardize(
        &structure,
        SymmetryOptions {
            symprec,
            ..SymmetryOptions::default()
        },
    )
    .map_err(crate::map_error)?;
    let dict = crate::structures::structure_payload(py, &result.structure)?;
    dict.set_item("origin_shift", result.origin_shift.to_vec())?;
    dict.set_item(
        "rotation_matrix",
        result
            .rotation_matrix
            .iter()
            .map(|row| row.to_vec())
            .collect::<Vec<_>>(),
    )?;
    Ok(dict)
}

/// Builds the standardized primitive cell.
#[pyfunction]
#[pyo3(signature = (text, symprec = 1e-3, format_hint = None))]
pub(crate) fn primitive_cell<'py>(
    py: Python<'py>,
    text: &str,
    symprec: f64,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyDict>> {
    let structure = load_structure(text, format_hint)?;
    let result = xtalforge::symmetry::primitive_cell(
        &structure,
        SymmetryOptions {
            symprec,
            ..SymmetryOptions::default()
        },
    )
    .map_err(crate::map_error)?;
    crate::structures::structure_payload(py, &result.structure)
}

/// Reduces a structure to a Niggli, Delaunay or Minkowski cell.
#[pyfunction]
#[pyo3(signature = (text, algorithm = "niggli".to_string(), format_hint = None))]
pub(crate) fn reduce_structure<'py>(
    py: Python<'py>,
    text: &str,
    algorithm: String,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyDict>> {
    let structure = load_structure(text, format_hint)?;
    let reduction = match algorithm.as_str() {
        "niggli" => xtalforge::symmetry::Reduction::Niggli,
        "delaunay" => xtalforge::symmetry::Reduction::Delaunay,
        "minkowski" => xtalforge::symmetry::Reduction::Minkowski,
        other => {
            return Err(crate::invalid_input(format!(
                "unknown reduction {other:?}; use niggli, delaunay or minkowski"
            )));
        }
    };
    let result =
        xtalforge::symmetry::reduce_structure(&structure, reduction).map_err(crate::map_error)?;
    let dict = crate::structures::structure_payload(py, &result.structure)?;
    dict.set_item(
        "transformation",
        result
            .transformation
            .iter()
            .map(|row| row.to_vec())
            .collect::<Vec<_>>(),
    )?;
    dict.set_item("determinant", result.determinant)?;
    Ok(dict)
}

/// Assigns Wyckoff letters, multiplicities and free parameters.
#[pyfunction]
#[pyo3(signature = (text, symprec = 1e-3, format_hint = None))]
pub(crate) fn wyckoff_assignments<'py>(
    py: Python<'py>,
    text: &str,
    symprec: f64,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyList>> {
    let structure = load_structure(text, format_hint)?;
    let assignments = xtalforge::symmetry::wyckoff_assignments(
        &structure,
        SymmetryOptions {
            symprec,
            ..SymmetryOptions::default()
        },
    )
    .map_err(crate::map_error)?;
    let items = PyList::empty(py);
    for assignment in assignments {
        let dict = PyDict::new(py);
        dict.set_item("atom", assignment.atom)?;
        dict.set_item("element", assignment.element)?;
        dict.set_item("label", assignment.label)?;
        dict.set_item("letter", assignment.letter.to_string())?;
        dict.set_item("site_symmetry", assignment.site_symmetry)?;
        dict.set_item("multiplicity", assignment.multiplicity)?;
        dict.set_item("coordinates", assignment.coordinates)?;
        dict.set_item(
            "free_parameters",
            assignment
                .free_parameters
                .iter()
                .map(|(name, value)| format!("{name}={value:.6}"))
                .collect::<Vec<_>>(),
        )?;
        items.append(dict)?;
    }
    Ok(items)
}

/// Returns the Wigner-Seitz Brillouin zone.
#[pyfunction]
#[pyo3(signature = (text, format_hint = None))]
pub(crate) fn brillouin_zone<'py>(
    py: Python<'py>,
    text: &str,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyDict>> {
    let structure = load_structure(text, format_hint)?;
    let zone =
        xtalforge::symmetry::brillouin_zone(structure.lattice()).map_err(crate::map_error)?;
    let dict = PyDict::new(py);
    dict.set_item(
        "vertices",
        zone.vertices
            .iter()
            .map(|point| point.to_vec())
            .collect::<Vec<_>>(),
    )?;
    dict.set_item("faces", zone.faces)?;
    dict.set_item(
        "edges",
        zone.edges
            .iter()
            .map(|(first, second)| vec![*first, *second])
            .collect::<Vec<_>>(),
    )?;
    Ok(dict)
}

/// Detects across the default tolerance ladder.
#[pyfunction]
#[pyo3(signature = (text, format_hint = None))]
pub(crate) fn tolerance_ladder<'py>(
    py: Python<'py>,
    text: &str,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyList>> {
    let structure = load_structure(text, format_hint)?;
    let ladder =
        xtalforge::symmetry::tolerance_ladder(&structure, &xtalforge::symmetry::DEFAULT_LADDER)
            .map_err(crate::map_error)?;
    let items = PyList::empty(py);
    for range in ladder {
        let dict = PyDict::new(py);
        dict.set_item("min", range.min)?;
        dict.set_item("max", range.max)?;
        dict.set_item("number", range.number)?;
        dict.set_item("hall_number", range.hall_number)?;
        dict.set_item("hm_symbol", range.hm_symbol)?;
        items.append(dict)?;
    }
    Ok(items)
}

/// Reports per-atom deviations from a space group.
#[pyfunction]
#[pyo3(signature = (text, symbol, format_hint = None))]
pub(crate) fn diagnose<'py>(
    py: Python<'py>,
    text: &str,
    symbol: &str,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyDict>> {
    let structure = load_structure(text, format_hint)?;
    let group = xtalforge::symmetry::SpaceGroup::parse(symbol).map_err(crate::map_error)?;
    let diagnostics =
        xtalforge::symmetry::diagnose(&structure, &group).map_err(crate::map_error)?;
    let dict = PyDict::new(py);
    dict.set_item("max_displacement", diagnostics.max_displacement)?;
    dict.set_item("mean_displacement", diagnostics.mean_displacement)?;
    let atoms = PyList::empty(py);
    for atom in diagnostics.atoms {
        let item = PyDict::new(py);
        item.set_item("atom", atom.atom)?;
        item.set_item("element", atom.element)?;
        item.set_item("label", atom.label)?;
        item.set_item("displacement", atom.displacement)?;
        atoms.append(item)?;
    }
    dict.set_item("atoms", atoms)?;
    Ok(dict)
}

/// Detects the layer group of a slab.
#[pyfunction]
#[pyo3(signature = (text, symprec = 1e-3, format_hint = None))]
pub(crate) fn detect_layer_symmetry<'py>(
    py: Python<'py>,
    text: &str,
    symprec: f64,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyDict>> {
    let structure = load_structure(text, format_hint)?;
    let dataset = xtalforge::symmetry::detect_layer_symmetry(&structure, symprec)
        .map_err(crate::map_error)?;
    let dict = PyDict::new(py);
    dict.set_item("number", dataset.number)?;
    dict.set_item("hall_number", dataset.hall_number)?;
    dict.set_item("hm_symbol", dataset.hm_symbol)?;
    dict.set_item("pearson_symbol", dataset.pearson_symbol)?;
    dict.set_item("orbits", dataset.orbits)?;
    dict.set_item(
        "wyckoffs",
        dataset
            .wyckoffs
            .iter()
            .map(char::to_string)
            .collect::<Vec<_>>(),
    )?;
    dict.set_item("site_symmetry_symbols", dataset.site_symmetry_symbols)?;
    dict.set_item(
        "operations",
        dataset
            .operations
            .iter()
            .map(|operation| operation.triplet())
            .collect::<Vec<_>>(),
    )?;
    Ok(dict)
}

/// Detects the collinear magnetic space group (experimental).
#[pyfunction]
#[pyo3(signature = (text, symprec = 1e-3, format_hint = None))]
pub(crate) fn detect_magnetic_symmetry<'py>(
    py: Python<'py>,
    text: &str,
    symprec: f64,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyDict>> {
    let structure = load_structure(text, format_hint)?;
    let dataset = xtalforge::symmetry::detect_magnetic_symmetry(&structure, symprec)
        .map_err(crate::map_error)?;
    let dict = PyDict::new(py);
    dict.set_item("uni_number", dataset.uni_number)?;
    dict.set_item("orbits", dataset.orbits)?;
    let operations = PyList::empty(py);
    for operation in dataset.operations {
        let item = PyDict::new(py);
        item.set_item(
            "rotation",
            operation
                .rotation
                .iter()
                .map(|row| row.to_vec())
                .collect::<Vec<_>>(),
        )?;
        item.set_item("translation", operation.translation.to_vec())?;
        item.set_item("time_reversal", operation.time_reversal)?;
        operations.append(item)?;
    }
    dict.set_item("operations", operations)?;
    Ok(dict)
}

/// Changes the setting of a structure to a Hall number.
#[pyfunction]
#[pyo3(signature = (text, hall_number, format_hint = None))]
pub(crate) fn change_setting<'py>(
    py: Python<'py>,
    text: &str,
    hall_number: u32,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyDict>> {
    let structure = load_structure(text, format_hint)?;
    let changed =
        xtalforge::symmetry::change_setting(&structure, hall_number).map_err(crate::map_error)?;
    crate::structures::structure_payload(py, &changed.structure)
}

/// Returns the Bravais type, special points and band path of a structure.
#[pyfunction]
#[pyo3(signature = (text, format_hint = None))]
pub(crate) fn band_path<'py>(
    py: Python<'py>,
    text: &str,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyDict>> {
    let structure = load_structure(text, format_hint)?;
    let path = xtalforge::symmetry::band_path_for_structure(
        &structure,
        SymmetryOptions {
            symprec: 1e-3,
            ..SymmetryOptions::default()
        },
    )
    .map_err(crate::map_error)?;
    let vertices = path
        .vertices(structure.lattice())
        .map_err(crate::map_error)?;
    let dict = PyDict::new(py);
    dict.set_item("bravais", &path.bravais)?;
    dict.set_item(
        "special_points",
        path.special_points
            .iter()
            .map(|(label, point)| (label.clone(), point.to_vec()))
            .collect::<std::collections::BTreeMap<_, _>>(),
    )?;
    dict.set_item("path", &path.path)?;
    dict.set_item(
        "vertices",
        vertices
            .iter()
            .map(|(label, point, distance)| (label.clone(), point.to_vec(), *distance))
            .collect::<Vec<_>>(),
    )?;
    Ok(dict)
}

/// Detects the molecular point group of a finite structure.
#[pyfunction]
#[pyo3(signature = (text, tolerance = 1e-3, format_hint = None))]
pub(crate) fn molecular_point_group<'py>(
    py: Python<'py>,
    text: &str,
    tolerance: f64,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyDict>> {
    let structure = load_structure(text, format_hint)?;
    let p1 = structure.p1().map_err(crate::map_error)?;
    let points: Vec<(String, [f64; 3])> = p1
        .atoms
        .iter()
        .map(|atom| (atom.element.clone(), p1.lattice.to_cart(atom.frac)))
        .collect();
    let group =
        xtalforge::symmetry::molecular_point_group(&points, tolerance).map_err(crate::map_error)?;
    let dict = PyDict::new(py);
    dict.set_item("schoenflies", group.schoenflies)?;
    dict.set_item("order", group.order)?;
    dict.set_item("proper_rotations", group.proper_rotations)?;
    dict.set_item("improper_operations", group.improper_operations)?;
    dict.set_item("has_inversion", group.has_inversion)?;
    dict.set_item("mirrors", group.mirrors)?;
    Ok(dict)
}
