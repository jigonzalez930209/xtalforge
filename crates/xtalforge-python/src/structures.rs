//! Structure reading, conversion and construction bindings (Task 1.17.2).

use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

use xtalforge::core::{Lattice, P1Cell, Site, Structure, StructureBuilder};

/// Reads a structure from text with a format hint.
///
/// Parameters
/// ----------
/// text:
///     File contents.
/// format_hint:
///     Format name or extension: ``cif``, ``poscar`` or ``xyz``.
///
/// Returns
/// -------
/// dict
///     Summary (``name``, ``formula``, ``hall``, ``n_sites``, ``volume``,
///     ``density``, ``total_charge``, ``is_p1``, ``n_atoms``), the ``lattice``
///     parameters and matrix, the asymmetric ``sites``, the expanded ``p1``
///     cell and the versioned structure JSON under ``json``.
#[pyfunction]
#[pyo3(signature = (text, format_hint))]
pub(crate) fn read_structure<'py>(
    py: Python<'py>,
    text: &str,
    format_hint: &str,
) -> PyResult<Bound<'py, PyDict>> {
    let structure = load_structure(text, Some(format_hint))?;
    structure_payload(py, &structure)
}

/// Reads CIF text and returns a structure dictionary.
///
/// The dictionary has the same shape as :func:`read_structure` with the
/// ``cif`` format hint.
#[pyfunction]
pub(crate) fn read_cif<'py>(py: Python<'py>, text: &str) -> PyResult<Bound<'py, PyDict>> {
    let text = strip_bom(text);
    let structure = xtalforge::io::read_cif(text).map_err(crate::map_error)?;
    structure_payload(py, &structure)
}

/// Builds ``na x nb x nc`` repeats and returns the P1 supercell dictionary.
///
/// The input is expanded to P1 before replication, so every atom of the
/// crystal appears in the result. The returned dictionary also carries
/// ``repeat``.
#[pyfunction]
#[pyo3(signature = (text, na, nb, nc, format_hint = None))]
pub(crate) fn supercell<'py>(
    py: Python<'py>,
    text: &str,
    na: i64,
    nb: i64,
    nc: i64,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyDict>> {
    let structure = load_structure(text, format_hint)?;
    let p1 = as_p1_structure(&structure)?;
    let cell = p1
        .supercell(repeat(na)?, repeat(nb)?, repeat(nc)?)
        .map_err(crate::map_error)?;
    let dict = structure_payload(py, &cell)?;
    dict.set_item("repeat", [na, nb, nc])?;
    Ok(dict)
}

/// Cuts a slab along the ``(h k l)`` plane.
///
/// ``layers`` is the number of interplanar spacings kept and ``vacuum`` the
/// empty space above the slab in Angstrom. The returned dictionary also
/// carries ``hkl``, ``layers``, ``thickness``, ``vacuum`` and ``cut``.
#[pyfunction]
#[pyo3(signature = (text, h, k, l, layers, vacuum, format_hint = None))]
#[allow(clippy::too_many_arguments)] // one Python argument per Miller index and slab option
pub(crate) fn slab<'py>(
    py: Python<'py>,
    text: &str,
    h: i32,
    k: i32,
    l: i32,
    layers: i64,
    vacuum: f64,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyDict>> {
    let structure = load_structure(text, format_hint)?;
    let slab = xtalforge::core::make_slab(&structure, [h, k, l], repeat(layers)?, vacuum, 0.0)
        .map_err(crate::map_error)?;
    let dict = structure_payload(py, &slab.structure)?;
    dict.set_item("hkl", slab.hkl.to_vec())?;
    dict.set_item("layers", slab.layers)?;
    dict.set_item("thickness", slab.thickness)?;
    dict.set_item("vacuum", slab.vacuum)?;
    dict.set_item("cut", slab.cut)?;
    Ok(dict)
}

pub(crate) fn load_structure(text: &str, format_hint: Option<&str>) -> PyResult<Structure> {
    let text = strip_bom(text);
    if text.trim_start().starts_with('{') {
        let structure = if text.contains("\"schema_version\"") {
            xtalforge::core::StructureDocument::from_json(text).map(|document| document.structure)
        } else {
            serde_json::from_str(text).map_err(|error| {
                xtalforge::core::Error::InvalidInput(format!("invalid structure JSON: {error}"))
            })
        };
        return structure.map_err(crate::map_error);
    }
    let hint = format_hint.unwrap_or_else(|| detect_format(text));
    xtalforge::io::read_structure(text, hint).map_err(crate::map_error)
}

fn strip_bom(text: &str) -> &str {
    text.strip_prefix('\u{feff}').unwrap_or(text)
}

/// Detects the format from the text: CIF tags, an XYZ atom count or POSCAR.
pub(crate) fn detect_format(text: &str) -> &'static str {
    if text.trim_start().starts_with('{') {
        return "json";
    }
    if text.trim_start().starts_with("data_")
        || text.contains("_cell_length_a")
        || text.contains("_atom_site_")
    {
        return "cif";
    }
    let mut lines = text.lines().filter(|line| !line.trim().is_empty());
    let first = lines.next().unwrap_or("").trim();
    let second = lines.next().unwrap_or("").trim();
    if first.parse::<usize>().is_ok() && second.parse::<f64>().is_err() {
        "xyz"
    } else {
        "poscar"
    }
}

fn repeat(value: i64) -> PyResult<u32> {
    u32::try_from(value).map_err(|_| {
        crate::invalid_input(format!(
            "repeat count must be a positive integer that fits in 32 bits, got {value}"
        ))
    })
}

/// Turns a structure into an explicit P1 structure without symmetry
/// operations, so a supercell contains every atom of the crystal.
pub(crate) fn as_p1_structure(structure: &Structure) -> PyResult<Structure> {
    if structure.is_p1() {
        return Ok(structure.clone());
    }
    let p1 = structure.p1().map_err(crate::map_error)?;
    let mut builder = StructureBuilder::new(structure.name(), *structure.lattice());
    for atom in &p1.atoms {
        let site = Site::new(atom.label.clone(), &atom.element, atom.frac)
            .and_then(|site| site.with_occupancy(atom.occupancy))
            .map_err(crate::map_error)?;
        builder = builder.site(site);
    }
    for (key, value) in structure.metadata() {
        builder = builder.metadata(key.clone(), value.clone());
    }
    builder.build().map_err(crate::map_error)
}

pub(crate) fn structure_payload<'py>(
    py: Python<'py>,
    structure: &Structure,
) -> PyResult<Bound<'py, PyDict>> {
    let summary = structure.summary().map_err(crate::map_error)?;
    let p1 = structure.p1().map_err(crate::map_error)?;
    let dict = PyDict::new(py);
    dict.set_item("name", summary.name)?;
    dict.set_item("formula", summary.formula)?;
    dict.set_item("hall", summary.hall)?;
    dict.set_item("n_sites", summary.n_sites)?;
    dict.set_item("n_atoms", p1.n_atoms())?;
    dict.set_item("volume", summary.volume)?;
    dict.set_item("density", summary.density)?;
    dict.set_item("total_charge", summary.total_charge)?;
    dict.set_item("is_p1", summary.is_p1)?;
    dict.set_item("lattice", lattice_payload(py, structure.lattice())?)?;
    dict.set_item("sites", sites_payload(py, structure)?)?;
    dict.set_item("p1", p1_payload(py, &p1)?)?;
    dict.set_item(
        "json",
        structure
            .to_document()
            .to_json()
            .map_err(crate::map_error)?,
    )?;
    Ok(dict)
}

fn lattice_payload<'py>(py: Python<'py>, lattice: &Lattice) -> PyResult<Bound<'py, PyDict>> {
    let [a, b, c, alpha, beta, gamma] = lattice.parameters();
    let dict = PyDict::new(py);
    dict.set_item("parameters", [a, b, c, alpha, beta, gamma].to_vec())?;
    dict.set_item("a", a)?;
    dict.set_item("b", b)?;
    dict.set_item("c", c)?;
    dict.set_item("alpha", alpha)?;
    dict.set_item("beta", beta)?;
    dict.set_item("gamma", gamma)?;
    dict.set_item("volume", lattice.volume())?;
    dict.set_item(
        "matrix",
        lattice
            .matrix()
            .iter()
            .map(|row| row.to_vec())
            .collect::<Vec<_>>(),
    )?;
    Ok(dict)
}

fn sites_payload<'py>(py: Python<'py>, structure: &Structure) -> PyResult<Bound<'py, PyList>> {
    let mut sites = Vec::with_capacity(structure.sites().len());
    for site in structure.sites() {
        let [x, y, z] = site.frac();
        let dict = PyDict::new(py);
        dict.set_item("label", site.label())?;
        dict.set_item("element", site.element())?;
        dict.set_item("frac", [x, y, z].to_vec())?;
        dict.set_item("cart", structure.lattice().to_cart([x, y, z]).to_vec())?;
        dict.set_item("occupancy", site.occupancy())?;
        dict.set_item("u_iso", site.u_iso())?;
        dict.set_item("charge", site.charge())?;
        dict.set_item("disorder_group", site.disorder_group())?;
        sites.push(dict);
    }
    PyList::new(py, sites)
}

fn p1_payload<'py>(py: Python<'py>, p1: &P1Cell) -> PyResult<Bound<'py, PyDict>> {
    let frac = p1.frac();
    let cart = p1.lattice.to_cart_many(&frac);
    let dict = PyDict::new(py);
    dict.set_item("n_atoms", p1.n_atoms())?;
    dict.set_item("elements", p1.elements())?;
    dict.set_item(
        "labels",
        p1.atoms
            .iter()
            .map(|atom| atom.label.clone())
            .collect::<Vec<_>>(),
    )?;
    dict.set_item(
        "frac",
        frac.iter().map(|point| point.to_vec()).collect::<Vec<_>>(),
    )?;
    dict.set_item(
        "cart",
        cart.iter().map(|point| point.to_vec()).collect::<Vec<_>>(),
    )?;
    dict.set_item("occupancy", p1.occupancy())?;
    dict.set_item(
        "site_idx",
        p1.atoms
            .iter()
            .map(|atom| atom.site_idx)
            .collect::<Vec<_>>(),
    )?;
    dict.set_item(
        "op_idx",
        p1.atoms.iter().map(|atom| atom.op_idx).collect::<Vec<_>>(),
    )?;
    Ok(dict)
}

/// Builds a graphene sheet or stack.
#[pyfunction]
#[pyo3(signature = (nx, ny, layers = 1, interlayer = 3.35, vacuum = 15.0))]
pub(crate) fn build_graphene<'py>(
    py: Python<'py>,
    nx: u32,
    ny: u32,
    layers: u32,
    interlayer: f64,
    vacuum: f64,
) -> PyResult<Bound<'py, PyDict>> {
    let structure =
        xtalforge::build::graphene(nx, ny, layers, interlayer, vacuum).map_err(crate::map_error)?;
    structure_payload(py, &structure)
}

/// Builds a chiral nanotube.
#[pyfunction]
#[pyo3(signature = (n, m, length = 1, vacuum = 10.0))]
pub(crate) fn build_nanotube<'py>(
    py: Python<'py>,
    n: i32,
    m: i32,
    length: u32,
    vacuum: f64,
) -> PyResult<Bound<'py, PyDict>> {
    let structure = xtalforge::build::nanotube(n, m, length, vacuum).map_err(crate::map_error)?;
    structure_payload(py, &structure)
}

/// Builds C60.
#[pyfunction]
pub(crate) fn build_c60<'py>(py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
    let structure = xtalforge::build::fullerene(60).map_err(crate::map_error)?;
    structure_payload(py, &structure)
}
