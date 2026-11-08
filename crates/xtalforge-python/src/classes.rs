//! Python classes for the core model (Task 1.17.2).

use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

use xtalforge::core::{
    Lattice as CoreLattice, P1Atom, P1Cell as CoreP1Cell, Site as CoreSite,
    Structure as CoreStructure,
};

fn numpy_array1(py: Python<'_>, values: Vec<f64>) -> PyResult<Py<PyAny>> {
    Ok(py
        .import("numpy")?
        .call_method1("array", (values,))?
        .unbind())
}

fn numpy_array2(py: Python<'_>, values: Vec<Vec<f64>>) -> PyResult<Py<PyAny>> {
    Ok(py
        .import("numpy")?
        .call_method1("array", (values,))?
        .unbind())
}

fn frac(values: &[f64]) -> PyResult<[f64; 3]> {
    let [x, y, z] = values else {
        return Err(crate::invalid_input(
            "expected three fractional coordinates",
        ));
    };
    Ok([*x, *y, *z])
}

/// A unit cell whose matrix rows are the vectors a, b and c.
#[pyclass(name = "Lattice", module = "xtalforge")]
#[derive(Clone)]
pub(crate) struct Lattice {
    inner: CoreLattice,
}

#[pymethods]
impl Lattice {
    #[new]
    #[pyo3(signature = (a, b, c, alpha, beta, gamma))]
    fn new(a: f64, b: f64, c: f64, alpha: f64, beta: f64, gamma: f64) -> PyResult<Self> {
        Ok(Self {
            inner: CoreLattice::from_parameters(a, b, c, alpha, beta, gamma)
                .map_err(crate::map_error)?,
        })
    }

    /// Builds a cell from a 3x3 matrix whose rows are a, b and c.
    #[staticmethod]
    fn from_matrix(matrix: Vec<Vec<f64>>) -> PyResult<Self> {
        if matrix.len() != 3 {
            return Err(crate::invalid_input("matrix must have three rows"));
        }
        let mut rows = [[0.0; 3]; 3];
        for (slot, row) in rows.iter_mut().zip(matrix.iter()) {
            *slot = frac(row)?;
        }
        Ok(Self {
            inner: CoreLattice::new(rows).map_err(crate::map_error)?,
        })
    }

    #[getter]
    fn a(&self) -> f64 {
        self.inner.lengths()[0]
    }

    #[getter]
    fn b(&self) -> f64 {
        self.inner.lengths()[1]
    }

    #[getter]
    fn c(&self) -> f64 {
        self.inner.lengths()[2]
    }

    #[getter]
    fn alpha(&self) -> f64 {
        self.inner.angles()[0]
    }

    #[getter]
    fn beta(&self) -> f64 {
        self.inner.angles()[1]
    }

    #[getter]
    fn gamma(&self) -> f64 {
        self.inner.angles()[2]
    }

    #[getter]
    fn volume(&self) -> f64 {
        self.inner.volume()
    }

    #[getter]
    fn parameters(&self) -> Vec<f64> {
        self.inner.parameters().to_vec()
    }

    #[getter]
    fn matrix(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        numpy_array2(
            py,
            self.inner.matrix().iter().map(|row| row.to_vec()).collect(),
        )
    }

    fn to_cart(&self, py: Python<'_>, point: Vec<f64>) -> PyResult<Py<PyAny>> {
        numpy_array1(py, self.inner.to_cart(frac(&point)?).to_vec())
    }

    fn to_frac(&self, py: Python<'_>, point: Vec<f64>) -> PyResult<Py<PyAny>> {
        numpy_array1(py, self.inner.to_frac(frac(&point)?).to_vec())
    }

    fn d_spacing(&self, hkl: Vec<f64>) -> PyResult<f64> {
        self.inner.d_spacing(frac(&hkl)?).map_err(crate::map_error)
    }

    fn distance(&self, first: Vec<f64>, second: Vec<f64>) -> PyResult<f64> {
        Ok(self.inner.distance(frac(&first)?, frac(&second)?))
    }

    fn minimum_image(&self, py: Python<'_>, delta: Vec<f64>) -> PyResult<(Py<PyAny>, f64)> {
        let (vector, distance) = self.inner.minimum_image(frac(&delta)?);
        Ok((numpy_array1(py, vector.to_vec())?, distance))
    }

    fn scaled_to_volume(&self, volume: f64) -> PyResult<Self> {
        Ok(Self {
            inner: self
                .inner
                .scaled_to_volume(volume)
                .map_err(crate::map_error)?,
        })
    }

    fn __repr__(&self) -> String {
        let [a, b, c, alpha, beta, gamma] = self.inner.parameters();
        format!(
            "Lattice(a={a:.4}, b={b:.4}, c={c:.4}, alpha={alpha:.2}, beta={beta:.2}, gamma={gamma:.2})"
        )
    }
}

/// One site of the asymmetric unit.
#[pyclass(name = "Site", module = "xtalforge")]
#[derive(Clone)]
pub(crate) struct Site {
    inner: CoreSite,
}

#[pymethods]
impl Site {
    #[new]
    fn new(label: String, element: &str, x: f64, y: f64, z: f64) -> PyResult<Self> {
        Ok(Self {
            inner: CoreSite::new(label, element, [x, y, z]).map_err(crate::map_error)?,
        })
    }

    #[getter]
    fn label(&self) -> &str {
        self.inner.label()
    }

    #[getter]
    fn element(&self) -> &str {
        self.inner.element()
    }

    #[getter]
    fn frac(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        numpy_array1(py, self.inner.frac().to_vec())
    }

    #[getter]
    fn occupancy(&self) -> f64 {
        self.inner.occupancy()
    }

    #[getter]
    fn u_iso(&self) -> Option<f64> {
        self.inner.u_iso()
    }

    #[getter]
    fn charge(&self) -> Option<f64> {
        self.inner.charge()
    }

    fn with_occupancy(&self, occupancy: f64) -> PyResult<Self> {
        Ok(Self {
            inner: self
                .inner
                .clone()
                .with_occupancy(occupancy)
                .map_err(crate::map_error)?,
        })
    }

    fn with_charge(&self, charge: f64) -> PyResult<Self> {
        Ok(Self {
            inner: self
                .inner
                .clone()
                .with_charge(charge)
                .map_err(crate::map_error)?,
        })
    }

    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);
        dict.set_item("label", self.inner.label())?;
        dict.set_item("element", self.inner.element())?;
        dict.set_item("frac", self.inner.frac().to_vec())?;
        dict.set_item("occupancy", self.inner.occupancy())?;
        dict.set_item("u_iso", self.inner.u_iso())?;
        dict.set_item("charge", self.inner.charge())?;
        Ok(dict)
    }
}

/// A crystal structure: lattice, Hall symbol and asymmetric sites.
#[pyclass(name = "Structure", module = "xtalforge")]
#[derive(Clone)]
pub(crate) struct Structure {
    inner: CoreStructure,
}

#[pymethods]
impl Structure {
    /// Reads a structure from text; an empty hint auto-detects the format.
    #[new]
    #[pyo3(signature = (text, format_hint = ""))]
    fn new(text: &str, format_hint: &str) -> PyResult<Self> {
        let hint = if format_hint.trim().is_empty() {
            crate::structures::detect_format(text)
        } else {
            format_hint
        };
        Ok(Self {
            inner: xtalforge::io::read_structure(text.trim_start_matches('\u{feff}'), hint)
                .map_err(crate::map_error)?,
        })
    }

    #[staticmethod]
    fn from_cif(text: &str) -> PyResult<Self> {
        Ok(Self {
            inner: xtalforge::io::read_cif(text.trim_start_matches('\u{feff}'))
                .map_err(crate::map_error)?,
        })
    }

    #[getter]
    fn name(&self) -> &str {
        self.inner.name()
    }

    #[getter]
    fn formula(&self) -> String {
        self.inner.formula()
    }

    #[getter]
    fn hall(&self) -> &str {
        self.inner.hall()
    }

    #[getter]
    fn lattice(&self) -> Lattice {
        Lattice {
            inner: *self.inner.lattice(),
        }
    }

    #[getter]
    fn sites(&self) -> Vec<Site> {
        self.inner
            .sites()
            .iter()
            .map(|site| Site {
                inner: site.clone(),
            })
            .collect()
    }

    #[getter]
    fn n_sites(&self) -> usize {
        self.inner.sites().len()
    }

    #[getter]
    fn volume(&self) -> f64 {
        self.inner.lattice().volume()
    }

    #[getter]
    fn density(&self) -> PyResult<f64> {
        self.inner.density().map_err(crate::map_error)
    }

    #[getter]
    fn is_p1(&self) -> bool {
        self.inner.is_p1()
    }

    fn p1(&self) -> PyResult<P1Cell> {
        Ok(P1Cell {
            inner: self.inner.p1().map_err(crate::map_error)?,
        })
    }

    #[pyo3(signature = (na, nb, nc))]
    fn supercell(&self, na: u32, nb: u32, nc: u32) -> PyResult<Self> {
        Ok(Self {
            inner: crate::structures::as_p1_structure(&self.inner)?
                .supercell(na, nb, nc)
                .map_err(crate::map_error)?,
        })
    }

    #[pyo3(signature = (h, k, l, layers, vacuum))]
    fn slab(&self, h: i32, k: i32, l: i32, layers: u32, vacuum: f64) -> PyResult<Self> {
        Ok(Self {
            inner: xtalforge::core::make_slab(&self.inner, [h, k, l], layers, vacuum, 0.0)
                .map_err(crate::map_error)?
                .structure,
        })
    }

    fn bonds<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        let graph = self.inner.bond_graph().map_err(crate::map_error)?;
        let items: Vec<Bound<'_, PyDict>> = graph
            .bonds
            .iter()
            .map(|bond| {
                let dict = PyDict::new(py);
                dict.set_item("i", bond.i)?;
                dict.set_item("j", bond.j)?;
                dict.set_item("image", bond.image.to_vec())?;
                dict.set_item("distance", bond.distance)?;
                dict.set_item("order", bond.order)?;
                dict.set_item("explicit", bond.explicit)?;
                Ok(dict)
            })
            .collect::<PyResult<_>>()?;
        PyList::new(py, items)
    }

    fn rings(&self, max_size: usize) -> PyResult<Vec<Vec<usize>>> {
        Ok(xtalforge::core::rings(&self.inner, max_size)
            .map_err(crate::map_error)?
            .into_iter()
            .map(|ring| ring.atoms)
            .collect())
    }

    fn summary<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let summary = self.inner.summary().map_err(crate::map_error)?;
        let p1 = self.inner.p1().map_err(crate::map_error)?;
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
        Ok(dict)
    }

    fn to_format(&self, format: &str) -> PyResult<String> {
        xtalforge::io::write_structure(&self.inner, format).map_err(crate::map_error)
    }

    fn to_json(&self) -> PyResult<String> {
        self.inner.to_document().to_json().map_err(crate::map_error)
    }

    fn __repr__(&self) -> String {
        format!(
            "Structure(name='{}', formula='{}', n_sites={})",
            self.inner.name(),
            self.inner.formula(),
            self.inner.sites().len()
        )
    }
}

/// The expanded P1 cell of a structure.
#[pyclass(name = "P1Cell", module = "xtalforge")]
#[derive(Clone)]
pub(crate) struct P1Cell {
    inner: CoreP1Cell,
}

#[pymethods]
impl P1Cell {
    #[getter]
    fn n_atoms(&self) -> usize {
        self.inner.n_atoms()
    }

    #[getter]
    fn lattice(&self) -> Lattice {
        Lattice {
            inner: self.inner.lattice,
        }
    }

    #[getter]
    fn elements(&self) -> Vec<String> {
        self.inner.elements()
    }

    #[getter]
    fn frac(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        numpy_array2(
            py,
            self.inner
                .frac()
                .iter()
                .map(|point| point.to_vec())
                .collect(),
        )
    }

    #[getter]
    fn occupancy(&self) -> Vec<f64> {
        self.inner.occupancy()
    }

    #[getter]
    fn atoms<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        let items: Vec<Bound<'_, PyDict>> = self
            .inner
            .atoms
            .iter()
            .map(|atom| atom_dict(py, atom))
            .collect::<PyResult<_>>()?;
        PyList::new(py, items)
    }

    fn __repr__(&self) -> String {
        format!("P1Cell(n_atoms={})", self.inner.n_atoms())
    }
}

fn atom_dict<'py>(py: Python<'py>, atom: &P1Atom) -> PyResult<Bound<'py, PyDict>> {
    let dict = PyDict::new(py);
    dict.set_item("element", &atom.element)?;
    dict.set_item("label", &atom.label)?;
    dict.set_item("frac", atom.frac.to_vec())?;
    dict.set_item("occupancy", atom.occupancy)?;
    dict.set_item("site_idx", atom.site_idx)?;
    dict.set_item("op_idx", atom.op_idx)?;
    dict.set_item("tau", atom.tau.to_vec())?;
    Ok(dict)
}
