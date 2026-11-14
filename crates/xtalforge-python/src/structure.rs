//! Python object API for `Structure` (Task 12.1.1).

use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

use xtalforge::core::{Lattice as CoreLattice, Site as CoreSite, Structure as CoreStructure};

use crate::classes::{Lattice, P1Cell, Site, numpy_array1, numpy_array2, triple};

type RingPath = Vec<Vec<(usize, (i32, i32, i32))>>;

/// Looks up one public function of the pure-Python adapter module.
fn adapter_function<'py>(py: Python<'py>, name: &str) -> PyResult<Bound<'py, PyAny>> {
    py.import("xtalforge.adapters")?.getattr(name)
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

    /// Enters a ``with`` block and returns the same structure.
    fn __enter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    /// Leaves a ``with`` block; exceptions are not suppressed.
    #[pyo3(signature = (_exc_type, _exc_value, _traceback))]
    fn __exit__(
        &self,
        _exc_type: Option<&Bound<'_, PyAny>>,
        _exc_value: Option<&Bound<'_, PyAny>>,
        _traceback: Option<&Bound<'_, PyAny>>,
    ) -> bool {
        false
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

    fn ring_paths(&self, max_size: usize) -> PyResult<RingPath> {
        Ok(xtalforge::core::ring_paths(&self.inner, max_size)
            .map_err(crate::map_error)?
            .into_iter()
            .map(|path| {
                path.into_iter()
                    .map(|step| (step.atom, (step.image[0], step.image[1], step.image[2])))
                    .collect()
            })
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

    /// The legacy dictionary shape for this structure.
    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        crate::structures::structure_payload(py, &self.inner)
    }

    /// Writes the structure to a path chosen by its extension.
    fn write(&self, path: &str) -> PyResult<()> {
        xtalforge::write(path, &self.inner).map_err(crate::map_error)
    }

    /// Cartesian positions of the P1 atoms as a NumPy array.
    #[getter]
    fn positions(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let p1 = self.inner.p1().map_err(crate::map_error)?;
        numpy_array2(
            py,
            p1.atoms
                .iter()
                .map(|atom| p1.lattice.to_cart(atom.frac).to_vec())
                .collect(),
        )
    }

    /// Fractional positions of the P1 atoms as a NumPy array.
    #[getter]
    fn frac_positions(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let p1 = self.inner.p1().map_err(crate::map_error)?;
        numpy_array2(py, p1.atoms.iter().map(|atom| atom.frac.to_vec()).collect())
    }

    /// Cell matrix as a NumPy array.
    #[getter]
    fn cell(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        numpy_array2(
            py,
            self.inner
                .lattice()
                .matrix()
                .iter()
                .map(|row| row.to_vec())
                .collect(),
        )
    }

    /// Atomic numbers of the P1 atoms as a NumPy array.
    #[getter]
    fn numbers(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let p1 = self.inner.p1().map_err(crate::map_error)?;
        let mut numbers = Vec::with_capacity(p1.atoms.len());
        for atom in &p1.atoms {
            let element = xtalforge::core::element_info(&atom.element).map_err(crate::map_error)?;
            numbers.push(f64::from(element.z));
        }
        numpy_array1(py, numbers)
    }

    /// Builds a P1 structure from cell, atomic numbers and positions.
    #[staticmethod]
    #[pyo3(signature = (cell, numbers, frac = None, cart = None))]
    fn from_arrays(
        cell: Vec<Vec<f64>>,
        numbers: Vec<u8>,
        frac: Option<Vec<Vec<f64>>>,
        cart: Option<Vec<Vec<f64>>>,
    ) -> PyResult<Self> {
        if cell.len() != 3 {
            return Err(crate::invalid_input("cell must have three rows"));
        }
        let mut rows = [[0.0; 3]; 3];
        for (slot, row) in rows.iter_mut().zip(cell.iter()) {
            *slot = triple(row)?;
        }
        let lattice = CoreLattice::new(rows).map_err(crate::map_error)?;
        let positions = match (frac, cart) {
            (Some(values), _) => values,
            (None, Some(values)) => {
                let mut positions = Vec::with_capacity(values.len());
                for point in values {
                    positions.push(lattice.to_frac(triple(&point)?).to_vec());
                }
                positions
            }
            (None, None) => {
                return Err(crate::invalid_input(
                    "from_arrays needs frac or cart positions",
                ));
            }
        };
        if positions.len() != numbers.len() {
            return Err(crate::invalid_input(
                "positions and numbers must have the same length",
            ));
        }
        let mut builder = xtalforge::core::StructureBuilder::new("structure", lattice);
        for (index, (&z, point)) in numbers.iter().zip(&positions).enumerate() {
            let symbol = xtalforge::core::elements::element_by_z(z)
                .map_err(crate::map_error)?
                .symbol
                .clone();
            builder = builder.site(
                CoreSite::new(format!("{symbol}{}", index + 1), &symbol, triple(point)?)
                    .map_err(crate::map_error)?,
            );
        }
        Ok(Self {
            inner: builder.build().map_err(crate::map_error)?,
        })
    }

    /// Returns a copy with new P1 fractional positions.
    fn set_positions(&self, frac_positions: Vec<Vec<f64>>) -> PyResult<Self> {
        let p1 = self.inner.p1().map_err(crate::map_error)?;
        if frac_positions.len() != p1.atoms.len() {
            return Err(crate::invalid_input(
                "set_positions needs one point per P1 atom",
            ));
        }
        let mut builder =
            xtalforge::core::StructureBuilder::new(self.inner.name(), *self.inner.lattice())
                .hall("P 1");
        for (atom, point) in p1.atoms.iter().zip(&frac_positions) {
            builder = builder.site(
                CoreSite::new(atom.label.clone(), &atom.element, triple(point)?)
                    .map_err(crate::map_error)?
                    .with_occupancy(atom.occupancy)
                    .map_err(crate::map_error)?,
            );
        }
        Ok(Self {
            inner: builder.build().map_err(crate::map_error)?,
        })
    }

    /// Returns the UFF energy in eV.
    fn ff_energy(&self) -> PyResult<f64> {
        let calculator =
            xtalforge::ff::UffCalculator::new(xtalforge::ff::CalculatorOptions::default())
                .map_err(|error| crate::invalid_input(error.to_string()))?;
        xtalforge::ff::Calculator::energy(&calculator, &self.inner).map_err(crate::map_ff_error)
    }

    /// Adds the missing hydrogens and returns a new structure.
    #[pyo3(signature = (distance = 1.09))]
    fn add_hydrogens(&self, distance: f64) -> PyResult<Self> {
        let options = xtalforge::build::HydrogenOptions {
            distance,
            ..xtalforge::build::HydrogenOptions::default()
        };
        Ok(Self {
            inner: xtalforge::build::add_hydrogens(&self.inner, &options)
                .map_err(crate::map_error)?,
        })
    }

    fn _repr_html_(&self) -> String {
        format!(
            "<b>{}</b><br>{}<br>{} sites, V={:.3} Å³",
            self.inner.name(),
            self.inner.formula(),
            self.inner.sites().len(),
            self.inner.lattice().volume()
        )
    }

    fn to_format(&self, format: &str) -> PyResult<String> {
        xtalforge::io::write_structure(&self.inner, format).map_err(crate::map_error)
    }

    fn to_json(&self) -> PyResult<String> {
        self.inner.to_document().to_json().map_err(crate::map_error)
    }

    /// A padded P1 box around a finite molecule (adapter support).
    #[staticmethod]
    #[pyo3(signature = (name, points, padding = 10.0))]
    fn _molecule(name: &str, points: Vec<(String, Vec<f64>)>, padding: f64) -> PyResult<Self> {
        let mut converted = Vec::with_capacity(points.len());
        for (element, point) in &points {
            converted.push((element.clone(), crate::classes::triple(point)?));
        }
        Ok(Self {
            inner: CoreStructure::molecule(name, &converted, padding).map_err(crate::map_error)?,
        })
    }

    /// Builds a structure from a core JSON document or a bare structure object.
    #[staticmethod]
    fn _from_json(text: &str) -> PyResult<Self> {
        Ok(Self {
            inner: crate::structures::load_structure(text, None)?,
        })
    }

    /// The explicit P1 cell as a versioned structure JSON document.
    fn _p1_json(&self) -> PyResult<String> {
        self.inner
            .as_p1()
            .map_err(crate::map_error)?
            .to_document()
            .to_json()
            .map_err(crate::map_error)
    }

    /// UFF energy in eV and forces in eV/Angstrom over the P1 atoms.
    fn _ff_energy_forces(&self) -> PyResult<(f64, Vec<Vec<f64>>)> {
        let calculator =
            xtalforge::ff::UffCalculator::new(xtalforge::ff::CalculatorOptions::default())
                .map_err(crate::map_ff_error)?;
        let (energy, forces) = xtalforge::ff::Calculator::energy_forces(&calculator, &self.inner)
            .map_err(crate::map_ff_error)?;
        Ok((
            energy,
            forces.into_iter().map(|force| force.to_vec()).collect(),
        ))
    }

    /// Converts an ASE ``Atoms`` object into a P1 structure (Task 13.1.1).
    #[staticmethod]
    fn from_ase(py: Python<'_>, atoms: &Bound<'_, PyAny>) -> PyResult<Self> {
        let text: String = adapter_function(py, "ase_to_structure_json")?
            .call1((atoms,))?
            .extract()?;
        Self::_from_json(&text)
    }

    /// Converts this structure into an ASE ``Atoms`` object (Task 13.1.1).
    fn to_ase(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        Ok(adapter_function(py, "structure_to_ase")?
            .call1((self._p1_json()?,))?
            .unbind())
    }

    /// Converts a pymatgen ``Structure`` or ``Molecule`` (Task 13.1.2).
    #[staticmethod]
    fn from_pymatgen(py: Python<'_>, structure_or_molecule: &Bound<'_, PyAny>) -> PyResult<Self> {
        let text: String = adapter_function(py, "pymatgen_to_structure_json")?
            .call1((structure_or_molecule,))?
            .extract()?;
        Self::_from_json(&text)
    }

    /// Converts this structure into a pymatgen ``Structure`` or ``Molecule``.
    fn to_pymatgen(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        Ok(adapter_function(py, "structure_to_pymatgen")?
            .call1((self._p1_json()?,))?
            .unbind())
    }

    /// Converts this structure into a phonopy ``PhonopyAtoms`` (Task 13.1.3).
    fn to_phonopy_atoms(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        Ok(adapter_function(py, "structure_to_phonopy_atoms")?
            .call1((self._p1_json()?,))?
            .unbind())
    }

    /// Converts an RDKit ``Mol`` with a conformer into a P1 structure.
    #[staticmethod]
    fn from_rdkit(py: Python<'_>, mol: &Bound<'_, PyAny>) -> PyResult<Self> {
        let text: String = adapter_function(py, "rdkit_to_structure_json")?
            .call1((mol,))?
            .extract()?;
        Self::_from_json(&text)
    }

    /// Converts this structure into an RDKit ``Mol`` with a conformer.
    fn to_rdkit(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        Ok(adapter_function(py, "structure_to_rdkit")?
            .call1((self._p1_json()?,))?
            .unbind())
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
