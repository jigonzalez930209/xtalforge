//! Extra carbon and polymer builder bindings (Tasks 6.5, 6.7-6.9).

use pyo3::prelude::*;
use pyo3::types::PyDict;

use xtalforge::build::{
    CarbonRecipe, PolymerOptions, ZtcOptions, build_carbon as carbon_build,
    build_polymer as polymer_build, build_schwarzite as schwarzite_build, build_ztc as ztc_build,
    commensurate_twist, face_spiral, nanoribbon,
};

use crate::structures::structure_payload;

fn ztc_options(samples: usize, cell: f64, isovalue: f64) -> ZtcOptions {
    ZtcOptions {
        samples,
        cell,
        isovalue,
    }
}

/// Builds a C20 dodecahedron.
#[pyfunction]
pub(crate) fn build_c20<'py>(py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
    let structure = face_spiral(20).map_err(crate::map_error)?;
    structure_payload(py, &structure)
}

/// Builds two twisted graphene layers.
#[pyfunction]
#[pyo3(signature = (n, m, theta, layers = 2))]
pub(crate) fn twist<'py>(
    py: Python<'py>,
    n: u32,
    m: u32,
    theta: f64,
    layers: u32,
) -> PyResult<Bound<'py, PyDict>> {
    let structure = commensurate_twist(n, m, theta, layers).map_err(crate::map_error)?;
    structure_payload(py, &structure)
}

/// Builds a honeycomb nanoribbon with vacuum.
#[pyfunction]
#[pyo3(signature = (width, length, vacuum = 12.0))]
pub(crate) fn ribbon<'py>(
    py: Python<'py>,
    width: u32,
    length: u32,
    vacuum: f64,
) -> PyResult<Bound<'py, PyDict>> {
    let structure = nanoribbon(width, length, vacuum).map_err(crate::map_error)?;
    structure_payload(py, &structure)
}

/// Builds a gyroid ZTC carbon model.
#[pyfunction]
#[pyo3(signature = (samples = 16, cell = 12.0, isovalue = 0.0))]
pub(crate) fn build_ztc<'py>(
    py: Python<'py>,
    samples: usize,
    cell: f64,
    isovalue: f64,
) -> PyResult<Bound<'py, PyDict>> {
    let structure = ztc_build(&ztc_options(samples, cell, isovalue)).map_err(crate::map_error)?;
    structure_payload(py, &structure)
}

/// Builds a Schwarz P carbon model.
#[pyfunction]
#[pyo3(signature = (samples = 16, cell = 12.0, isovalue = 0.0))]
pub(crate) fn build_schwarzite<'py>(
    py: Python<'py>,
    samples: usize,
    cell: f64,
    isovalue: f64,
) -> PyResult<Bound<'py, PyDict>> {
    let structure =
        schwarzite_build(&ztc_options(samples, cell, isovalue)).map_err(crate::map_error)?;
    structure_payload(py, &structure)
}

/// Builds a polymer chain and records the build steps.
#[pyfunction]
#[pyo3(signature = (smiles, monomers, density = 1.0))]
pub(crate) fn build_polymer<'py>(
    py: Python<'py>,
    smiles: &str,
    monomers: u32,
    density: f64,
) -> PyResult<Bound<'py, PyDict>> {
    let build =
        polymer_build(smiles, monomers, &PolymerOptions { density }).map_err(crate::map_error)?;
    let dict = structure_payload(py, &build.structure)?;
    dict.set_item("steps", build.steps)?;
    Ok(dict)
}

/// Builds a disordered carbon from a recipe.
///
/// Returns the structure summary plus the density, mean C-C bond, edge
/// ratios, ring census and the builder notes.
#[pyfunction]
#[allow(clippy::too_many_arguments)] // one keyword per recipe field
#[pyo3(signature = (net = "dia", repeat = (1, 1, 1), density = 0.42, radius_ratio = 0.25, coverage = 0.38, seed = 0, hydrogen = 0.07, fluorine = 0.29, oxygen = 0.044))]
pub(crate) fn build_carbon<'py>(
    py: Python<'py>,
    net: &str,
    repeat: (usize, usize, usize),
    density: f64,
    radius_ratio: f64,
    coverage: f64,
    seed: u64,
    hydrogen: f64,
    fluorine: f64,
    oxygen: f64,
) -> PyResult<Bound<'py, PyDict>> {
    let recipe = CarbonRecipe {
        net: net.to_owned(),
        repeat: [repeat.0, repeat.1, repeat.2],
        density,
        radius_ratio,
        coverage,
        seed,
        hydrogen,
        fluorine,
        oxygen,
        ..CarbonRecipe::default()
    };
    let build = carbon_build(&recipe).map_err(crate::map_error)?;
    let dict = structure_payload(py, &build.structure)?;
    dict.set_item("density", build.density)?;
    dict.set_item("mean_bond", build.mean_bond)?;
    dict.set_item("edge_fraction", build.edge_fraction)?;
    dict.set_item("pieces", build.pieces)?;
    dict.set_item("periodicity", build.periodicity)?;
    dict.set_item("ratios", build.ratios.to_vec())?;
    dict.set_item("census", build.census)?;
    dict.set_item("notes", build.notes)?;
    Ok(dict)
}
