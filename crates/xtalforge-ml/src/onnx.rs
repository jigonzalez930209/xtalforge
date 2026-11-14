#![allow(clippy::indexing_slicing, clippy::needless_range_loop)] // fixed 3x3 strain and Cartesian axes

//! Native ONNX inference through tract (Task 21.3.1).
//!
//! The `onnx` feature adds a pure-Rust ONNX runtime. [`OnnxModel`] wraps a
//! single-input, single-output graph whose input is a feature vector of shape
//! `[1, n_features]` and whose output is one or more scalars. [`OnnxCalculator`]
//! turns such a model into a [`Calculator`] for
//! descriptor-based ML potentials: the energy is the first output, forces are
//! central differences in Cartesian coordinates, and the stress is a central
//! difference in the six symmetric strain components.
//!
//! Message-passing models (MACE, ORB, MatterSim) are not single-vector graphs;
//! the Python ML bridge of Phase 8 remains the route for those.

use std::sync::Arc;

use serde::{Deserialize, Serialize};
use tract_onnx::prelude::*;
use xtalforge_core::{Error, Site, Structure, StructureBuilder};
use xtalforge_ff::Calculator;

use crate::soap::{SoapOptions, soap_power_spectrum};

type Runner = dyn Fn(&[f64]) -> Result<Vec<f64>, Error> + Send + Sync;

/// Options for ONNX loading.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OnnxOptions {
    /// Number of input features.
    pub n_features: usize,
    /// Input tensor name; `None` uses the first model input.
    pub input: Option<String>,
}

/// A loaded ONNX graph.
#[derive(Clone)]
pub struct OnnxModel {
    runner: Arc<Runner>,
    n_features: usize,
    n_outputs: usize,
}

impl std::fmt::Debug for OnnxModel {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("OnnxModel")
            .field("n_features", &self.n_features)
            .field("n_outputs", &self.n_outputs)
            .finish_non_exhaustive()
    }
}

impl OnnxModel {
    /// Loads a model from a file.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for I/O, parse or shape errors.
    pub fn from_path(
        path: impl AsRef<std::path::Path>,
        options: &OnnxOptions,
    ) -> Result<Self, Error> {
        let text = std::fs::read(path.as_ref()).map_err(|error| {
            Error::InvalidInput(format!(
                "could not read {}: {error}",
                path.as_ref().display()
            ))
        })?;
        Self::from_bytes(&text, options)
    }

    /// Loads a model from bytes.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for parse or shape errors.
    pub fn from_bytes(bytes: &[u8], options: &OnnxOptions) -> Result<Self, Error> {
        if options.n_features == 0 {
            return Err(Error::InvalidInput(
                "the ONNX input width must be at least one".to_owned(),
            ));
        }
        let mut cursor = std::io::Cursor::new(bytes);
        let model = tract_onnx::onnx()
            .model_for_read(&mut cursor)
            .map_err(|error| Error::InvalidInput(format!("ONNX parse: {error}")))?;
        let model = if let Some(name) = &options.input {
            model
                .with_input_fact(0, f32::fact([1, options.n_features as i64]).into())
                .map_err(|error| Error::InvalidInput(format!("ONNX input fact: {error}")))?
                .with_input_names([name.clone()])
                .map_err(|error| Error::InvalidInput(format!("ONNX input name: {error}")))?
        } else {
            model
                .with_input_fact(0, f32::fact([1, options.n_features as i64]).into())
                .map_err(|error| Error::InvalidInput(format!("ONNX input fact: {error}")))?
        };
        let plan = model
            .into_optimized()
            .and_then(|model| model.into_runnable())
            .map_err(|error| Error::InvalidInput(format!("ONNX optimize: {error}")))?;
        let n_features = options.n_features;
        let runner = move |features: &[f64]| -> Result<Vec<f64>, Error> {
            if features.len() != n_features {
                return Err(Error::InvalidInput(format!(
                    "the ONNX model expects {n_features} features, got {}",
                    features.len()
                )));
            }
            let input = tract_ndarray::Array2::from_shape_vec(
                (1, n_features),
                features.iter().map(|value| *value as f32).collect(),
            )
            .map_err(|error| Error::InvalidInput(format!("ONNX input array: {error}")))?;
            let tensor: Tensor = input.into();
            let outputs = plan
                .run(tvec!(tensor.into()))
                .map_err(|error| Error::InvalidInput(format!("ONNX run: {error}")))?;
            let Some(output) = outputs.first() else {
                return Err(Error::InvalidInput(
                    "the ONNX model returned no output".to_owned(),
                ));
            };
            let view = output
                .to_plain_array_view::<f32>()
                .map_err(|error| Error::InvalidInput(format!("ONNX output: {error}")))?;
            Ok(view.iter().map(|value| f64::from(*value)).collect())
        };
        Ok(Self {
            runner: Arc::new(runner),
            n_features,
            n_outputs: 1,
        })
    }

    /// Number of input features.
    #[must_use]
    pub fn n_features(&self) -> usize {
        self.n_features
    }

    /// Number of output scalars.
    #[must_use]
    pub fn n_outputs(&self) -> usize {
        self.n_outputs
    }

    /// Runs the model on one feature vector.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a wrong input width or an inference
    /// failure.
    pub fn predict(&self, features: &[f64]) -> Result<Vec<f64>, Error> {
        (self.runner)(features)
    }
}

/// Which descriptor the calculator feeds to the model.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum OnnxDescriptor {
    /// Average SOAP power spectrum.
    Soap(SoapOptions),
    /// Flattened Cartesian coordinates of the P1 atoms.
    Cartesian,
}

/// A descriptor-based ONNX potential.
pub struct OnnxCalculator {
    model: OnnxModel,
    descriptor: OnnxDescriptor,
    displacement: f64,
    strain: f64,
}

impl OnnxCalculator {
    /// Builds a calculator.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a non-positive finite difference.
    pub fn new(
        model: OnnxModel,
        descriptor: OnnxDescriptor,
        displacement: f64,
        strain: f64,
    ) -> Result<Self, Error> {
        if !displacement.is_finite() || displacement <= 0.0 || !strain.is_finite() || strain <= 0.0
        {
            return Err(Error::InvalidInput(
                "ONNX finite-difference steps must be positive and finite".to_owned(),
            ));
        }
        Ok(Self {
            model,
            descriptor,
            displacement,
            strain,
        })
    }

    /// The underlying model.
    #[must_use]
    pub fn model(&self) -> &OnnxModel {
        &self.model
    }

    fn features(&self, structure: &Structure) -> Result<Vec<f64>, Error> {
        match &self.descriptor {
            OnnxDescriptor::Soap(options) => {
                let descriptor = soap_power_spectrum(structure, options)?;
                descriptor
                    .rows
                    .first()
                    .cloned()
                    .ok_or_else(|| Error::InvalidInput("SOAP produced no row".to_owned()))
            }
            OnnxDescriptor::Cartesian => {
                // Keep the coordinates as supplied for a P1 structure so a
                // finite-difference step across the cell boundary does not
                // wrap by a full lattice vector.
                if structure.is_p1() {
                    let mut features = Vec::with_capacity(3 * structure.sites().len());
                    for atom in structure.sites() {
                        features.extend(structure.lattice().to_cart(atom.frac()));
                    }
                    Ok(features)
                } else {
                    let p1 = structure.p1()?;
                    let mut features = Vec::with_capacity(3 * p1.atoms.len());
                    for atom in &p1.atoms {
                        features.extend(p1.lattice.to_cart(atom.frac));
                    }
                    Ok(features)
                }
            }
        }
    }

    fn energy(&self, structure: &Structure) -> Result<f64, Error> {
        let features = self.features(structure)?;
        let outputs = self.model.predict(&features)?;
        outputs.first().copied().ok_or_else(|| {
            Error::InvalidInput("the ONNX model returned an empty output".to_owned())
        })
    }
}

impl Calculator for OnnxCalculator {
    fn name(&self) -> &'static str {
        "onnx"
    }

    fn energy(&self, structure: &Structure) -> Result<f64, xtalforge_ff::Error> {
        self.energy_checked(structure)
            .map_err(|error| xtalforge_ff::Error::InvalidInput(error.message()))
    }

    fn energy_forces(
        &self,
        structure: &Structure,
    ) -> Result<(f64, Vec<[f64; 3]>), xtalforge_ff::Error> {
        let p1 = structure
            .p1()
            .map_err(|error| xtalforge_ff::Error::InvalidInput(error.message()))?;
        let base = self
            .energy_checked(structure)
            .map_err(|error| xtalforge_ff::Error::InvalidInput(error.message()))?;
        let mut forces = Vec::with_capacity(p1.atoms.len());
        for index in 0..p1.atoms.len() {
            let mut gradient = [0.0; 3];
            for axis in 0..3 {
                let plus = build_p1(structure, &p1, index, axis, self.displacement)?;
                let minus = build_p1(structure, &p1, index, axis, -self.displacement)?;
                let plus_energy = self.energy(&plus)?;
                let minus_energy = self.energy(&minus)?;
                gradient[axis] = (plus_energy - minus_energy) / (2.0 * self.displacement);
            }
            forces.push([-gradient[0], -gradient[1], -gradient[2]]);
        }
        Ok((base, forces))
    }

    fn provides_forces(&self) -> bool {
        true
    }

    fn provides_stress(&self) -> bool {
        true
    }

    fn stress(&self, structure: &Structure) -> Result<Option<[[f64; 3]; 3]>, xtalforge_ff::Error> {
        let base = self
            .energy_checked(structure)
            .map_err(|error| xtalforge_ff::Error::InvalidInput(error.message()))?;
        let matrix = structure.lattice().matrix();
        let volume = structure.lattice().volume();
        let mut tensor = [[0.0; 3]; 3];
        for (row, column) in [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)] {
            let mut plus = matrix;
            plus[row][column] += self.strain * plus[row][column];
            let plus_structure = replace_lattice(structure, plus)?;
            let plus_energy = self.energy(&plus_structure)?;
            let mut minus = matrix;
            minus[row][column] -= self.strain * minus[row][column];
            let minus_structure = replace_lattice(structure, minus)?;
            let minus_energy = self.energy(&minus_structure)?;
            let derivative = (plus_energy - minus_energy) / (2.0 * self.strain);
            let value = derivative / volume * 160.217_663_4;
            if row == column {
                tensor[row][column] = value;
            } else {
                tensor[row][column] = value;
                tensor[column][row] = value;
            }
        }
        let _ = base;
        Ok(Some(tensor))
    }
}

impl OnnxCalculator {
    fn energy_checked(&self, structure: &Structure) -> Result<f64, Error> {
        self.energy(structure)
    }
}

fn build_p1(
    structure: &Structure,
    p1: &xtalforge_core::P1Cell,
    index: usize,
    axis: usize,
    offset: f64,
) -> Result<Structure, xtalforge_ff::Error> {
    let mut sites = Vec::with_capacity(p1.atoms.len());
    for (atom_index, atom) in p1.atoms.iter().enumerate() {
        let mut frac = atom.frac;
        if atom_index == index {
            let cart = p1.lattice.to_cart(frac);
            let mut moved = cart;
            moved[axis] += offset;
            frac = p1.lattice.to_frac(moved);
        }
        sites.push(
            Site::new(&atom.label, &atom.element, frac)
                .map_err(|error| xtalforge_ff::Error::InvalidInput(error.message()))?,
        );
    }
    let mut builder = StructureBuilder::new(structure.name(), p1.lattice);
    for site in sites {
        builder = builder.site(site);
    }
    builder
        .build()
        .map_err(|error| xtalforge_ff::Error::InvalidInput(error.message()))
}

fn replace_lattice(
    structure: &Structure,
    matrix: [[f64; 3]; 3],
) -> Result<Structure, xtalforge_ff::Error> {
    let lattice = xtalforge_core::Lattice::new(matrix)
        .map_err(|error| xtalforge_ff::Error::InvalidInput(error.message()))?;
    Ok(structure.clone().with_lattice(lattice))
}
