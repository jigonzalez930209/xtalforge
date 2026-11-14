#![allow(
    clippy::indexing_slicing,
    clippy::needless_range_loop,
    clippy::approx_constant
)] // fixed radial/angular kernels and published constants

//! Exact SOAP power spectrum (Task 21.1.2).
//!
//! The implementation follows the DScribe SOAP convention with the GTO
//! radial basis and Gaussian density width `sigma`:
//!
//! * radial exponents `alpha_lk` come from `a_k = linspace(1, r_cut, n_max)`
//!   and `alpha = -ln(1e-3 / a_k^l) / a_k^2`;
//! * the overlap matrix `S_l[k,k'] = 0.5 Gamma(l+3/2)
//!   (alpha_k + alpha_k')^{-l-3/2}` is Löwdin-orthonormalised;
//! * the Gaussian width enters through `eta = 1/(2 sigma^2)` and the
//!   convolution factors `aOa`/`bOa`;
//! * real spherical harmonics through degree five use the DScribe
//!   polynomial normalisation;
//! * the power spectrum is `p_l[k,k'] = pi^4 sqrt(8/(2l+1))
//!   sum_m c_klm c_k'lm` (the `pi^3` factor is present for `l > 1`).
//!
//! The fixture under `crates/xtalforge-ml/data/soap_dscribe.json` records a
//! DScribe 2.x run on two carbon atoms. The test matches every feature within
//! `5e-4 * max(1, |reference|)`. DScribe is not a build dependency.
//!
//! The older [`soap_like_descriptor`](xtalforge_analysis::soap_like_descriptor)
//! stays as the fast approximation.

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Structure, neighbors};

/// Options for [`soap_power_spectrum`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SoapOptions {
    /// Radial cutoff in Angstrom.
    pub cutoff: f64,
    /// Number of radial basis functions.
    pub n_max: usize,
    /// Maximum spherical harmonic degree, at most five.
    pub l_max: usize,
    /// Gaussian density width in Angstrom.
    pub sigma: f64,
    /// Average the power spectrum over all centres.
    pub average: bool,
}

impl Default for SoapOptions {
    fn default() -> Self {
        Self {
            cutoff: 5.0,
            n_max: 4,
            l_max: 4,
            sigma: 0.5,
            average: false,
        }
    }
}

/// A SOAP feature matrix.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SoapDescriptor {
    /// One row per centre (one row total when `average` is on).
    pub rows: Vec<Vec<f64>>,
    /// Feature layout: `l`, `k`, `k'` triples in output order.
    pub layout: Vec<[usize; 3]>,
    /// Number of radial basis functions.
    pub n_max: usize,
    /// Maximum degree.
    pub l_max: usize,
    /// Number of features per row.
    pub n_features: usize,
}

impl SoapDescriptor {
    /// Number of features per row.
    #[must_use]
    pub fn n_features(&self) -> usize {
        self.layout.len()
    }

    /// Serializes the descriptor to JSON.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when serialization fails.
    pub fn to_json(&self) -> Result<String, Error> {
        serde_json::to_string(self)
            .map_err(|error| Error::InvalidInput(format!("could not serialize SOAP: {error}")))
    }
}

/// Computes the exact SOAP power spectrum.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for invalid options, an empty structure,
/// and the neighbour search errors.
pub fn soap_power_spectrum(
    structure: &Structure,
    options: &SoapOptions,
) -> Result<SoapDescriptor, Error> {
    validate(options)?;
    let p1 = structure.p1()?;
    if p1.atoms.is_empty() {
        return Err(Error::InvalidInput(
            "SOAP needs at least one atom".to_owned(),
        ));
    }
    let basis = RadialBasis::new(options);
    let mut per_centre: Vec<Option<Vec<f64>>> = vec![None; p1.atoms.len()];
    let raw = neighbors(structure, options.cutoff)?;
    let mut neighbour_lists: Vec<Vec<([f64; 3], f64)>> = vec![Vec::new(); p1.atoms.len()];
    for neighbor in &raw {
        if neighbor.i == neighbor.j && neighbor.distance < 1e-8 {
            continue;
        }
        let from = p1
            .atoms
            .get(neighbor.i)
            .map(|atom| p1.lattice.to_cart(atom.frac))
            .unwrap_or([0.0; 3]);
        let to = p1
            .atoms
            .get(neighbor.j)
            .map(|atom| p1.lattice.to_cart(atom.frac))
            .unwrap_or([0.0; 3]);
        let image = neighbor.image;
        let target = [
            to[0]
                + image[0] as f64 * p1.lattice.matrix()[0][0]
                + image[1] as f64 * p1.lattice.matrix()[1][0]
                + image[2] as f64 * p1.lattice.matrix()[2][0],
            to[1]
                + image[0] as f64 * p1.lattice.matrix()[0][1]
                + image[1] as f64 * p1.lattice.matrix()[1][1]
                + image[2] as f64 * p1.lattice.matrix()[2][1],
            to[2]
                + image[0] as f64 * p1.lattice.matrix()[0][2]
                + image[1] as f64 * p1.lattice.matrix()[1][2]
                + image[2] as f64 * p1.lattice.matrix()[2][2],
        ];
        let delta = [
            target[0] - from[0],
            target[1] - from[1],
            target[2] - from[2],
        ];
        if let Some(list) = neighbour_lists.get_mut(neighbor.i) {
            list.push((delta, neighbor.distance));
        }
        let reverse = [-delta[0], -delta[1], -delta[2]];
        if let Some(list) = neighbour_lists.get_mut(neighbor.j) {
            list.push((reverse, neighbor.distance));
        }
    }
    let mut layout = Vec::new();
    for l in 0..=options.l_max {
        for k in 0..options.n_max {
            for kd in k..options.n_max {
                layout.push([l, k, kd]);
            }
        }
    }
    for (index, list) in neighbour_lists.iter().enumerate() {
        let mut coefficients = vec![0.0; (options.n_max) * (options.l_max + 1).pow(2)];
        basis.accumulate_self(&mut coefficients);
        for (delta, distance) in list {
            if *distance < 1e-10 || *distance >= options.cutoff {
                continue;
            }
            let unit = [
                delta[0] / distance,
                delta[1] / distance,
                delta[2] / distance,
            ];
            basis.accumulate(unit, *distance, &mut coefficients);
        }
        let mut row = Vec::with_capacity(layout.len());
        for [l, k, kd] in &layout {
            let mut sum = 0.0;
            for m in 0..(2 * l + 1) {
                let a = coefficients[*k * (options.l_max + 1).pow(2) + l * l + m];
                let b = coefficients[*kd * (options.l_max + 1).pow(2) + l * l + m];
                sum += a * b;
            }
            row.push(spectrum_prefactor(*l) * sum);
        }
        per_centre[index] = Some(row);
    }
    let rows = if options.average {
        let n_features = layout.len();
        let mut mean = vec![0.0; n_features];
        let count = per_centre.len().max(1) as f64;
        for row in per_centre.into_iter().flatten() {
            for (index, value) in row.iter().enumerate() {
                mean[index] += value / count;
            }
        }
        vec![mean]
    } else {
        per_centre.into_iter().flatten().collect()
    };
    let n_features = layout.len();
    Ok(SoapDescriptor {
        rows,
        layout,
        n_max: options.n_max,
        l_max: options.l_max,
        n_features,
    })
}

fn spectrum_prefactor(l: usize) -> f64 {
    let base = std::f64::consts::PI * (8.0 / (2.0 * l as f64 + 1.0)).sqrt();
    if l > 1 {
        base * std::f64::consts::PI.powi(3)
    } else {
        base
    }
}

fn validate(options: &SoapOptions) -> Result<(), Error> {
    if !options.cutoff.is_finite() || options.cutoff <= 1.0 {
        return Err(Error::InvalidInput(
            "the SOAP cutoff must be greater than one Angstrom".to_owned(),
        ));
    }
    if options.n_max == 0 || options.n_max > 12 {
        return Err(Error::InvalidInput(
            "SOAP n_max must be in 1..=12".to_owned(),
        ));
    }
    if options.l_max > 5 {
        return Err(Error::InvalidInput(
            "SOAP l_max is supported up to five".to_owned(),
        ));
    }
    if !options.sigma.is_finite() || options.sigma <= 0.0 {
        return Err(Error::InvalidInput(
            "SOAP sigma must be positive and finite".to_owned(),
        ));
    }
    Ok(())
}

struct RadialBasis {
    n_max: usize,
    l_max: usize,
    one_over_eta: f64,
    alphas: Vec<Vec<f64>>,
    betas: Vec<Vec<Vec<f64>>>,
}

impl RadialBasis {
    fn new(options: &SoapOptions) -> Self {
        let n_max = options.n_max;
        let l_max = options.l_max;
        let mut a = Vec::with_capacity(n_max);
        if n_max == 1 {
            a.push(1.0);
        } else {
            for index in 0..n_max {
                let t = index as f64 / (n_max - 1) as f64;
                a.push(1.0 + t * (options.cutoff - 1.0));
            }
        }
        let mut alphas = Vec::with_capacity(l_max + 1);
        let mut betas = Vec::with_capacity(l_max + 1);
        for l in 0..=l_max {
            let mut alpha_row = Vec::with_capacity(n_max);
            for value in &a {
                let ratio = 1e-3 / value.powi(l as i32);
                alpha_row.push(-ratio.ln() / (value * value));
            }
            let mut overlap = vec![vec![0.0; n_max]; n_max];
            let gamma = gamma_half(l);
            for i in 0..n_max {
                for j in 0..n_max {
                    overlap[i][j] =
                        0.5 * gamma * (alpha_row[i] + alpha_row[j]).powf(-(l as f64) - 1.5);
                }
            }
            betas.push(inverse_sqrt(&overlap));
            alphas.push(alpha_row);
        }
        let eta = 1.0 / (2.0 * options.sigma * options.sigma);
        Self {
            n_max,
            l_max,
            one_over_eta: 1.0 / eta,
            alphas,
            betas,
        }
    }

    /// Adds the centre atom itself at `r = 0`.
    ///
    /// DScribe includes the centre in the neighbour list. Only the `l = 0`
    /// block is non-zero because every real harmonic of degree one or more
    /// vanishes at the origin.
    fn accumulate_self(&self, coefficients: &mut [f64]) {
        let degree = (self.l_max + 1).pow(2);
        let mut radial = vec![0.0; self.n_max];
        for n in 0..self.n_max {
            for k in 0..self.n_max {
                let one_over = 1.0 / (1.0 + self.one_over_eta * self.alphas[0][k]);
                radial[n] += self.betas[0][n][k] * self.one_over_eta.powf(1.5) * one_over.powf(1.5);
            }
        }
        for n in 0..self.n_max {
            coefficients[n * degree] += radial[n] * 1.570_796_326_794_896_6;
        }
    }

    fn accumulate(&self, unit: [f64; 3], distance: f64, coefficients: &mut [f64]) {
        let r2 = distance * distance;
        let x = unit[0] * distance;
        let y = unit[1] * distance;
        let z = unit[2] * distance;
        let degree = (self.l_max + 1).pow(2);
        for l in 0..=self.l_max {
            let mut exp_terms = Vec::with_capacity(self.n_max);
            let mut aoa = Vec::with_capacity(self.n_max);
            for k in 0..self.n_max {
                let alpha = self.alphas[l][k];
                let one_over = 1.0 / (1.0 + self.one_over_eta * alpha);
                aoa.push(-alpha * one_over);
                exp_terms.push(
                    self.one_over_eta.powf(1.5)
                        * one_over.powf(l as f64 + 1.5)
                        * (aoa[k] * r2).exp(),
                );
            }
            let harmonics = dscribe_harmonics(l, x, y, z, r2);
            for (m, harmonic) in harmonics.iter().enumerate() {
                for n in 0..self.n_max {
                    let mut radial = 0.0;
                    for k in 0..self.n_max {
                        radial += self.betas[l][n][k] * exp_terms[k];
                    }
                    coefficients[n * degree + l * l + m] += radial * harmonic;
                }
            }
        }
    }
}

fn gamma_half(l: usize) -> f64 {
    // Gamma(l + 3/2) for integer l.
    let mut value = std::f64::consts::PI.sqrt() / 2.0; // Gamma(3/2)
    for index in 0..l {
        value *= index as f64 + 1.5;
    }
    value
}

fn inverse_sqrt(matrix: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = matrix.len();
    let (values, vectors) = jacobi(matrix);
    let mut inverse = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            let mut sum = 0.0;
            for k in 0..n {
                if values[k] > 1e-14 {
                    sum += vectors[i][k] * vectors[j][k] / values[k].sqrt();
                }
            }
            inverse[i][j] = sum;
        }
    }
    inverse
}

fn jacobi(matrix: &[Vec<f64>]) -> (Vec<f64>, Vec<Vec<f64>>) {
    let n = matrix.len();
    let mut a = matrix.to_vec();
    let mut v = vec![vec![0.0; n]; n];
    for (i, row) in v.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    for _ in 0..100 {
        let mut off = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                off += a[i][j] * a[i][j];
            }
        }
        if off < 1e-24 {
            break;
        }
        for p in 0..n {
            for q in (p + 1)..n {
                if a[p][q].abs() < 1e-18 {
                    continue;
                }
                let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                let app = a[p][p];
                let aqq = a[q][q];
                let apq = a[p][q];
                a[p][p] = c * c * app - 2.0 * s * c * apq + s * s * aqq;
                a[q][q] = s * s * app + 2.0 * s * c * apq + c * c * aqq;
                a[p][q] = 0.0;
                a[q][p] = 0.0;
                for r in 0..n {
                    if r != p && r != q {
                        let arp = a[r][p];
                        let arq = a[r][q];
                        a[r][p] = c * arp - s * arq;
                        a[p][r] = a[r][p];
                        a[r][q] = s * arp + c * arq;
                        a[q][r] = a[r][q];
                    }
                }
                for r in 0..n {
                    let vrp = v[r][p];
                    let vrq = v[r][q];
                    v[r][p] = c * vrp - s * vrq;
                    v[r][q] = s * vrp + c * vrq;
                }
            }
        }
    }
    ((0..n).map(|index| a[index][index]).collect(), v)
}

/// DScribe real-harmonic polynomials through degree five.
#[allow(clippy::many_single_char_names)] // x, y, z are the coordinates
fn dscribe_harmonics(l: usize, x: f64, y: f64, z: f64, r2: f64) -> Vec<f64> {
    let x2 = x * x;
    let y2 = y * y;
    let z2 = z * z;
    let r4 = r2 * r2;
    match l {
        0 => vec![1.570_796_326_794_896_6],
        1 => vec![
            2.720_699_046_384_954_3 * z,
            2.720_699_046_384_954_3 * x,
            2.720_699_046_384_954_3 * y,
        ],
        2 => vec![
            1.092_548_430_592_08 * x * y,
            1.092_548_430_592_08 * y * z,
            -0.315_391_565_252_52 * x2 - 0.315_391_565_252_52 * y2 + 0.630_783_130_505_04 * z2,
            1.092_548_430_592_08 * x * z,
            0.546_274_215_296_04 * x2 - 0.546_274_215_296_04 * y2,
        ],
        3 => vec![
            0.590_043_589_926_644 * y * (3.0 * x2 - y2),
            2.890_611_442_640_55 * x * y * z,
            -0.457_045_799_464_466 * y * (x2 + y2 - 4.0 * z2),
            0.373_176_332_590_115 * z * (-3.0 * x2 - 3.0 * y2 + 2.0 * z2),
            -0.457_045_799_464_466 * x * (x2 + y2 - 4.0 * z2),
            1.445_305_721_320_28 * z * (x2 - y2),
            0.590_043_589_926_644 * x * (x2 - 3.0 * y2),
        ],
        4 => vec![
            2.503_342_941_796_7 * x * y * (x2 - y2),
            1.770_130_769_779_93 * y * z * (3.0 * x2 - y2),
            -0.946_174_695_757_56 * x * y * (x2 + y2 - 6.0 * z2),
            0.669_046_543_557_289 * y * z * (-3.0 * x2 - 3.0 * y2 + 4.0 * z2),
            3.702_494_142_032_15 * z2 * z2 - 3.173_566_407_456_13 * z2 * r2
                + 0.317_356_640_745_613 * r4,
            0.669_046_543_557_289 * x * z * (-3.0 * x2 - 3.0 * y2 + 4.0 * z2),
            -0.473_087_347_878_78 * (x2 - y2) * (x2 + y2 - 6.0 * z2),
            1.770_130_769_779_93 * x * z * (x2 - 3.0 * y2),
            0.625_835_735_449_176 * x2 * x2 - 3.755_014_412_695_06 * x2 * y2
                + 0.625_835_735_449_176 * y2 * y2,
        ],
        _ => vec![
            0.656_382_056_840_17 * y * (5.0 * x2 * x2 - 10.0 * x2 * y2 + y2 * y2),
            8.302_649_259_524_16 * x * y * z * (x2 - y2),
            -0.489_238_299_435_25 * y * (3.0 * x2 - y2) * (x2 + y2 - 8.0 * z2),
            4.793_536_784_973_32 * x * y * z * (-x2 - y2 + 2.0 * z2),
            0.452_946_651_195_697 * y * (21.0 * z2 * z2 - 14.0 * z2 * r2 + r4),
            0.116_950_322_453_424 * z * (63.0 * z2 * z2 - 70.0 * z2 * r2 + 15.0 * r4),
            0.452_946_651_195_697 * x * (21.0 * z2 * z2 - 14.0 * z2 * r2 + r4),
            2.396_768_392_486_66 * z * (x2 - y2) * (-x2 - y2 + 2.0 * z2),
            -0.489_238_299_435_25 * x * (x2 - 3.0 * y2) * (x2 + y2 - 8.0 * z2),
            2.075_662_314_881_04 * z * (x2 * x2 - 6.0 * x2 * y2 + y2 * y2),
            0.656_382_056_840_17 * x * (x2 * x2 - 10.0 * x2 * y2 + 5.0 * y2 * y2),
        ],
    }
}

/// The number of features for one centre.
#[must_use]
pub fn soap_feature_count(n_max: usize, l_max: usize) -> usize {
    (n_max * (n_max + 1) / 2) * (l_max + 1)
}
