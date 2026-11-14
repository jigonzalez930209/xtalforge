//! Band structures and band gaps on a k-path.

use serde::{Deserialize, Serialize};
use xtalforge_core::Error;

use super::ElectronicStructure;
use super::validate::{self, invalid};

/// Fraction of an electron that marks a state as occupied in a gap analysis.
///
/// Smearing spreads an occupation over neighbouring states; the threshold at
/// half filling keeps the classification stable for standard semiconductors.
const OCCUPIED_THRESHOLD: f64 = 0.5;

/// Spin channel of a band structure.
#[derive(
    Debug, Clone, Copy, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize,
)]
#[serde(rename_all = "snake_case")]
pub enum Spin {
    /// Spin up, or the only channel of an unpolarized calculation.
    #[default]
    Up,
    /// Spin down.
    Down,
}

impl Spin {
    /// Stable label used in reports and serialized output.
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Up => "up",
            Self::Down => "down",
        }
    }
}

/// High-symmetry label attached to one k-point index.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct KPointLabel {
    /// Index into the k-point list.
    pub index: usize,
    /// Label such as `Gamma`, `X` or `L`.
    pub label: String,
}

/// One band along the path.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Band {
    /// Energies in eV, one per k-point.
    pub energies: Vec<f64>,
}

/// Band energies on a high-symmetry k-path for one spin channel.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BandStructure {
    /// Spin channel.
    pub spin: Spin,
    /// High-symmetry labels.
    #[serde(default)]
    pub labels: Vec<KPointLabel>,
    /// Fractional k-point coordinates.
    pub kpoints: Vec<[f64; 3]>,
    /// Path coordinate for plotting, in the unit supplied by the input.
    pub distances: Vec<f64>,
    /// Bands in file order.
    pub bands: Vec<Band>,
    /// Occupation numbers with the same `[band][kpoint]` shape, when reported.
    #[serde(default)]
    pub occupations: Option<Vec<Vec<f64>>>,
    /// Fermi or reference energy in eV.
    #[serde(default)]
    pub fermi: Option<f64>,
}

impl BandStructure {
    /// Number of bands on the path.
    #[must_use]
    pub fn band_count(&self) -> usize {
        self.bands.len()
    }

    /// Number of k-points on the path.
    #[must_use]
    pub fn kpoint_count(&self) -> usize {
        self.kpoints.len()
    }

    /// Label of one k-point, when the path carries one.
    #[must_use]
    pub fn label_at(&self, index: usize) -> Option<&str> {
        self.labels
            .iter()
            .find(|label| label.index == index)
            .map(|label| label.label.as_str())
    }

    /// Energy of one state in eV.
    #[must_use]
    pub fn energy(&self, band: usize, kpoint: usize) -> Option<f64> {
        let band = self.bands.get(band)?;
        band.energies.get(kpoint).copied()
    }

    /// Checks shapes, finiteness and label indices.
    ///
    /// This is public because the fields are public: a band structure built by
    /// hand must satisfy the same invariants as one returned by
    /// [`band_structure`].
    pub fn validate(&self) -> Result<(), Error> {
        if self.kpoints.is_empty() {
            return Err(invalid("a band structure needs at least one k-point"));
        }
        validate::kpoints(&self.kpoints)?;
        validate::distances(&self.distances, self.kpoints.len(), "distances")?;
        validate::bands(&self.bands, self.kpoints.len(), "bands")?;
        if let Some(occupations) = &self.occupations {
            validate::occupations_against_bands(occupations, &self.bands, "occupations")?;
        }
        validate::labels(&self.labels, self.kpoints.len())?;
        validate::fermi(self.fermi)?;
        Ok(())
    }

    /// Direct and indirect gap of this channel.
    ///
    /// A state is occupied when its occupation exceeds one half. Without
    /// occupations, states below the Fermi energy are occupied. The method
    /// returns an error when neither classification is available or when the
    /// channel has no occupied or no empty state.
    pub fn gap(&self) -> Result<BandGap, Error> {
        self.validate()?;
        if self.occupations.is_none() && self.fermi.is_none() {
            return Err(invalid(
                "a band gap needs occupations or a Fermi energy to classify states",
            ));
        }
        let mut vbm = f64::NEG_INFINITY;
        let mut cbm = f64::INFINITY;
        let mut vbm_kpoint = 0_usize;
        let mut cbm_kpoint = 0_usize;
        let mut direct_energy = f64::INFINITY;
        let mut direct_points = 0_usize;
        for kpoint in 0..self.kpoints.len() {
            let mut valence = f64::NEG_INFINITY;
            let mut conduction = f64::INFINITY;
            for (band_index, band) in self.bands.iter().enumerate() {
                let energy = band.energies.get(kpoint).copied().ok_or_else(|| {
                    invalid(format!("band {band_index} is missing k-point {kpoint}"))
                })?;
                let occupied = if let Some(occupations) = &self.occupations {
                    let value = occupations
                        .get(band_index)
                        .and_then(|row| row.get(kpoint))
                        .copied()
                        .ok_or_else(|| {
                            invalid(format!(
                                "occupations are missing band {band_index}, k-point {kpoint}"
                            ))
                        })?;
                    value > OCCUPIED_THRESHOLD
                } else if let Some(fermi) = self.fermi {
                    energy < fermi
                } else {
                    return Err(invalid(
                        "a band gap needs occupations or a Fermi energy to classify states",
                    ));
                };
                if occupied {
                    valence = valence.max(energy);
                } else {
                    conduction = conduction.min(energy);
                }
            }
            if valence > vbm {
                vbm = valence;
                vbm_kpoint = kpoint;
            }
            if conduction < cbm {
                cbm = conduction;
                cbm_kpoint = kpoint;
            }
            if valence.is_finite() && conduction.is_finite() {
                direct_energy = direct_energy.min((conduction - valence).max(0.0));
                direct_points += 1;
            }
        }
        if !vbm.is_finite() {
            return Err(invalid("the band structure has no occupied state"));
        }
        if !cbm.is_finite() {
            return Err(invalid("the band structure has no empty state"));
        }
        if direct_points == 0 {
            return Err(invalid(
                "no k-point has both an occupied and an empty state",
            ));
        }
        let energy = (cbm - vbm).max(0.0);
        let metallic = cbm < vbm;
        Ok(BandGap {
            energy,
            direct: !metallic && direct_energy <= energy + 1e-12,
            direct_energy,
            vbm,
            cbm,
            vbm_kpoint,
            cbm_kpoint,
            vbm_label: self.label_at(vbm_kpoint).map(str::to_owned),
            cbm_label: self.label_at(cbm_kpoint).map(str::to_owned),
            metallic,
        })
    }
}

/// Direct and indirect band gap of one spin channel.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BandGap {
    /// Fundamental (indirect) gap in eV, clamped at zero for an overlap.
    pub energy: f64,
    /// Whether the fundamental gap is direct: the smallest vertical gap
    /// equals `energy` and the bands do not overlap.
    pub direct: bool,
    /// Smallest vertical gap in eV, clamped at zero.
    pub direct_energy: f64,
    /// Valence-band maximum in eV.
    pub vbm: f64,
    /// Conduction-band minimum in eV.
    pub cbm: f64,
    /// K-point index of the valence-band maximum.
    pub vbm_kpoint: usize,
    /// K-point index of the conduction-band minimum.
    pub cbm_kpoint: usize,
    /// Label of the valence-band maximum, when the path carries one.
    pub vbm_label: Option<String>,
    /// Label of the conduction-band minimum, when the path carries one.
    pub cbm_label: Option<String>,
    /// True when the bands overlap and the fundamental gap is zero.
    pub metallic: bool,
}

/// Gap of one spin channel.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SpinGap {
    /// Spin channel.
    pub spin: Spin,
    /// Gap of that channel.
    pub gap: BandGap,
}

/// Cumulative Euclidean distance between fractional k-points.
///
/// This is a plotting coordinate. Effective masses need Cartesian distances in
/// inverse Angstrom, which the caller supplies through
/// [`ElectronicStructure::distances`](super::ElectronicStructure::distances).
#[must_use]
pub fn path_distances(kpoints: &[[f64; 3]]) -> Vec<f64> {
    let mut distances = Vec::with_capacity(kpoints.len());
    let mut previous: Option<[f64; 3]> = None;
    let mut total = 0.0;
    for point in kpoints {
        if let Some(previous) = previous {
            let [x, y, z] = *point;
            let [px, py, pz] = previous;
            total += ((x - px).powi(2) + (y - py).powi(2) + (z - pz).powi(2)).sqrt();
        }
        distances.push(total);
        previous = Some(*point);
    }
    distances
}

/// Builds the spin-up (or unpolarized) band structure of an input.
///
/// Returns [`Error::InvalidInput`] when the input fails
/// [`ElectronicStructure::validate`].
pub fn band_structure(input: &ElectronicStructure) -> Result<BandStructure, Error> {
    input.validate()?;
    build_channel(
        input,
        Spin::Up,
        &input.energies,
        input.occupations.as_deref(),
    )
}

/// Builds one band structure per spin channel, spin up first.
pub fn spin_band_structures(input: &ElectronicStructure) -> Result<Vec<BandStructure>, Error> {
    input.validate()?;
    let mut structures = Vec::with_capacity(input.spin_count());
    structures.push(build_channel(
        input,
        Spin::Up,
        &input.energies,
        input.occupations.as_deref(),
    )?);
    if let Some(spin_down) = &input.spin_down {
        structures.push(build_channel(
            input,
            Spin::Down,
            spin_down,
            input.spin_down_occupations.as_deref(),
        )?);
    }
    Ok(structures)
}

/// Builds the band structure of one channel.
fn build_channel(
    input: &ElectronicStructure,
    spin: Spin,
    energies: &[Vec<f64>],
    occupations: Option<&[Vec<f64>]>,
) -> Result<BandStructure, Error> {
    let distances = match &input.distances {
        Some(distances) => distances.clone(),
        None => path_distances(&input.kpoints),
    };
    let bands = energies
        .iter()
        .map(|row| Band {
            energies: row.clone(),
        })
        .collect();
    Ok(BandStructure {
        spin,
        labels: input.labels.clone(),
        kpoints: input.kpoints.clone(),
        distances,
        bands,
        occupations: occupations.map(<[Vec<f64>]>::to_vec),
        fermi: input.fermi,
    })
}

/// Direct and indirect gap of the first spin channel.
pub fn band_gap(input: &ElectronicStructure) -> Result<BandGap, Error> {
    band_structure(input)?.gap()
}

/// Gaps of every spin channel, spin up first.
pub fn spin_gaps(input: &ElectronicStructure) -> Result<Vec<SpinGap>, Error> {
    spin_band_structures(input)?
        .into_iter()
        .map(|structure| {
            let spin = structure.spin;
            Ok(SpinGap {
                spin,
                gap: structure.gap()?,
            })
        })
        .collect()
}
