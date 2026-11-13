//! Machine-learning descriptors (Task 16.3.3).
//!
//! Four descriptors are provided:
//!
//! * [`composition_features`] — element fractions, counts, mean atomic number,
//!   mean mass, mean electronegativity and metal fraction.
//! * [`rdf_descriptor`] — the total and partial radial distribution function
//!   from Task 9.3.4, ready for a feature vector.
//! * [`coulomb_matrix`] — the standard molecular Coulomb matrix evaluated with
//!   minimum-image distances.
//! * [`soap_like_descriptor`] — a documented, SOAP-like local descriptor.
//!
//! Every descriptor is a `serde` value with a stable field order, so the JSON
//! output is deterministic. Descriptors work on the P1 expansion of the
//! structure: symmetry operations are applied first, occupancies are treated
//! as one, and sites are sorted by atomic number, then fractional position,
//! then label. The canonical order makes the output independent of the input
//! site order.
//!
//! # SOAP-like descriptor
//!
//! The local descriptor is a rotation-invariant expansion of the neighbour
//! density around each centre. For radial basis functions `g_n(r)` and
//! Legendre polynomials `P_l(cos theta)`, the feature blocks are
//!
//! ```text
//! radial(s, n)               = sum_{j in s} g_n(r_j) / max(1, N_s)
//! angular(s, t, n, n', l)    = sum_{j in s, k in t, j != k}
//!                              g_n(r_j) g_n'(r_k) P_l(cos theta_jk)
//!                              / max(1, N_s N_t)
//! ```
//!
//! with `g_n(r) = exp(-0.5 ((r - r_n) / sigma)^2)`. This is the power-spectrum
//! idea of SOAP without the spherical-harmonic projection: it is invariant to
//! rotations and neighbour permutations, uses a Gaussian radial basis and
//! Legendre angular moments. It is **not** the exact SOAP power spectrum, so
//! values are not comparable with a full SOAP implementation.

mod composition;
mod coulomb;
mod rdf;
mod soap;

pub use composition::{CompositionFeatures, composition_features};
pub use coulomb::{CoulombMatrix, CoulombMatrixOptions, coulomb_matrix};
pub use rdf::{RdfDescriptor, rdf_descriptor};
pub use soap::{
    MAX_SOAP_FEATURES, MAX_SOAP_L, MAX_SOAP_NEIGHBORS, MAX_SOAP_RADIAL, SoapLikeDescriptor,
    SoapLikeOptions, soap_like_descriptor,
};

use xtalforge_core::{Error, Lattice, Structure};

/// One atom in the canonical descriptor order.
pub(crate) struct CanonicalSite {
    pub(crate) element: String,
    pub(crate) z: u8,
    pub(crate) mass: f64,
    pub(crate) electronegativity: Option<f64>,
    pub(crate) is_metal: bool,
    pub(crate) frac: [f64; 3],
    pub(crate) label: String,
}

/// Expands the structure to P1 and sorts its atoms deterministically.
pub(crate) fn canonical_cell(
    structure: &Structure,
) -> Result<(Lattice, Vec<CanonicalSite>), Error> {
    let p1 = structure.p1()?;
    if p1.atoms.is_empty() {
        return Err(Error::InvalidInput(
            "a descriptor needs at least one atom".to_owned(),
        ));
    }
    let mut sites = Vec::with_capacity(p1.atoms.len());
    for atom in &p1.atoms {
        let info = xtalforge_core::elements::element_by_symbol(&atom.element)?;
        sites.push(CanonicalSite {
            element: info.symbol.clone(),
            z: info.z,
            mass: info.mass,
            electronegativity: info.electronegativity,
            is_metal: info.is_metal,
            frac: atom.frac,
            label: atom.label.clone(),
        });
    }
    sites.sort_by(|left, right| {
        left.z
            .cmp(&right.z)
            .then_with(|| left.frac[0].total_cmp(&right.frac[0]))
            .then_with(|| left.frac[1].total_cmp(&right.frac[1]))
            .then_with(|| left.frac[2].total_cmp(&right.frac[2]))
            .then_with(|| left.label.cmp(&right.label))
    });
    Ok((p1.lattice, sites))
}

/// Builds a P1 structure in the canonical order for the neighbour search.
pub(crate) fn canonical_snapshot(
    lattice: Lattice,
    sites: &[CanonicalSite],
) -> Result<Structure, Error> {
    let mut builder = xtalforge_core::StructureBuilder::new("descriptor", lattice);
    for site in sites {
        builder = builder.site(xtalforge_core::Site::new(
            site.label.clone(),
            &site.element,
            site.frac,
        )?);
    }
    builder.build()
}

/// Serializes a descriptor as pretty, deterministic JSON.
pub(crate) fn to_json<T: serde::Serialize>(value: &T, label: &str) -> Result<String, Error> {
    serde_json::to_string_pretty(value)
        .map_err(|error| Error::InvalidInput(format!("{label} encoding: {error}")))
}
