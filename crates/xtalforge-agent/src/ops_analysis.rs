//! Analysis verbs: symmetry, bonds, rings, PXRD and porosity.

use serde_json::Value;
use xtalforge_analysis::{PxrdOptions, Radiation, calculate_pxrd, default_radii, distance_grid};
use xtalforge_symmetry::{SymmetryOptions, detect_symmetry};

use crate::answers::{
    Answer, BondSummary, PatternSample, PorositySummary, RingCount, SymmetrySummary,
};
use crate::run::{document_argument, number_or, optional_integer, optional_string};
use crate::{Error, Session};

pub(crate) fn symmetry(session: &mut Session, arguments: &Value) -> Result<Answer, Error> {
    let name = session.resolve(document_argument(arguments)?)?;
    let symprec = number_or(arguments, "symprec", 1e-5)?;
    let document = session.document(Some(&name))?;
    let options = SymmetryOptions {
        symprec,
        ..Default::default()
    };
    let dataset = detect_symmetry(document.structure(), options)?;
    let n_atoms = dataset.n_atoms();
    Ok(Answer::Symmetry {
        document: name,
        symmetry: SymmetrySummary {
            number: dataset.number,
            hm_symbol: dataset.hm_symbol,
            pearson_symbol: dataset.pearson_symbol,
            n_operations: dataset.operations.len(),
            n_atoms,
            wyckoffs: dataset.wyckoffs,
            site_symmetries: dataset.site_symmetry_symbols,
        },
    })
}

pub(crate) fn bonds(session: &mut Session, arguments: &Value) -> Result<Answer, Error> {
    let name = session.resolve(document_argument(arguments)?)?;
    let document = session.document(Some(&name))?;
    let p1 = document.structure().p1()?;
    let graph = document.structure().bond_graph()?;
    let bonds = graph
        .bonds
        .iter()
        .map(|bond| BondSummary {
            i: bond.i,
            j: bond.j,
            image: bond.image,
            distance: bond.distance,
            order: bond.order,
            explicit: bond.explicit,
        })
        .collect();
    Ok(Answer::Bonds {
        document: name,
        n_atoms: p1.atoms.len(),
        bonds,
    })
}

pub(crate) fn rings(session: &mut Session, arguments: &Value) -> Result<Answer, Error> {
    let name = session.resolve(document_argument(arguments)?)?;
    let max_size = optional_integer(arguments, "max_size")?.unwrap_or(12);
    let max_size = usize::try_from(max_size).map_err(|_| Error::Argument {
        name: "max_size".to_owned(),
        reason: "does not fit in usize".to_owned(),
    })?;
    let document = session.document(Some(&name))?;
    let counts = xtalforge_core::ring_counts(document.structure(), max_size)?;
    let total: usize = counts.iter().map(|(_, count)| *count).sum();
    Ok(Answer::Rings {
        document: name,
        max_size,
        counts: counts
            .into_iter()
            .map(|(size, count)| RingCount { size, count })
            .collect(),
        total,
    })
}

fn radiation(name: &str) -> Result<Radiation, Error> {
    match name.trim().to_ascii_lowercase().as_str() {
        "cu-ka" | "cuka" => Ok(Radiation::CuKa),
        "cu-ka1" | "cuka1" => Ok(Radiation::CuKa1),
        "mo-ka" | "moka" => Ok(Radiation::MoKa),
        "co-ka" | "coka" => Ok(Radiation::CoKa),
        "fe-ka" | "feka" => Ok(Radiation::FeKa),
        "cr-ka" | "crka" => Ok(Radiation::CrKa),
        "ag-ka" | "agka" => Ok(Radiation::AgKa),
        other => Err(Error::Argument {
            name: "radiation".to_owned(),
            reason: format!("unknown radiation {other:?}"),
        }),
    }
}

pub(crate) fn pxrd(session: &mut Session, arguments: &Value) -> Result<Answer, Error> {
    let name = session.resolve(document_argument(arguments)?)?;
    let radiation = radiation(optional_string(arguments, "radiation")?.unwrap_or("cu-ka"))?;
    let minimum = number_or(arguments, "two_theta_min", 5.0)?;
    let maximum = number_or(arguments, "two_theta_max", 80.0)?;
    let fwhm = number_or(arguments, "fwhm", 0.1)?;
    let document = session.document(Some(&name))?;
    let options = PxrdOptions {
        radiation,
        two_theta_range: [minimum, maximum],
        step: 0.02,
        fwhm,
    };
    let pattern = calculate_pxrd(document.structure(), &options)?;
    let maximum_intensity = pattern
        .reflections
        .iter()
        .map(|reflection| reflection.intensity)
        .fold(0.0_f64, f64::max);
    let reflections = pattern
        .reflections
        .iter()
        .take(240)
        .map(|reflection| {
            [
                reflection.two_theta,
                if maximum_intensity > 0.0 {
                    reflection.intensity * 100.0 / maximum_intensity
                } else {
                    0.0
                },
            ]
        })
        .collect();
    let points = decimate(&pattern, 2000);
    Ok(Answer::Pxrd {
        document: name,
        radiation: radiation.label(),
        n_reflections: pattern.reflections.len(),
        reflections,
        points,
    })
}

fn decimate(pattern: &xtalforge_analysis::Pattern, limit: usize) -> Vec<PatternSample> {
    let normalised = pattern.normalised();
    let stride = normalised.len().div_ceil(limit).max(1);
    normalised
        .iter()
        .step_by(stride)
        .map(|point| PatternSample {
            two_theta: point.two_theta,
            intensity: point.intensity,
        })
        .collect()
}

pub(crate) fn porosity(session: &mut Session, arguments: &Value) -> Result<Answer, Error> {
    let name = session.resolve(document_argument(arguments)?)?;
    let spacing = number_or(arguments, "spacing", 0.4)?;
    let probe = number_or(arguments, "probe_radius", 1.2)?;
    let document = session.document(Some(&name))?;
    let structure = document.structure();
    let radii = default_radii(structure)?;
    let grid = distance_grid(structure, &radii, spacing)?;
    let accessible = grid.accessible_volume(probe)?;
    let cell_volume = grid.cell_volume();
    Ok(Answer::Porosity {
        document: name,
        porosity: PorositySummary {
            spacing,
            probe_radius: probe,
            cell_volume,
            accessible_volume: accessible,
            accessible_fraction: if cell_volume > 0.0 {
                accessible / cell_volume
            } else {
                0.0
            },
            surface_area: grid.surface_area(),
            largest_included_diameter: grid.largest_included_diameter(),
        },
    })
}
