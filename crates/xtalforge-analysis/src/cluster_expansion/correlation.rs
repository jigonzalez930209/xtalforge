//! Correlation evaluation of one cluster-model assignment.

use xtalforge_core::Error;

use super::{ClusterModel, ClusterShellData};

/// Correlations of one assignment, one value per shell.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ClusterCorrelations {
    /// Pair correlations.
    pub(crate) pair: Vec<f64>,
    /// Triplet correlations.
    pub(crate) triplet: Vec<f64>,
    /// Quadruplet correlations.
    pub(crate) quadruplet: Vec<f64>,
}

/// Maps species labels of one decoration to indices.
pub(crate) fn species_indices(
    species: &[String],
    decoration: &[String],
) -> Result<Vec<usize>, Error> {
    let mut indices = Vec::with_capacity(decoration.len());
    for label in decoration {
        let index = species.iter().position(|candidate| candidate == label);
        let Some(index) = index else {
            return Err(Error::InvalidInput(format!(
                "the decoration contains unknown species {label}"
            )));
        };
        indices.push(index);
    }
    Ok(indices)
}

/// Evaluates the pair, triplet and quadruplet correlations of one assignment.
pub(crate) fn correlations_for_assignment(
    model: &ClusterModel,
    assignment: &[usize],
) -> Result<ClusterCorrelations, Error> {
    if assignment.len() != model.n_sites {
        return Err(Error::InvalidInput(format!(
            "the assignment needs {} sites, got {}",
            model.n_sites,
            assignment.len()
        )));
    }
    let species_count = model.fractions.len();
    if species_count == 0 {
        return Err(Error::InvalidInput(
            "the cluster model has no species".to_owned(),
        ));
    }
    let mut centered = vec![vec![0.0; species_count]; model.n_sites];
    for (site, &species) in assignment.iter().enumerate() {
        if species >= species_count {
            return Err(Error::InvalidInput(format!(
                "species index {species} is outside the model"
            )));
        }
        let Some(row) = centered.get_mut(site) else {
            return Err(Error::InvalidInput(
                "the assignment index is outside the model".to_owned(),
            ));
        };
        for (index, fraction) in model.fractions.iter().enumerate() {
            if let Some(value) = row.get_mut(index) {
                *value = if index == species {
                    1.0 - fraction
                } else {
                    -fraction
                };
            }
        }
    }
    Ok(ClusterCorrelations {
        pair: shell_correlations(&model.pair, &centered),
        triplet: shell_correlations(&model.triplet, &centered),
        quadruplet: shell_correlations(&model.quadruplet, &centered),
    })
}

pub(crate) fn feature_names(model: &ClusterModel) -> Vec<String> {
    let mut names =
        Vec::with_capacity(model.pair.len() + model.triplet.len() + model.quadruplet.len());
    for index in 0..model.pair.len() {
        names.push(format!("pair[{index}]"));
    }
    for index in 0..model.triplet.len() {
        names.push(format!("triplet[{index}]"));
    }
    for index in 0..model.quadruplet.len() {
        names.push(format!("quadruplet[{index}]"));
    }
    names
}

fn shell_correlations(shells: &[ClusterShellData], centered: &[Vec<f64>]) -> Vec<f64> {
    let species_count = centered.first().map_or(0, Vec::len);
    shells
        .iter()
        .map(|shell| {
            if shell.clusters.is_empty() {
                return 0.0;
            }
            let mut total = 0.0;
            for cluster in &shell.clusters {
                for species in 0..species_count {
                    let mut product = 1.0;
                    for &site in cluster {
                        let value = centered
                            .get(site)
                            .and_then(|row| row.get(species))
                            .copied()
                            .unwrap_or(0.0);
                        product *= value;
                    }
                    total += product;
                }
            }
            total / shell.clusters.len() as f64
        })
        .collect()
}
