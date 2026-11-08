//! Copy and paste as a Cartesian, cell-free fragment (Task 3.10.6).

use std::collections::BTreeSet;

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Site, Structure};

use crate::command::Command;
use crate::document::Document;
use crate::edit::rebuild;

/// Atoms lifted out of a structure, in Cartesian coordinates.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct Fragment {
    /// `(element symbol, Cartesian position)` in Angstrom.
    pub points: Vec<(String, [f64; 3])>,
}

impl Fragment {
    /// Creates a fragment from Cartesian points.
    #[must_use]
    pub fn new(points: Vec<(String, [f64; 3])>) -> Self {
        Self { points }
    }

    /// Number of atoms.
    #[must_use]
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// Whether the fragment is empty.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }

    /// Cartesian centroid, or the origin for an empty fragment.
    #[must_use]
    pub fn centroid(&self) -> [f64; 3] {
        if self.points.is_empty() {
            return [0.0; 3];
        }
        let mut sum = [0.0; 3];
        for (_, point) in &self.points {
            sum[0] += point[0];
            sum[1] += point[1];
            sum[2] += point[2];
        }
        let count = self.points.len() as f64;
        [sum[0] / count, sum[1] / count, sum[2] / count]
    }
}

/// Copies the selected sites as a fragment recentred on its centroid.
#[must_use]
pub fn copy_selection(document: &Document) -> Fragment {
    let structure = document.structure();
    let mut selected: Vec<(String, [f64; 3])> = Vec::new();
    for index in document.selection() {
        let Some(site) = structure.sites().get(*index) else {
            continue;
        };
        selected.push((
            site.element().to_owned(),
            structure.lattice().to_cart(site.frac()),
        ));
    }
    if selected.is_empty() {
        return Fragment::default();
    }
    let mut centroid = [0.0; 3];
    for (_, point) in &selected {
        centroid[0] += point[0];
        centroid[1] += point[1];
        centroid[2] += point[2];
    }
    let count = selected.len() as f64;
    centroid = [
        centroid[0] / count,
        centroid[1] / count,
        centroid[2] / count,
    ];
    let points = selected
        .into_iter()
        .map(|(element, point)| {
            (
                element,
                [
                    point[0] - centroid[0],
                    point[1] - centroid[1],
                    point[2] - centroid[2],
                ],
            )
        })
        .collect();
    Fragment { points }
}

/// Adds a fragment's atoms to the asymmetric unit at a Cartesian offset.
#[derive(Debug, Clone)]
pub struct PasteFragment {
    fragment: Fragment,
    offset: [f64; 3],
}

impl PasteFragment {
    /// Creates the command.
    #[must_use]
    pub fn new(fragment: Fragment, offset: [f64; 3]) -> Self {
        Self { fragment, offset }
    }

    /// Fragment being pasted.
    #[must_use]
    pub fn fragment(&self) -> &Fragment {
        &self.fragment
    }

    /// Cartesian offset applied to every point.
    #[must_use]
    pub fn offset(&self) -> [f64; 3] {
        self.offset
    }
}

impl Command for PasteFragment {
    fn apply(&self, structure: &Structure) -> Result<Structure, Error> {
        if !self.offset.iter().all(|value| value.is_finite()) {
            return Err(Error::InvalidInput(
                "paste offset must be finite".to_owned(),
            ));
        }
        let mut used: BTreeSet<String> = structure
            .sites()
            .iter()
            .map(|site| site.label().to_owned())
            .collect();
        let mut counter = structure.sites().len() + 1;
        let mut sites = structure.sites().to_vec();
        for (element, point) in &self.fragment.points {
            if !point.iter().all(|value| value.is_finite()) {
                return Err(Error::InvalidInput(
                    "fragment coordinates must be finite".to_owned(),
                ));
            }
            let cart = [
                point[0] + self.offset[0],
                point[1] + self.offset[1],
                point[2] + self.offset[2],
            ];
            let frac = structure.lattice().to_frac(cart);
            let label = unique_label(element, &mut used, &mut counter);
            sites.push(Site::new(label, element, frac)?);
        }
        rebuild(structure, sites, *structure.lattice())
    }

    fn description(&self) -> String {
        format!("Paste {} atoms", self.fragment.points.len())
    }
}

fn unique_label(element: &str, used: &mut BTreeSet<String>, counter: &mut usize) -> String {
    loop {
        let label = format!("{element}{counter}");
        *counter += 1;
        if used.insert(label.clone()) {
            return label;
        }
    }
}
