//! Scene data any viewer can consume (Task 3.14.4).

use serde::{Deserialize, Serialize};
use xtalforge_core::elements::element_by_symbol;
use xtalforge_core::{Error, Structure};

use crate::{DrawStyle, ViewSettings};

/// One atom ready to draw.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SceneAtom {
    /// Element symbol.
    pub element: String,
    /// Cartesian position in Angstrom.
    pub position: [f64; 3],
    /// Display radius in Angstrom.
    pub radius: f64,
    /// 0-255 RGB.
    pub color: [u8; 3],
    /// Index in the P1 cell.
    pub index: usize,
    /// Hidden by the current settings.
    pub hidden: bool,
}

/// One bond ready to draw.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SceneBond {
    /// First atom index.
    pub i: usize,
    /// Second atom index.
    pub j: usize,
    /// Cartesian start.
    pub start: [f64; 3],
    /// Cartesian end.
    pub end: [f64; 3],
    /// Display radius in Angstrom.
    pub radius: f64,
    /// 0-255 RGB.
    pub color: [u8; 3],
    /// Bond order.
    pub order: f64,
}

/// Everything a viewer needs for one frame.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SceneData {
    /// Atom spheres.
    pub atoms: Vec<SceneAtom>,
    /// Bonds.
    pub bonds: Vec<SceneBond>,
    /// Cell edge segments as Cartesian pairs.
    pub cell_edges: Vec<([f64; 3], [f64; 3])>,
    /// Settings used to build the scene.
    pub settings: ViewSettings,
}

/// Builds a scene from a structure, its P1 bonds and the view settings.
///
/// # Errors
///
/// Returns the errors of the P1 expansion or bond perception.
pub fn build_scene(structure: &Structure, settings: &ViewSettings) -> Result<SceneData, Error> {
    let mut settings = settings.clone();
    settings.sanitize();
    let p1 = structure.p1()?;
    let graph = structure.bond_graph()?;
    let mut atoms = Vec::with_capacity(p1.atoms.len());
    for (index, atom) in p1.atoms.iter().enumerate() {
        let info = element_by_symbol(&atom.element)?;
        let radius = match settings.style {
            DrawStyle::SpaceFill => info.vdw_radius * settings.atom_scale,
            DrawStyle::Wireframe | DrawStyle::Sticks => 0.0,
            DrawStyle::BallStick => 0.35 * info.vdw_radius * settings.atom_scale,
        };
        atoms.push(SceneAtom {
            element: atom.element.clone(),
            position: p1.lattice.to_cart(atom.frac),
            radius,
            color: match settings.style {
                DrawStyle::Wireframe => [0, 0, 0],
                _ => info.color,
            },
            index,
            hidden: settings.hide_hydrogens && matches!(atom.element.as_str(), "H" | "D"),
        });
    }
    let mut bonds = Vec::with_capacity(graph.bonds.len());
    for bond in &graph.bonds {
        let Some(first) = atoms.get(bond.i) else {
            continue;
        };
        let Some(second) = atoms.get(bond.j) else {
            continue;
        };
        let image = p1.lattice.to_cart(bond.image.map(f64::from));
        let start = first.position;
        let end = [
            second.position[0] + image[0],
            second.position[1] + image[1],
            second.position[2] + image[2],
        ];
        bonds.push(SceneBond {
            i: bond.i,
            j: bond.j,
            start,
            end,
            radius: if settings.style == DrawStyle::Wireframe {
                0.02
            } else {
                settings.bond_radius
            },
            color: [90, 90, 90],
            order: bond.order,
        });
    }
    let matrix = structure.lattice().matrix();
    let mut cell_edges = Vec::new();
    if settings.show_cell {
        let corners = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 1.0, 0.0],
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
            [1.0, 1.0, 1.0],
        ];
        let cart = corners.map(|point| {
            [
                point[0] * matrix[0][0] + point[1] * matrix[1][0] + point[2] * matrix[2][0],
                point[0] * matrix[0][1] + point[1] * matrix[1][1] + point[2] * matrix[2][1],
                point[0] * matrix[0][2] + point[1] * matrix[1][2] + point[2] * matrix[2][2],
            ]
        });
        for (first, second) in [
            (0, 1),
            (0, 2),
            (0, 3),
            (1, 4),
            (1, 5),
            (2, 4),
            (2, 6),
            (3, 5),
            (3, 6),
            (4, 7),
            (5, 7),
            (6, 7),
        ] {
            if let (Some(first), Some(second)) = (cart.get(first), cart.get(second)) {
                cell_edges.push((*first, *second));
            }
        }
    }
    Ok(SceneData {
        atoms,
        bonds,
        cell_edges,
        settings,
    })
}
