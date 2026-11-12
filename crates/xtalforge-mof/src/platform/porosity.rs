//! Native pore descriptors: PLD, LCD, accessible surface and volume
//! (Tasks 10.3.1-10.3.3).

use serde::{Deserialize, Serialize};
use xtalforge_analysis::{DistanceGrid, PoreSizeBin, default_radii, distance_grid};
use xtalforge_core::{Error, Structure};

/// Options for [`pore_descriptors`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PoreOptions {
    /// Grid spacing in Angstrom.
    pub spacing: f64,
    /// Probe radius for accessible surface and volume.
    pub probe_radius: f64,
    /// Probe radius for the channel-dimension test.
    pub channel_probe: f64,
}

impl Default for PoreOptions {
    fn default() -> Self {
        Self {
            spacing: 0.5,
            probe_radius: 1.2,
            channel_probe: 1.2,
        }
    }
}

/// Standard geometric descriptors of a porous framework.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PoreDescriptors {
    /// Pore-limiting diameter in Angstrom.
    pub pld: f64,
    /// Largest cavity diameter in Angstrom.
    pub lcd: f64,
    /// Accessible surface area in m2/g.
    pub asa_m2_g: f64,
    /// Accessible surface area in m2/cm3.
    pub asa_m2_cm3: f64,
    /// Accessible volume in cubic Angstrom.
    pub accessible_volume: f64,
    /// Void fraction.
    pub void_fraction: f64,
    /// Number of percolating directions at `channel_probe` (0 to 3).
    pub channel_dimension: usize,
    /// Density in g/cm3.
    pub density: f64,
    /// Probe radius used for the accessible quantities.
    pub probe_radius: f64,
    /// Grid spacing used.
    pub spacing: f64,
}

/// A pore descriptor bundle with the geometric pore-size distribution.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PoreReport {
    /// Descriptor table.
    pub descriptors: PoreDescriptors,
    /// Geometric pore-size histogram.
    pub distribution: Vec<PoreSizeBin>,
}

/// Computes native PLD, LCD, ASA, accessible volume, void fraction and
/// channel dimensionality from a distance grid.
///
/// The PLD is the largest probe radius at which a connected free component
/// spans the cell in at least one direction (periodic union-find). The LCD is
/// the largest included sphere. The accessible area is the probe isosurface.
///
/// # Errors
///
/// Returns the errors of the distance grid or a non-positive option.
pub fn pore_descriptors(structure: &Structure, options: &PoreOptions) -> Result<PoreReport, Error> {
    if !options.probe_radius.is_finite() || options.probe_radius < 0.0 {
        return Err(Error::InvalidInput(
            "the probe radius must be non-negative and finite".to_owned(),
        ));
    }
    if !options.channel_probe.is_finite() || options.channel_probe < 0.0 {
        return Err(Error::InvalidInput(
            "the channel probe must be non-negative and finite".to_owned(),
        ));
    }
    let p1 = structure.p1()?;
    let snapshot = crate::platform::util::p1_structure(&p1)?;
    let radii = default_radii(&snapshot)?;
    let grid = distance_grid(&snapshot, &radii, options.spacing)?;
    let cell_volume = grid.cell_volume();
    let accessible_volume = grid.accessible_volume(options.probe_radius)?;
    let void_fraction = accessible_volume / cell_volume;
    let area = accessible_surface_area(&grid, options.probe_radius);
    let mass = snapshot.cell_mass()?;
    let asa_m2_g = if mass > 0.0 {
        area / (mass * 1.660_539_066_60) * 1.0e4
    } else {
        0.0
    };
    let asa_m2_cm3 = area / cell_volume * 1.0e4;
    let pld = 2.0 * percolation_radius(&grid, None);
    let lcd = grid.largest_included_diameter();
    let channel_dimension = percolating_axes(&grid, options.channel_probe);
    let distribution = grid.pore_diameter(options.spacing.max(0.1))?;
    Ok(PoreReport {
        descriptors: PoreDescriptors {
            pld,
            lcd,
            asa_m2_g,
            asa_m2_cm3,
            accessible_volume,
            void_fraction,
            channel_dimension,
            density: snapshot.density()?,
            probe_radius: options.probe_radius,
            spacing: options.spacing,
        },
        distribution,
    })
}

/// Accessible surface area for a probe radius in square Angstrom.
#[must_use]
pub fn accessible_surface_area(grid: &DistanceGrid, probe_radius: f64) -> f64 {
    let [n_x, n_y, n_z] = grid.grid.shape;
    let step = grid.grid.step;
    let area_x = cross_length(step[1], step[2]);
    let area_y = cross_length(step[0], step[2]);
    let area_z = cross_length(step[0], step[1]);
    let value = |x: usize, y: usize, z: usize| {
        grid.grid
            .get(x % n_x, y % n_y, z % n_z)
            .unwrap_or(f64::NEG_INFINITY)
    };
    let mut area = 0.0;
    for z in 0..n_z {
        for y in 0..n_y {
            for x in 0..n_x {
                let inside = value(x, y, z) >= probe_radius;
                if inside != (value(x + 1, y, z) >= probe_radius) {
                    area += area_x;
                }
                if inside != (value(x, y + 1, z) >= probe_radius) {
                    area += area_y;
                }
                if inside != (value(x, y, z + 1) >= probe_radius) {
                    area += area_z;
                }
            }
        }
    }
    area
}

/// Largest probe radius that percolates along one axis.
///
/// When `axis` is `None`, any axis counts; otherwise only that axis is
/// tested. Returns zero when nothing percolates.
#[must_use]
pub fn percolation_radius(grid: &DistanceGrid, axis: Option<usize>) -> f64 {
    let [n_x, n_y, n_z] = grid.grid.shape;
    let len = n_x * n_y * n_z;
    let mut parent: Vec<usize> = (0..len).collect();
    let mut rank = vec![0u8; len];
    let mut mask = vec![0u8; len];
    let mut active = vec![false; len];
    let index = |x: usize, y: usize, z: usize| (z * n_y + y) * n_x + x;
    let mut values: Vec<(f64, usize)> = grid
        .grid
        .data
        .iter()
        .copied()
        .enumerate()
        .map(|(index, value)| (value, index))
        .collect();
    values.sort_by(|first, second| {
        second
            .0
            .partial_cmp(&first.0)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| first.1.cmp(&second.1))
    });
    let mut position = 0;
    while position < values.len() {
        let threshold = values.get(position).map(|entry| entry.0).unwrap_or(0.0);
        while position < values.len()
            && values
                .get(position)
                .is_some_and(|entry| entry.0 == threshold)
        {
            if let Some((_, flat)) = values.get(position) {
                let opened = activate(
                    *flat,
                    n_x,
                    n_y,
                    n_z,
                    &index,
                    &mut parent,
                    &mut rank,
                    &mut mask,
                    &mut active,
                    axis,
                );
                if opened {
                    return threshold.max(0.0);
                }
            }
            position += 1;
        }
    }
    0.0
}

fn spans(bits: u8, axis: Option<usize>) -> bool {
    let pairs = [(0, 1), (2, 3), (4, 5)];
    pairs.iter().enumerate().any(|(position, (low, high))| {
        if axis.is_some_and(|axis| axis != position) {
            return false;
        }
        bits & (1 << low) != 0 && bits & (1 << high) != 0
    })
}

fn percolating_axes(grid: &DistanceGrid, probe: f64) -> usize {
    let [n_x, n_y, n_z] = grid.grid.shape;
    let len = n_x * n_y * n_z;
    let mut parent: Vec<usize> = (0..len).collect();
    let mut rank = vec![0u8; len];
    let mut mask = vec![0u8; len];
    let mut active = vec![false; len];
    let index = |x: usize, y: usize, z: usize| (z * n_y + y) * n_x + x;
    for z in 0..n_z {
        for y in 0..n_y {
            for x in 0..n_x {
                let flat = index(x, y, z);
                if grid.grid.get(x, y, z).is_some_and(|value| value >= probe) {
                    activate(
                        flat,
                        n_x,
                        n_y,
                        n_z,
                        &index,
                        &mut parent,
                        &mut rank,
                        &mut mask,
                        &mut active,
                        None,
                    );
                }
            }
        }
    }
    let mut dimensions = 0;
    for axis in 0..3 {
        let found = (0..len).any(|flat| {
            if !active.get(flat).copied().unwrap_or(false) {
                return false;
            }
            let root = find(&mut parent, flat);
            spans(mask.get(root).copied().unwrap_or(0), Some(axis))
        });
        if found {
            dimensions += 1;
        }
    }
    dimensions
}

#[allow(clippy::too_many_arguments)]
fn activate(
    flat: usize,
    n_x: usize,
    n_y: usize,
    n_z: usize,
    index: &impl Fn(usize, usize, usize) -> usize,
    parent: &mut [usize],
    rank: &mut [u8],
    mask: &mut [u8],
    active: &mut [bool],
    axis: Option<usize>,
) -> bool {
    if active.get(flat).copied().unwrap_or(true) {
        return false;
    }
    if let Some(slot) = active.get_mut(flat) {
        *slot = true;
    }
    let x = flat % n_x;
    let y = (flat / n_x) % n_y;
    let z = flat / (n_x * n_y);
    let mut bits = 0u8;
    if x == 0 {
        bits |= 1 << 0;
    }
    if x + 1 == n_x {
        bits |= 1 << 1;
    }
    if y == 0 {
        bits |= 1 << 2;
    }
    if y + 1 == n_y {
        bits |= 1 << 3;
    }
    if z == 0 {
        bits |= 1 << 4;
    }
    if z + 1 == n_z {
        bits |= 1 << 5;
    }
    if let Some(slot) = mask.get_mut(flat) {
        *slot = bits;
    }
    let neighbours = [
        index((x + 1) % n_x, y, z),
        index((x + n_x - 1) % n_x, y, z),
        index(x, (y + 1) % n_y, z),
        index(x, (y + n_y - 1) % n_y, z),
        index(x, y, (z + 1) % n_z),
        index(x, y, (z + n_z - 1) % n_z),
    ];
    let mut opened = false;
    for neighbour in neighbours {
        if active.get(neighbour).copied().unwrap_or(false)
            && let Some(root) = union(flat, neighbour, parent, rank, mask)
            && spans(mask.get(root).copied().unwrap_or(0), axis)
        {
            opened = true;
        }
    }
    opened
}

fn find(parent: &mut [usize], node: usize) -> usize {
    let mut root = node;
    while parent.get(root).copied().unwrap_or(root) != root {
        root = parent.get(root).copied().unwrap_or(root);
    }
    let mut current = node;
    while parent.get(current).copied().unwrap_or(current) != root {
        let next = parent.get(current).copied().unwrap_or(root);
        if let Some(slot) = parent.get_mut(current) {
            *slot = root;
        }
        current = next;
    }
    root
}

fn union(
    first: usize,
    second: usize,
    parent: &mut [usize],
    rank: &mut [u8],
    mask: &mut [u8],
) -> Option<usize> {
    let first_root = find(parent, first);
    let second_root = find(parent, second);
    if first_root == second_root {
        return None;
    }
    let (low, high) = if rank.get(first_root).copied().unwrap_or(0)
        < rank.get(second_root).copied().unwrap_or(0)
    {
        (second_root, first_root)
    } else {
        (first_root, second_root)
    };
    if let Some(slot) = parent.get_mut(high) {
        *slot = low;
    }
    if rank.get(low).copied().unwrap_or(0) == rank.get(high).copied().unwrap_or(0)
        && let Some(slot) = rank.get_mut(low)
    {
        *slot = slot.saturating_add(1);
    }
    let merged = mask.get(low).copied().unwrap_or(0) | mask.get(high).copied().unwrap_or(0);
    if let Some(slot) = mask.get_mut(low) {
        *slot = merged;
    }
    Some(low)
}

fn cross_length(first: [f64; 3], second: [f64; 3]) -> f64 {
    let x = first[1] * second[2] - first[2] * second[1];
    let y = first[2] * second[0] - first[0] * second[2];
    let z = first[0] * second[1] - first[1] * second[0];
    (x * x + y * y + z * z).sqrt()
}
