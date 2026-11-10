//! Pictures for the agent (Task 8.1.6).
//!
//! Crystal-Builder renders with VTK. This module keeps the intent — an agent
//! can look at a structure — with a small pure-Rust orthographic renderer of
//! the scene data (Task 3.14.4). It writes SVG for documents and PNG for
//! clients that only accept raster images. There is no shading beyond one
//! fixed light and no transparency; the picture is good enough to see a
//! folded linker or a misplaced cluster.

use serde::{Deserialize, Serialize};
use xtalforge_core::Structure;
use xtalforge_view::{ViewSettings, build_scene};

use crate::Error;
use crate::png::RgbImage;

/// Smallest accepted image dimension.
pub const MIN_DIMENSION: u32 = 32;
/// Largest accepted image dimension.
pub const MAX_DIMENSION: u32 = 4096;

/// How to draw a picture.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RenderOptions {
    /// Width in pixels.
    pub width: u32,
    /// Height in pixels.
    pub height: u32,
    /// Camera direction: the vector from the structure to the camera.
    pub direction: [f64; 3],
    /// Up hint in world coordinates.
    pub up: [f64; 3],
    /// Margin around the structure in pixels.
    pub margin: f64,
    /// View settings for radii, colours and the cell.
    pub settings: ViewSettings,
}

impl Default for RenderOptions {
    fn default() -> Self {
        Self {
            width: 800,
            height: 600,
            direction: [1.0, -1.4, 0.9],
            up: [0.0, 0.0, 1.0],
            margin: 24.0,
            settings: ViewSettings::default(),
        }
    }
}

impl RenderOptions {
    /// Checks the dimensions and normalises the camera.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Argument`] for a dimension outside the accepted range
    /// or a degenerate camera.
    pub fn validate(&self) -> Result<(), Error> {
        if !(MIN_DIMENSION..=MAX_DIMENSION).contains(&self.width)
            || !(MIN_DIMENSION..=MAX_DIMENSION).contains(&self.height)
        {
            return Err(Error::Argument {
                name: "width/height".to_owned(),
                reason: format!("must be between {MIN_DIMENSION} and {MAX_DIMENSION} pixels"),
            });
        }
        if !self.direction.iter().all(|value| value.is_finite()) || norm(self.direction) < 1e-9 {
            return Err(Error::Argument {
                name: "direction".to_owned(),
                reason: "must be a finite non-zero vector".to_owned(),
            });
        }
        Ok(())
    }
}

fn norm(vector: [f64; 3]) -> f64 {
    (vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]).sqrt()
}

fn normalised(vector: [f64; 3]) -> [f64; 3] {
    let length = norm(vector);
    if length == 0.0 {
        [0.0, 0.0, 0.0]
    } else {
        [vector[0] / length, vector[1] / length, vector[2] / length]
    }
}

fn dot(left: [f64; 3], right: [f64; 3]) -> f64 {
    left[0] * right[0] + left[1] * right[1] + left[2] * right[2]
}

fn cross(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
    [
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    ]
}

fn subtract(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
    [left[0] - right[0], left[1] - right[1], left[2] - right[2]]
}

struct Camera {
    center: [f64; 3],
    right: [f64; 3],
    up: [f64; 3],
    forward: [f64; 3],
    scale: f64,
    width: f64,
    height: f64,
}

struct Projected {
    x: f64,
    y: f64,
    depth: f64,
    radius: f64,
}

impl Camera {
    fn project(&self, point: [f64; 3], radius_world: f64) -> Projected {
        let relative = subtract(point, self.center);
        Projected {
            x: self.width / 2.0 + dot(relative, self.right) * self.scale,
            y: self.height / 2.0 - dot(relative, self.up) * self.scale,
            depth: dot(relative, self.forward),
            radius: radius_world * self.scale,
        }
    }
}

fn camera(structure: &Structure, options: &RenderOptions) -> Result<Camera, Error> {
    let scene = build_scene(structure, &options.settings)?;
    let forward = normalised(options.direction);
    let mut right = cross(options.up, forward);
    if norm(right) < 1e-9 {
        right = cross([1.0, 0.0, 0.0], forward);
    }
    let right = normalised(right);
    let up = normalised(cross(forward, right));

    let mut minimum = [f64::INFINITY; 3];
    let mut maximum = [f64::NEG_INFINITY; 3];
    let mut include = |point: [f64; 3], radius: f64| {
        for ((low, high), value) in minimum.iter_mut().zip(maximum.iter_mut()).zip(point.iter()) {
            *low = low.min(*value - radius);
            *high = high.max(*value + radius);
        }
    };
    for atom in &scene.atoms {
        if !atom.hidden {
            include(atom.position, atom.radius.max(0.2));
        }
    }
    for (start, end) in &scene.cell_edges {
        include(*start, 0.0);
        include(*end, 0.0);
    }
    if !minimum[0].is_finite() {
        return Err(Error::Core(xtalforge_core::Error::InvalidInput(
            "a picture needs at least one visible atom".to_owned(),
        )));
    }
    let center = [
        (minimum[0] + maximum[0]) / 2.0,
        (minimum[1] + maximum[1]) / 2.0,
        (minimum[2] + maximum[2]) / 2.0,
    ];
    let mut half_width: f64 = 0.0;
    let mut half_height: f64 = 0.0;
    let mut measure = |point: [f64; 3], radius: f64| {
        let relative = subtract(point, center);
        half_width = half_width.max(dot(relative, right).abs() + radius);
        half_height = half_height.max(dot(relative, up).abs() + radius);
    };
    for atom in &scene.atoms {
        if !atom.hidden {
            measure(atom.position, atom.radius.max(0.2));
        }
    }
    for (start, end) in &scene.cell_edges {
        measure(*start, 0.0);
        measure(*end, 0.0);
    }
    let usable_width = (f64::from(options.width) - 2.0 * options.margin).max(1.0);
    let usable_height = (f64::from(options.height) - 2.0 * options.margin).max(1.0);
    let scale = (usable_width / (2.0 * half_width.max(1e-6)))
        .min(usable_height / (2.0 * half_height.max(1e-6)));
    Ok(Camera {
        center,
        right,
        up,
        forward,
        scale,
        width: f64::from(options.width),
        height: f64::from(options.height),
    })
}

fn shade(color: [u8; 3], normal: [f64; 3]) -> [u8; 3] {
    let light = normalised([-0.4, 0.5, 0.75]);
    let diffuse = dot(normal, light).max(0.0);
    let intensity = 0.35 + 0.65 * diffuse;
    [
        (f64::from(color[0]) * intensity).clamp(0.0, 255.0) as u8,
        (f64::from(color[1]) * intensity).clamp(0.0, 255.0) as u8,
        (f64::from(color[2]) * intensity).clamp(0.0, 255.0) as u8,
    ]
}

/// Renders a structure to SVG.
///
/// # Errors
///
/// Returns the scene, P1 expansion or argument error.
pub fn render_svg(structure: &Structure, options: &RenderOptions) -> Result<String, Error> {
    options.validate()?;
    let camera = camera(structure, options)?;
    let scene = build_scene(structure, &options.settings)?;
    let mut svg = String::with_capacity(16_384);
    svg.push_str(&format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{}\" height=\"{}\" viewBox=\"0 0 {} {}\">\n",
        options.width, options.height, options.width, options.height
    ));
    svg.push_str(&format!(
        "<rect width=\"{}\" height=\"{}\" fill=\"#{:02x}{:02x}{:02x}\"/>\n",
        options.width,
        options.height,
        options.settings.background[0],
        options.settings.background[1],
        options.settings.background[2]
    ));

    // Cell edges first: they are context, not geometry.
    let mut cell = String::new();
    for (start, end) in &scene.cell_edges {
        let first = camera.project(*start, 0.0);
        let second = camera.project(*end, 0.0);
        cell.push_str(&format!(
            "<line x1=\"{:.1}\" y1=\"{:.1}\" x2=\"{:.1}\" y2=\"{:.1}\" stroke=\"#9aa0a6\" stroke-width=\"1\"/>\n",
            first.x, first.y, second.x, second.y
        ));
    }

    // Bonds and atoms share one depth order so a bond behind a sphere stays
    // behind it.
    let mut primitives: Vec<(f64, String)> = Vec::new();
    for bond in &scene.bonds {
        let first = camera.project(bond.start, 0.0);
        let second = camera.project(bond.end, 0.0);
        let radius = (bond.radius * camera.scale).max(1.0);
        primitives.push((
            (first.depth + second.depth) / 2.0 - radius,
            format!(
                "<line x1=\"{:.1}\" y1=\"{:.1}\" x2=\"{:.1}\" y2=\"{:.1}\" stroke=\"#{:02x}{:02x}{:02x}\" stroke-width=\"{:.1}\" stroke-linecap=\"round\"/>\n",
                first.x, first.y, second.x, second.y, bond.color[0], bond.color[1], bond.color[2], radius * 2.0
            ),
        ));
    }
    for atom in &scene.atoms {
        if atom.hidden || atom.radius <= 0.0 {
            continue;
        }
        let projected = camera.project(atom.position, atom.radius);
        primitives.push((
            projected.depth - projected.radius,
            format!(
                "<circle cx=\"{:.1}\" cy=\"{:.1}\" r=\"{:.1}\" fill=\"#{:02x}{:02x}{:02x}\"/>\n",
                projected.x,
                projected.y,
                projected.radius,
                atom.color[0],
                atom.color[1],
                atom.color[2]
            ),
        ));
    }
    primitives.sort_by(|left, right| left.0.total_cmp(&right.0));
    svg.push_str(&cell);
    for (_, primitive) in primitives {
        svg.push_str(&primitive);
    }
    svg.push_str("</svg>\n");
    Ok(svg)
}

/// Renders a structure to PNG.
///
/// # Errors
///
/// Returns the scene, P1 expansion, argument or encoder error.
pub fn render_png(structure: &Structure, options: &RenderOptions) -> Result<Vec<u8>, Error> {
    options.validate()?;
    let camera = camera(structure, options)?;
    let scene = build_scene(structure, &options.settings)?;
    let mut image = RgbImage::filled(options.width, options.height, options.settings.background);
    let mut depth_buffer =
        vec![f64::NEG_INFINITY; options.width as usize * options.height as usize];

    for (start, end) in &scene.cell_edges {
        let first = camera.project(*start, 0.0);
        let second = camera.project(*end, 0.0);
        draw_capsule(
            &mut image,
            &mut depth_buffer,
            first,
            second,
            0.75,
            [154, 160, 166],
        );
    }
    for bond in &scene.bonds {
        let first = camera.project(bond.start, 0.0);
        let second = camera.project(bond.end, 0.0);
        draw_capsule(
            &mut image,
            &mut depth_buffer,
            first,
            second,
            (bond.radius * camera.scale).max(0.75),
            bond.color,
        );
    }
    for atom in &scene.atoms {
        if atom.hidden || atom.radius <= 0.0 {
            continue;
        }
        let projected = camera.project(atom.position, atom.radius);
        draw_sphere(&mut image, &mut depth_buffer, projected, atom.color);
    }
    image
        .to_png()
        .map_err(|error| Error::Io(format!("PNG encoder: {error}")))
}

fn draw_sphere(image: &mut RgbImage, depth: &mut [f64], projected: Projected, color: [u8; 3]) {
    let radius = projected.radius.max(0.75);
    let left = (projected.x - radius).floor().max(0.0) as u32;
    let right = (projected.x + radius).ceil().min(f64::from(image.width())) as u32;
    let top = (projected.y - radius).floor().max(0.0) as u32;
    let bottom = (projected.y + radius).ceil().min(f64::from(image.height())) as u32;
    for y in top..bottom {
        for x in left..right {
            let dx = f64::from(x) + 0.5 - projected.x;
            let dy = f64::from(y) + 0.5 - projected.y;
            let distance_squared = dx * dx + dy * dy;
            if distance_squared > radius * radius {
                continue;
            }
            let offset = (radius * radius - distance_squared).sqrt();
            let z = projected.depth + offset;
            let index = y as usize * image.width() as usize + x as usize;
            if let Some(slot) = depth.get_mut(index)
                && z > *slot
            {
                *slot = z;
                let normal = normalised([dx / radius, -dy / radius, offset / radius]);
                image.set(x, y, shade(color, normal));
            }
        }
    }
}

fn draw_capsule(
    image: &mut RgbImage,
    depth: &mut [f64],
    start: Projected,
    end: Projected,
    radius: f64,
    color: [u8; 3],
) {
    let left = (start.x.min(end.x) - radius).floor().max(0.0) as u32;
    let right = (start.x.max(end.x) + radius)
        .ceil()
        .min(f64::from(image.width())) as u32;
    let top = (start.y.min(end.y) - radius).floor().max(0.0) as u32;
    let bottom = (start.y.max(end.y) + radius)
        .ceil()
        .min(f64::from(image.height())) as u32;
    let axis = [end.x - start.x, end.y - start.y];
    let length_squared = axis[0] * axis[0] + axis[1] * axis[1];
    for y in top..bottom {
        for x in left..right {
            let point = [f64::from(x) + 0.5 - start.x, f64::from(y) + 0.5 - start.y];
            let t = if length_squared > 0.0 {
                ((point[0] * axis[0] + point[1] * axis[1]) / length_squared).clamp(0.0, 1.0)
            } else {
                0.0
            };
            let nearest = [axis[0] * t, axis[1] * t];
            let dx = point[0] - nearest[0];
            let dy = point[1] - nearest[1];
            let distance_squared = dx * dx + dy * dy;
            if distance_squared > radius * radius {
                continue;
            }
            let z = start.depth + (end.depth - start.depth) * t;
            let index = y as usize * image.width() as usize + x as usize;
            if let Some(slot) = depth.get_mut(index)
                && z > *slot
            {
                *slot = z;
                let normal = normalised([dx / radius, -dy / radius, 0.35]);
                image.set(x, y, shade(color, normal));
            }
        }
    }
}
