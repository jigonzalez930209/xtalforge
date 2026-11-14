//! Rendered SAED image buffer (Task 19.3.3).

use serde::{Deserialize, Serialize};

/// Smallest accepted image edge in pixels.
pub(crate) const MIN_IMAGE_SIZE: usize = 16;
/// Largest accepted image edge in pixels.
pub(crate) const MAX_IMAGE_SIZE: usize = 2048;

/// A square rendered diffraction image, normalized to a maximum of one.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SaedImage {
    /// Image width in pixels.
    pub width: usize,
    /// Image height in pixels.
    pub height: usize,
    /// Row-major pixels in `[0, 1]`.
    pub pixels: Vec<f32>,
}

impl SaedImage {
    /// Largest pixel value.
    #[must_use]
    pub fn max_intensity(&self) -> f32 {
        self.pixels.iter().copied().fold(0.0_f32, f32::max)
    }

    /// Pixel value at `(x, y)`, or `None` outside the image.
    #[must_use]
    pub fn pixel(&self, x: usize, y: usize) -> Option<f32> {
        if x >= self.width || y >= self.height {
            return None;
        }
        let index = y.checked_mul(self.width)?.checked_add(x)?;
        self.pixels.get(index).copied()
    }

    /// Converts the normalized image to 8-bit grey levels.
    #[must_use]
    pub fn to_u8(&self) -> Vec<u8> {
        let maximum = self.max_intensity();
        if maximum <= 0.0 {
            return vec![0; self.pixels.len()];
        }
        self.pixels
            .iter()
            .map(|pixel| (pixel / maximum * 255.0).round().clamp(0.0, 255.0) as u8)
            .collect()
    }
}

/// Adds one Gaussian spot to a square buffer.
///
/// The Gaussian is truncated at four standard deviations. The truncation is
/// symmetric in the pixel centres so a mirrored spot renders a mirrored disc.
pub(crate) fn add_spot(
    pixels: &mut [f32],
    size: usize,
    x: f64,
    y: f64,
    sigma: f64,
    amplitude: f64,
) {
    let reach = 4.0 * sigma;
    let start_x = (x - reach - 0.5).floor() as i64;
    let end_x = (x + reach - 0.5).ceil() as i64;
    let start_y = (y - reach - 0.5).floor() as i64;
    let end_y = (y + reach - 0.5).ceil() as i64;
    let inverse = 1.0 / (2.0 * sigma * sigma);
    for row in start_y..=end_y {
        if row < 0 || row >= size as i64 {
            continue;
        }
        for column in start_x..=end_x {
            if column < 0 || column >= size as i64 {
                continue;
            }
            let dx = column as f64 + 0.5 - x;
            let dy = row as f64 + 0.5 - y;
            if dx * dx + dy * dy > reach * reach {
                continue;
            }
            let value = amplitude * (-(dx * dx + dy * dy) * inverse).exp();
            if let Some(pixel) = pixels.get_mut(row as usize * size + column as usize) {
                *pixel += value as f32;
            }
        }
    }
}
