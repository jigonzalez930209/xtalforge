//! Sequential and diverging colour maps (Task 3.14.3).

use serde::{Deserialize, Serialize};

/// Built-in colour map.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub enum Colormap {
    /// Perceptually uniform blue-green-yellow.
    #[default]
    Viridis,
    /// Blue-purple-orange.
    Plasma,
    /// Black-red-yellow.
    Inferno,
    /// Black-purple-orange-white.
    Magma,
    /// Blue-white-red diverging map.
    CoolWarm,
    /// Black to white.
    Grayscale,
}

impl Colormap {
    /// Samples the map at `t` in `[0, 1]` and returns 0-255 RGB.
    #[must_use]
    pub fn sample(self, t: f64) -> [u8; 3] {
        let t = t.clamp(0.0, 1.0);
        match self {
            Self::Viridis => interpolate(
                &[
                    [68, 1, 84],
                    [59, 82, 139],
                    [33, 145, 140],
                    [94, 201, 98],
                    [253, 231, 37],
                ],
                t,
            ),
            Self::Plasma => interpolate(
                &[
                    [13, 8, 135],
                    [126, 3, 168],
                    [204, 71, 120],
                    [248, 149, 64],
                    [240, 249, 33],
                ],
                t,
            ),
            Self::Inferno => interpolate(
                &[
                    [0, 0, 4],
                    [87, 16, 110],
                    [188, 55, 84],
                    [249, 142, 9],
                    [252, 255, 164],
                ],
                t,
            ),
            Self::Magma => interpolate(
                &[
                    [0, 0, 4],
                    [81, 18, 124],
                    [183, 55, 121],
                    [252, 137, 97],
                    [252, 253, 191],
                ],
                t,
            ),
            Self::CoolWarm => interpolate(&[[59, 76, 192], [221, 221, 221], [180, 4, 38]], t),
            Self::Grayscale => {
                let value = (t * 255.0).round() as u8;
                [value, value, value]
            }
        }
    }
}

fn interpolate(colors: &[[u8; 3]], t: f64) -> [u8; 3] {
    if colors.len() < 2 {
        return colors.first().copied().unwrap_or([0, 0, 0]);
    }
    let scaled = t * (colors.len() - 1) as f64;
    let index = (scaled.floor() as usize).min(colors.len().saturating_sub(2));
    let fraction = scaled - index as f64;
    let (Some(first), Some(second)) = (colors.get(index), colors.get(index + 1)) else {
        return [0, 0, 0];
    };
    std::array::from_fn(|axis| {
        let first = f64::from(first.get(axis).copied().unwrap_or(0));
        let second = f64::from(second.get(axis).copied().unwrap_or(0));
        (first + fraction * (second - first)).round() as u8
    })
}
