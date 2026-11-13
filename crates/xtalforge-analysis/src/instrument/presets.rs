//! Named instrument presets (Task 14.1.5).
//!
//! A preset is a starting profile for a common diffractometer. The numbers are
//! typical literature values, not a calibration: a real instrument replaces
//! them with its own Caglioti and asymmetry coefficients after a standard
//! refinement. The preset table is deterministic and available on every
//! surface.

use xtalforge_core::Error;

use super::Geometry;
use super::profile::{
    Caglioti, FcjAsymmetry, InstrumentProfile, LineShape, Monochromator, RadiationKind,
    RadiationSource, Tch,
};

/// One named instrument preset.
///
/// Presets are a static lookup table rather than a serialized data struct; the
/// embedded [`InstrumentProfile`] is serializable on its own.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InstrumentPreset {
    /// Stable lowercase identifier, for example `lab-cu-ka`.
    pub name: &'static str,
    /// Human-readable label.
    pub label: &'static str,
    /// Diffractometer geometry.
    pub geometry: Geometry,
    /// Profile seed for the preset.
    pub profile: InstrumentProfile,
}

impl InstrumentPreset {
    /// Validates the profile of the preset.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the preset profile is invalid,
    /// which can only happen after the table was edited incorrectly.
    pub fn validate(&self) -> Result<(), Error> {
        self.profile.validate()
    }
}

/// Number of entries returned by [`instrument_presets`].
pub const INSTRUMENT_PRESET_COUNT: usize = 7;

/// The preset table, in declaration order.
#[must_use]
pub fn instrument_presets() -> [InstrumentPreset; INSTRUMENT_PRESET_COUNT] {
    [
        laboratory(
            "lab-cu-ka",
            "Laboratory Cu K alpha",
            1.541_84,
            Geometry::BraggBrentano,
            Caglioti {
                u: 0.004,
                v: -0.002,
                w: 0.012,
            },
        ),
        laboratory(
            "lab-co-ka",
            "Laboratory Co K alpha",
            1.788_97,
            Geometry::BraggBrentano,
            Caglioti {
                u: 0.006,
                v: -0.004,
                w: 0.018,
            },
        ),
        laboratory(
            "lab-mo-ka",
            "Laboratory Mo K alpha",
            0.710_73,
            Geometry::BraggBrentano,
            Caglioti {
                u: 0.002,
                v: -0.001,
                w: 0.006,
            },
        ),
        laboratory(
            "capillary-cu-ka",
            "Cu K alpha capillary",
            1.541_84,
            Geometry::DebyeScherrer,
            Caglioti {
                u: 0.003,
                v: -0.001,
                w: 0.008,
            },
        ),
        InstrumentPreset {
            name: "synchrotron",
            label: "Synchrotron, 0.5 Angstrom",
            geometry: Geometry::DebyeScherrer,
            profile: InstrumentProfile {
                radiation: RadiationSource {
                    kind: RadiationKind::XRay,
                    wavelength: 0.5,
                },
                shape: LineShape::PseudoVoigt,
                fwhm: 0.01,
                caglioti: Some(Caglioti {
                    u: 0.0005,
                    v: -0.0002,
                    w: 0.001,
                }),
                tch: None,
                eta: 0.7,
                asymmetry: None,
                monochromator: Monochromator { polarization: 0.95 },
            },
        },
        InstrumentPreset {
            name: "neutron-cw",
            label: "Constant-wavelength neutron, 1.909 Angstrom",
            geometry: Geometry::DebyeScherrer,
            profile: InstrumentProfile {
                radiation: RadiationSource {
                    kind: RadiationKind::Neutron,
                    wavelength: 1.909,
                },
                shape: LineShape::Gaussian,
                fwhm: 0.4,
                caglioti: Some(Caglioti {
                    u: 0.02,
                    v: -0.01,
                    w: 0.03,
                }),
                tch: None,
                eta: 0.0,
                asymmetry: None,
                monochromator: Monochromator { polarization: 0.5 },
            },
        },
        InstrumentPreset {
            name: "electron-200kv",
            label: "Electron microscope, 200 kV",
            geometry: Geometry::DebyeScherrer,
            profile: InstrumentProfile {
                radiation: RadiationSource {
                    kind: RadiationKind::Electron,
                    wavelength: 0.025_08,
                },
                shape: LineShape::Gaussian,
                fwhm: 0.02,
                caglioti: None,
                tch: None,
                eta: 0.0,
                asymmetry: None,
                monochromator: Monochromator { polarization: 0.5 },
            },
        },
    ]
}

/// Looks a preset up by name, case-insensitively.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] naming the available presets when the name
/// is unknown.
pub fn preset_by_name(name: &str) -> Result<InstrumentPreset, Error> {
    let wanted = name.trim().to_ascii_lowercase();
    for preset in instrument_presets() {
        if preset.name == wanted {
            return Ok(preset);
        }
    }
    Err(Error::InvalidInput(format!(
        "unknown instrument preset {name:?}; available: {}",
        INSTRUMENT_PRESET_NAMES.join(", ")
    )))
}

/// Every preset identifier, in table order.
pub const INSTRUMENT_PRESET_NAMES: [&str; INSTRUMENT_PRESET_COUNT] = [
    "lab-cu-ka",
    "lab-co-ka",
    "lab-mo-ka",
    "capillary-cu-ka",
    "synchrotron",
    "neutron-cw",
    "electron-200kv",
];

/// Builds one laboratory preset with a pseudo-Voigt shape and FCJ asymmetry.
fn laboratory(
    name: &'static str,
    label: &'static str,
    wavelength: f64,
    geometry: Geometry,
    caglioti: Caglioti,
) -> InstrumentPreset {
    InstrumentPreset {
        name,
        label,
        geometry,
        profile: InstrumentProfile {
            radiation: RadiationSource {
                kind: RadiationKind::XRay,
                wavelength,
            },
            shape: LineShape::PseudoVoigt,
            fwhm: 0.1,
            caglioti: Some(caglioti),
            tch: Some(Tch {
                gaussian: caglioti,
                x: 0.002,
                y: 0.004,
            }),
            eta: 0.5,
            asymmetry: Some(FcjAsymmetry {
                fcj1: 0.02,
                fcj2: 0.03,
            }),
            monochromator: Monochromator { polarization: 0.5 },
        },
    }
}

/// TCH-aware variant of [`preset_by_name`] used by callers that want the TCH
/// width model of a laboratory preset.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an unknown name.
pub fn tch_preset_by_name(name: &str) -> Result<InstrumentProfile, Error> {
    let mut preset = preset_by_name(name)?;
    if let Some(tch) = preset.profile.tch {
        preset.profile.shape = LineShape::Tch;
        preset.profile.caglioti = None;
        preset.profile.tch = Some(tch);
    }
    Ok(preset.profile)
}
