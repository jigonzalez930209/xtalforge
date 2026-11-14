//! Integration tests for magnetic neutron powder diffraction (Task 19.3.2).
//!
//! # Provenance of the hand calculation
//!
//! Body-centred iron is written by hand with `a = 2.8665 Angstrom`. The
//! antiferromagnetic ordering has `+2.2 mu_B` at the corner and `-2.2 mu_B` at
//! the body centre along `c`. The expected magnetic structure factor follows
//! the module formula with `2.695e-5 Angstrom / mu_B` and the module's
//! normalized free-atom form factor, which the test evaluates through the same
//! public `magnetic_form_factor` function. No measured magnetic form factor or
//! magnetic structure from a database is used.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_analysis::{
    InstrumentProfile, MagneticOptions, MagneticPattern, NeutronOptions, NeutronPattern,
    RadiationSource, magnetic_form_factor, magnetic_moment_cartesian, magnetic_powder,
    neutron_powder, site_moment,
};
use xtalforge_core::{Lattice, MagneticAxis, MagneticMoment, Site, Structure, StructureBuilder};

const FE_A: f64 = 2.8665;
const MOMENT: f64 = 2.2;
const WAVELENGTH: f64 = 2.5;

fn iron_site(label: &str, frac: [f64; 3], moment: f64) -> Site {
    Site::new(label, "Fe", frac)
        .unwrap()
        .with_collinear_moment(MagneticAxis::C, moment)
        .unwrap()
}

fn bcc_iron(moment: f64) -> Structure {
    StructureBuilder::new("bcc Fe", Lattice::cubic(FE_A).unwrap())
        .site(iron_site("Fe1", [0.0, 0.0, 0.0], moment))
        .site(iron_site("Fe2", [0.5, 0.5, 0.5], -moment))
        .build()
        .unwrap()
}

fn bcc_iron_without_moments() -> Structure {
    StructureBuilder::new("bcc Fe", Lattice::cubic(FE_A).unwrap())
        .site(Site::new("Fe1", "Fe", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("Fe2", "Fe", [0.5, 0.5, 0.5]).unwrap())
        .build()
        .unwrap()
}

fn profile() -> InstrumentProfile {
    InstrumentProfile {
        radiation: RadiationSource::neutron(WAVELENGTH).unwrap(),
        shape: xtalforge_analysis::LineShape::Gaussian,
        fwhm: 0.2,
        ..InstrumentProfile::default()
    }
}

fn magnetic_options() -> MagneticOptions {
    MagneticOptions {
        profile: profile(),
        two_theta_range: [40.0, 90.0],
        step: 0.02,
        ..MagneticOptions::default()
    }
}

fn nuclear_pattern(structure: &Structure) -> NeutronPattern {
    neutron_powder(
        structure,
        &NeutronOptions {
            geometry: xtalforge_analysis::NeutronGeometry::ConstantWavelength {
                wavelength: WAVELENGTH,
            },
            profile: profile(),
            range: [40.0, 90.0],
            step: 0.02,
        },
    )
    .unwrap()
}

fn lorentz(two_theta: f64) -> f64 {
    let theta = 0.5 * two_theta.to_radians();
    1.0 / (theta.sin().powi(2) * theta.cos())
}

fn two_theta(d: f64) -> f64 {
    2.0 * (0.5 * WAVELENGTH / d).asin().to_degrees()
}

fn nearest(pattern: &MagneticPattern, x: f64) -> f64 {
    pattern
        .points
        .iter()
        .min_by(|left, right| (left.x - x).abs().total_cmp(&(right.x - x).abs()))
        .map_or(0.0, |point| point.intensity)
}

#[test]
fn antiferromagnetic_bcc_iron_has_a_magnetic_100_and_no_magnetic_110() {
    let structure = bcc_iron(MOMENT);
    let magnetic = magnetic_powder(&structure, &magnetic_options()).unwrap();
    let nuclear = nuclear_pattern(&structure);

    let d100 = FE_A;
    let s = 0.5 / d100;
    let form = magnetic_form_factor("Fe", s).unwrap();
    let per_atom = 2.695e-5 * MOMENT * form;
    let expected_squared = (2.0 * per_atom).powi(2) * 2.0 / 3.0;
    let reflection = magnetic
        .reflections
        .iter()
        .find(|reflection| reflection.hkl == [1, 0, 0])
        .expect("magnetic (100)");
    assert!((reflection.d - d100).abs() < 1e-12);
    assert_eq!(reflection.label(), "(1 0 0)");
    assert!((reflection.two_theta - two_theta(d100)).abs() < 1e-12);
    assert!(
        (reflection.magnetic_amplitude_squared - expected_squared).abs() < 1e-12 * expected_squared
    );
    assert_eq!(reflection.multiplicity, 1);
    assert!(
        (reflection.intensity - expected_squared * lorentz(reflection.two_theta)).abs()
            < 1e-12 * reflection.intensity
    );

    // The body-centring phase makes the magnetic (110) factor vanish while
    // the nuclear (100) reflection is systematically absent.
    assert!(
        !magnetic
            .reflections
            .iter()
            .any(|reflection| reflection.hkl == [1, 1, 0])
    );
    assert!(
        !nuclear
            .reflections
            .iter()
            .any(|reflection| reflection.hkl == [1, 0, 0])
    );
    assert!(
        nuclear
            .reflections
            .iter()
            .any(|reflection| reflection.hkl == [1, 1, 0])
    );

    assert_eq!(magnetic.magnetic_sites.len(), 2);
    assert_eq!(magnetic.magnetic_sites[0].label, "Fe1");
    assert_eq!(magnetic.magnetic_sites[0].moment, [0.0, 0.0, MOMENT]);
    assert_eq!(magnetic.magnetic_sites[1].moment, [0.0, 0.0, -MOMENT]);
    assert!(magnetic.missing_moments.is_empty());

    let magnetic_100 = nearest(&magnetic, two_theta(d100));
    let magnetic_110 = nearest(&magnetic, two_theta(FE_A / 2.0_f64.sqrt()));
    assert!(magnetic_100 > 0.0);
    assert!(magnetic_110 < 1e-6 * magnetic.max_intensity());

    let total = magnetic.with_nuclear(&nuclear).unwrap();
    assert_eq!(total.len(), nuclear.points.len());
    let total_100 = total
        .iter()
        .min_by(|left, right| {
            (left.x - two_theta(d100))
                .abs()
                .total_cmp(&(right.x - two_theta(d100)).abs())
        })
        .unwrap()
        .intensity;
    assert!(total_100 > 0.0);
}

#[test]
fn fixed_domain_uses_the_perpendicular_moment_component() {
    let structure = StructureBuilder::new("ferro Fe", Lattice::cubic(FE_A).unwrap())
        .site(iron_site("Fe1", [0.0, 0.0, 0.0], MOMENT))
        .build()
        .unwrap();
    let options = MagneticOptions {
        domain_average: false,
        ..magnetic_options()
    };
    let magnetic = magnetic_powder(&structure, &options).unwrap();
    // The moment is along c: (001) has g parallel to the moment and vanishes,
    // while (100) has the full perpendicular component.
    assert!(
        !magnetic
            .reflections
            .iter()
            .any(|reflection| reflection.hkl == [0, 0, 1])
    );
    let reflection = magnetic
        .reflections
        .iter()
        .find(|reflection| reflection.hkl == [1, 0, 0])
        .expect("fixed-domain magnetic (100)");
    let form = magnetic_form_factor("Fe", 0.5 / FE_A).unwrap();
    let expected = (2.695e-5 * MOMENT * form).powi(2);
    assert!((reflection.magnetic_amplitude_squared - expected).abs() < 1e-12 * expected);
}

#[test]
fn zero_moments_give_an_empty_magnetic_pattern_and_equal_the_nuclear_one() {
    let structure = bcc_iron(0.0);
    let magnetic = magnetic_powder(&structure, &magnetic_options()).unwrap();
    let nuclear = nuclear_pattern(&structure);
    assert!(magnetic.reflections.is_empty());
    assert_eq!(magnetic.max_intensity(), 0.0);
    assert_eq!(magnetic.magnetic_sites.len(), 2);
    assert!(
        magnetic
            .magnetic_sites
            .iter()
            .all(|site| site.moment == [0.0; 3])
    );
    let total = magnetic.with_nuclear(&nuclear).unwrap();
    for (combined, original) in total.iter().zip(&nuclear.points) {
        assert_eq!(combined.intensity, original.intensity);
        assert_eq!(combined.x, original.x);
    }
}

#[test]
fn missing_moments_are_refused_or_reported() {
    let structure = bcc_iron_without_moments();
    let error = magnetic_powder(&structure, &magnetic_options()).unwrap_err();
    assert_eq!(error.code(), "invalid_input");
    assert!(error.to_string().contains("Fe1"));

    let options = MagneticOptions {
        allow_missing_moments: true,
        ..magnetic_options()
    };
    let pattern = magnetic_powder(&structure, &options).unwrap();
    assert_eq!(pattern.missing_moments, vec!["Fe1", "Fe2"]);
    assert!(pattern.reflections.is_empty());
    assert!(pattern.magnetic_sites.is_empty());
    assert_eq!(pattern.max_intensity(), 0.0);
}

#[test]
fn form_factor_and_moment_conversion_follow_the_documented_model() {
    let one = magnetic_form_factor("Fe", 0.0).unwrap();
    assert!((one - 1.0).abs() < 1e-12);
    let low = magnetic_form_factor("Fe", 0.1).unwrap();
    let high = magnetic_form_factor("Fe", 0.5).unwrap();
    assert!(low > high);
    assert!(high > 0.0);
    assert!(magnetic_form_factor("Fe", -0.1).is_err());
    assert!(magnetic_form_factor("X", 0.1).is_err());

    let structure = bcc_iron(MOMENT);
    let moment = site_moment(&structure.sites()[0]).unwrap();
    assert_eq!(moment.components(), [0.0, 0.0, MOMENT]);
    assert!(site_moment(&bcc_iron_without_moments().sites()[0]).is_none());

    let cubic = Lattice::cubic(4.0).unwrap();
    let collinear = MagneticMoment::collinear(MagneticAxis::C, MOMENT).unwrap();
    let cartesian = magnetic_moment_cartesian(&cubic, &collinear);
    for (actual, expected) in cartesian.iter().zip([0.0, 0.0, MOMENT]) {
        assert!((actual - expected).abs() < 1e-12);
    }
    let vector = MagneticMoment::vector([1.0, -2.0, 0.5]).unwrap();
    let cartesian = magnetic_moment_cartesian(&cubic, &vector);
    for (actual, expected) in cartesian.iter().zip([1.0, -2.0, 0.5]) {
        assert!((actual - expected).abs() < 1e-12);
    }
}

#[test]
fn invalid_magnetic_requests_return_typed_errors() {
    let structure = bcc_iron(MOMENT);
    let mut bad = magnetic_options();
    bad.profile = InstrumentProfile::default();
    assert!(magnetic_powder(&structure, &bad).is_err());
    bad = magnetic_options();
    bad.two_theta_range = [90.0, 40.0];
    assert!(magnetic_powder(&structure, &bad).is_err());
    bad = magnetic_options();
    bad.step = 0.0;
    assert!(magnetic_powder(&structure, &bad).is_err());
    bad = magnetic_options();
    bad.profile.radiation = RadiationSource {
        kind: xtalforge_analysis::RadiationKind::Neutron,
        wavelength: 0.0,
    };
    assert!(magnetic_powder(&structure, &bad).is_err());
    bad = magnetic_options();
    bad.two_theta_range = [5.0, 179.0];
    bad.step = 1e-6;
    assert!(magnetic_powder(&structure, &bad).is_err());

    let wrong_grid = neutron_powder(
        &structure,
        &NeutronOptions {
            geometry: xtalforge_analysis::NeutronGeometry::ConstantWavelength {
                wavelength: WAVELENGTH,
            },
            profile: profile(),
            range: [40.0, 90.0],
            step: 0.03,
        },
    )
    .unwrap();
    let magnetic = magnetic_powder(&structure, &magnetic_options()).unwrap();
    assert!(magnetic.with_nuclear(&wrong_grid).is_err());
}
