//! UFF/UFF4MOF over the tmQM metals (Task 0.9).
//!
//! The tmQM subset spans Sc through Ir, one complex per metal. A successful
//! build means the typer found a UFF type for the metal, its donor atoms and
//! the ligand environment; the energy and forces must be finite. The energies
//! are not compared with a published value: tmQM reports GFN2-xTB and DFT
//! numbers, not UFF.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_ff::{CalculatorOptions, UffModel};
use xtalforge_io::read_xyz;

macro_rules! xyz {
    ($name:literal) => {
        (
            $name,
            include_str!(concat!("../../../tests/external/tmqm/", $name)),
        )
    };
}

const FILES: &[(&str, &str)] = &[
    xyz!("ANAXOE_tmQM.xyz"),
    xyz!("CEYJEZ_tmQM.xyz"),
    xyz!("CILXOO_tmQM.xyz"),
    xyz!("CPMINC_tmQM.xyz"),
    xyz!("DAQHAI_tmQM.xyz"),
    xyz!("DASFIQ_tmQM.xyz"),
    xyz!("HEWTEM_tmQM.xyz"),
    xyz!("IPOSAK_tmQM.xyz"),
    xyz!("LAXLAB_tmQM.xyz"),
    xyz!("PITRAM_tmQM.xyz"),
    xyz!("QOSJIR_tmQM.xyz"),
    xyz!("REJFUL_tmQM.xyz"),
    xyz!("UMAHUO_tmQM.xyz"),
    xyz!("UTANIP_tmQM.xyz"),
    xyz!("VUDZEC_tmQM.xyz"),
    xyz!("WUGVEZ_tmQM.xyz"),
    xyz!("XARWUJ_tmQM.xyz"),
    xyz!("ZOYPEL_tmQM.xyz"),
];

#[test]
fn uff_types_and_evaluates_every_metal() {
    for (name, text) in FILES {
        let structure = read_xyz(text).unwrap_or_else(|error| panic!("{name}: {error}"));
        let metal = structure.sites().first().expect("first atom").element();
        let model = UffModel::build(&structure, CalculatorOptions::default())
            .unwrap_or_else(|error| panic!("{name}: {error}"));
        assert_eq!(model.types().len(), structure.sites().len(), "{name}");
        assert!(
            model
                .types()
                .first()
                .is_some_and(|kind| kind.starts_with(metal)),
            "{name}: first type {:?} is not a {metal} type",
            model.types().first()
        );
        let (energy, forces) = model
            .energy_forces()
            .unwrap_or_else(|error| panic!("{name}: {error}"));
        assert!(energy.is_finite(), "{name} energy {energy}");
        assert_eq!(forces.len(), structure.sites().len(), "{name} forces");
        assert!(
            forces
                .iter()
                .all(|force| force.iter().all(|value| value.is_finite())),
            "{name} has a non-finite force"
        );
        let largest = forces
            .iter()
            .map(|force| (force[0] * force[0] + force[1] * force[1] + force[2] * force[2]).sqrt())
            .fold(0.0_f64, f64::max);
        assert!(largest.is_finite() && largest >= 0.0, "{name} max force");
    }
}
