//! Integration tests for the Zeo++ bridge (Task 5.5).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_analysis::zeo::{
    ZeoOptions, parse_zeo_channels, parse_zeo_volume, run_zeo, write_zeo_inputs,
};
use xtalforge_core::{Error, Lattice, Site, Structure, StructureBuilder, SymOp};

fn inversion_structure() -> Structure {
    let inversion = SymOp::new([[-1, 0, 0], [0, -1, 0], [0, 0, -1]], [0.0; 3]).unwrap();
    StructureBuilder::new("NaCl", Lattice::cubic(4.0).unwrap())
        .site(Site::new("Na1", "Na", [0.25, 0.25, 0.25]).unwrap())
        .operation(SymOp::identity())
        .operation(inversion)
        .build()
        .unwrap()
}

const CHANNELS: &str = "\
Channels: 1
Pore Volume (A^3): 1234.5
Accessible volume (A^3): 1000.25
Surface Area (A^2): 2500.75
largest included sphere (A): 12.5
Largest Free Sphere: 9.25
Pore Limiting Diameter (A): 3.5
";

#[test]
fn channels_parser_reads_every_labelled_line() {
    let result = parse_zeo_channels(CHANNELS).unwrap();
    assert!((result.pore_volume - 1234.5).abs() < 1e-9);
    assert!((result.accessible_volume - 1000.25).abs() < 1e-9);
    assert!((result.surface_area - 2500.75).abs() < 1e-9);
    assert!((result.largest_included_sphere - 12.5).abs() < 1e-9);
    assert!((result.largest_free_sphere - 9.25).abs() < 1e-9);
    assert!((result.pore_limiting_diameter - 3.5).abs() < 1e-9);
}

#[test]
fn channels_parser_accepts_the_old_poav_spelling() {
    let text = "\
@ sample.vol Unitcell_volume: 1000 Density: 0.5 POAV_A^3: 600 ASA_A^2: 200
Largest included sphere: 10
Largest free sphere: 8
Pore limiting diameter: 6
";
    let result = parse_zeo_channels(text).unwrap();
    assert!((result.pore_volume - 600.0).abs() < 1e-12);
    assert!((result.accessible_volume - 600.0).abs() < 1e-12);
    assert!((result.surface_area - 200.0).abs() < 1e-12);
    assert!((result.largest_included_sphere - 10.0).abs() < 1e-12);
}

#[test]
fn missing_lines_are_a_typed_error() {
    let error = parse_zeo_channels("Channels: 1\nPore volume: 4.0\n").unwrap_err();
    assert!(matches!(error, Error::InvalidInput(_)));
    let Error::InvalidInput(message) = error else {
        panic!("wrong error variant");
    };
    assert!(message.contains("missing"));
    assert!(message.contains("largest free sphere"));
}

#[test]
fn volume_parser_accepts_both_spellings() {
    let old = "@ sample.vol Unitcell_volume: 29955.3 Density: 0.559382 AV_A^3: 13827.4 \
               AV_Volume_fraction: 0.4616\n";
    assert!((parse_zeo_volume(old).unwrap() - 13827.4).abs() < 1e-9);
    let occupiable = "@ sample.vol POAV_A^3: 22202.9 POAV_Volume_fraction: 0.7412\n";
    assert!((parse_zeo_volume(occupiable).unwrap() - 22202.9).abs() < 1e-9);
    assert!(parse_zeo_volume("nothing here\n").is_err());
}

#[test]
fn writers_expand_to_p1_and_document_the_arguments() {
    let structure = inversion_structure();
    let options = ZeoOptions::default();
    let (cssr, arguments) = write_zeo_inputs(&structure, &options).unwrap();
    let read_back = xtalforge_io::read_cssr(&cssr).unwrap();
    assert_eq!(read_back.sites().len(), 2);
    let atom_lines = cssr.lines().skip(4).filter(|line| !line.is_empty()).count();
    assert_eq!(atom_lines, 2);
    assert!(arguments.contains("-res result.res"));
    assert!(arguments.contains("-probe_radius 1.2"));
    assert!(arguments.contains("-samples 1000"));
    assert!(!arguments.contains("-chan"));

    let channels = ZeoOptions {
        channels: true,
        ..options
    };
    let (_, arguments) = write_zeo_inputs(&structure, &channels).unwrap();
    assert!(arguments.contains("-chan result.chan"));
}

#[test]
fn invalid_options_and_empty_structures_are_rejected() {
    let structure = inversion_structure();
    assert!(
        write_zeo_inputs(
            &structure,
            &ZeoOptions {
                probe_radius: -1.0,
                ..ZeoOptions::default()
            }
        )
        .is_err()
    );
    assert!(
        write_zeo_inputs(
            &structure,
            &ZeoOptions {
                samples: 0,
                ..ZeoOptions::default()
            }
        )
        .is_err()
    );
    let empty = StructureBuilder::new("empty", Lattice::cubic(4.0).unwrap())
        .build()
        .unwrap();
    assert!(write_zeo_inputs(&empty, &ZeoOptions::default()).is_err());
}

#[test]
fn probe_without_the_binary_returns_a_typed_error() {
    let error = run_zeo(&inversion_structure(), &ZeoOptions::default()).unwrap_err();
    assert!(matches!(error, Error::InvalidInput(_)));
    let Error::InvalidInput(message) = error else {
        panic!("wrong error variant");
    };
    assert!(message.contains("zeo++ is not available"));
}
