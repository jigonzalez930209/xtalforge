//! Integration tests for the streaming readers (Task 15.2.2).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::io::{BufRead, Cursor, Read};

use xtalforge_core::Lattice;
use xtalforge_io::{
    LammpsStream, StreamFrame, XyzStream, read_cif, read_cif_stream, read_lammps_dump,
    read_trajectory_xyz,
};

const FRAMES: usize = 40;
const ATOMS: usize = 64;

fn lattice() -> Lattice {
    Lattice::new([[12.0, 0.0, 0.0], [1.5, 11.0, 0.0], [0.5, 1.0, 10.0]]).unwrap()
}

fn frame_text(lattice: &Lattice, frame: usize) -> String {
    let mut text = String::new();
    for atom in 0..ATOMS {
        let frac = [
            ((atom * 7 + frame * 3) % 97) as f64 / 97.0,
            ((atom * 11 + frame * 5) % 89) as f64 / 89.0,
            ((atom * 13 + frame * 7) % 83) as f64 / 83.0,
        ];
        let cart = lattice.to_cart(frac);
        text.push_str(&format!(
            "C {:.10} {:.10} {:.10}\n",
            cart[0], cart[1], cart[2]
        ));
    }
    text
}

fn xyz_text() -> String {
    let lattice = lattice();
    let matrix = lattice
        .matrix()
        .iter()
        .flatten()
        .map(|value| format!("{value:.10}"))
        .collect::<Vec<_>>()
        .join(" ");
    let mut text = String::new();
    for frame in 0..FRAMES {
        text.push_str(&format!(
            "{ATOMS}\nLattice=\"{matrix}\" energy={:.4} step={frame}\n",
            frame as f64 * 1.5
        ));
        text.push_str(&frame_text(&lattice, frame));
    }
    text
}

fn lammps_text() -> String {
    let lattice = lattice();
    let matrix = lattice.matrix();
    let lengths = lattice.lengths();
    let mut text = String::new();
    for frame in 0..6 {
        text.push_str("ITEM: TIMESTEP\n");
        text.push_str(&format!("{frame}\n"));
        text.push_str("ITEM: NUMBER OF ATOMS\n");
        text.push_str(&format!("{ATOMS}\n"));
        text.push_str("ITEM: BOX BOUNDS xy xz yz pp pp pp\n");
        text.push_str(&format!("0 {}\n", lengths[0]));
        text.push_str(&format!("0 {}\n", lengths[1]));
        text.push_str(&format!("0 {}\n", lengths[2]));
        text.push_str("ITEM: ATOMS id element x y z\n");
        for atom in 0..ATOMS {
            let frac = [
                ((atom * 5 + frame) % 71) as f64 / 71.0,
                ((atom * 3 + frame * 2) % 67) as f64 / 67.0,
                ((atom * 2 + frame * 3) % 61) as f64 / 61.0,
            ];
            let cart = [
                frac[0] * matrix[0][0] + frac[1] * matrix[1][0] + frac[2] * matrix[2][0],
                frac[0] * matrix[0][1] + frac[1] * matrix[1][1] + frac[2] * matrix[2][1],
                frac[0] * matrix[0][2] + frac[1] * matrix[1][2] + frac[2] * matrix[2][2],
            ];
            text.push_str(&format!(
                "{} C {:.10} {:.10} {:.10}\n",
                atom + 1,
                cart[0],
                cart[1],
                cart[2]
            ));
        }
    }
    text
}

#[test]
fn xyz_stream_matches_the_whole_file_reader() {
    let text = xyz_text();
    let expected = read_trajectory_xyz(&text).unwrap();
    let mut stream = XyzStream::new(Cursor::new(text.as_bytes()));
    let mut frames = Vec::new();
    while let Some(frame) = stream.next_frame().unwrap() {
        frames.push(frame);
    }
    assert_eq!(frames.len(), FRAMES);
    assert_eq!(frames.len(), expected.frames.len());
    for (index, (frame, expected)) in frames.iter().zip(expected.frames.iter()).enumerate() {
        assert_eq!(frame.len(), ATOMS);
        assert_eq!(frame.step, Some(index as u64));
        assert_eq!(frame.energy, Some(index as f64 * 1.5));
        assert_eq!(frame.lattice, expected.lattice);
        assert_eq!(frame.frac.len(), expected.frac.len());
        for (actual, expected) in frame.frac.iter().zip(expected.frac.iter()) {
            assert_eq!(actual, expected);
        }
        assert!(frame.elements.iter().all(|element| element == "C"));
    }
    assert_eq!(frames[0].structure("frame").unwrap().sites().len(), ATOMS);
}

#[test]
fn xyz_stream_reads_one_frame_at_a_time() {
    let text = xyz_text();
    let total = text.len();
    let mut reader = CountingReader::new(Cursor::new(text.into_bytes()));
    let mut stream = XyzStream::new(&mut reader);
    let frame = stream.next_frame().unwrap().expect("one frame");
    assert_eq!(frame.len(), ATOMS);
    let consumed = reader.consumed();
    assert!(consumed < total / 4, "consumed {consumed} of {total} bytes");
}

#[test]
fn plain_xyz_frames_stay_cartesian() {
    let mut text = String::new();
    for frame in 0..3 {
        text.push_str(&format!("2\nframe {frame}\n"));
        text.push_str(&format!("C {frame}.0 0.0 0.0\nO 0.0 1.0 0.0\n"));
    }
    let mut stream = XyzStream::new(Cursor::new(text.as_bytes()));
    let mut count = 0;
    while let Some(frame) = stream.next_frame().unwrap() {
        assert_eq!(frame.lattice, None);
        assert_eq!(frame.frac[0], [count as f64, 0.0, 0.0]);
        assert!(frame.structure("no-cell").is_err());
        count += 1;
    }
    assert_eq!(count, 3);
}

#[test]
fn lammps_stream_matches_the_whole_file_reader() {
    let text = lammps_text();
    let expected = read_lammps_dump(&text).unwrap();
    let mut stream = LammpsStream::new(Cursor::new(text.as_bytes()));
    let mut frames = Vec::new();
    while let Some(frame) = stream.next_frame().unwrap() {
        frames.push(frame);
    }
    assert_eq!(frames.len(), expected.frames.len());
    for (index, frame) in frames.iter().enumerate() {
        let Some(reference) = expected.frames.get(index) else {
            continue;
        };
        assert_eq!(frame.len(), ATOMS);
        assert_eq!(frame.step, Some(index as u64));
        assert_eq!(frame.lattice, expected.lattice);
        for (actual, expected) in frame.frac.iter().zip(reference.frac.iter()) {
            assert!((actual[0] - expected[0]).abs() < 1e-9);
            assert!((actual[1] - expected[1]).abs() < 1e-9);
            assert!((actual[2] - expected[2]).abs() < 1e-9);
        }
    }
}

#[test]
fn lammps_stream_accepts_a_type_mapping() {
    let text = "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n1\nITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\nITEM: ATOMS id type xs ys zs\n1 1 0.25 0.5 0.75\n";
    let mut stream = LammpsStream::with_types(Cursor::new(text.as_bytes()), &["Zn"]);
    let frame = stream.next_frame().unwrap().expect("one frame");
    assert_eq!(frame.elements, vec!["Zn".to_owned()]);
    assert_eq!(frame.frac[0], [0.25, 0.5, 0.75]);
    assert!(stream.next_frame().unwrap().is_none());
}

#[test]
fn lammps_stream_reads_one_frame_at_a_time() {
    let text = lammps_text();
    let total = text.len();
    let mut reader = CountingReader::new(Cursor::new(text.into_bytes()));
    let mut stream = LammpsStream::new(&mut reader);
    let frame = stream.next_frame().unwrap().expect("one frame");
    assert_eq!(frame.len(), ATOMS);
    assert!(reader.consumed() < total / 4);
}

#[test]
fn cif_stream_matches_the_document_reader() {
    for text in [
        include_str!("../../../docs/public/samples/Si.cif"),
        include_str!("../../../tests/corpus/cod/MOF-5.cif"),
        include_str!("../../../tests/corpus/cod/MIL-53.cif"),
    ] {
        let expected = read_cif(text).unwrap();
        let actual = read_cif_stream(Cursor::new(text.as_bytes())).unwrap();
        assert_eq!(actual.name(), expected.name());
        assert_eq!(actual.hall(), expected.hall());
        assert_eq!(actual.lattice(), expected.lattice());
        assert_eq!(actual.operations(), expected.operations());
        assert_eq!(actual.sites().len(), expected.sites().len());
        for (actual, expected) in actual.sites().iter().zip(expected.sites().iter()) {
            assert_eq!(actual.label(), expected.label());
            assert_eq!(actual.element(), expected.element());
            assert_eq!(actual.frac(), expected.frac());
            assert_eq!(actual.occupancy(), expected.occupancy());
            assert_eq!(actual.u_iso(), expected.u_iso());
            assert_eq!(actual.u_aniso(), expected.u_aniso());
            assert_eq!(actual.charge(), expected.charge());
            assert_eq!(actual.props(), expected.props());
            assert_eq!(actual.disorder_group(), expected.disorder_group());
        }
        assert_eq!(actual.metadata(), expected.metadata());
    }
}

#[test]
fn cif_stream_handles_a_generated_atom_loop() {
    let mut text = String::from(
        "data_generated\n_cell_length_a 20\n_cell_length_b 20\n_cell_length_c 20\n_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\nloop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n_atom_site_occupancy\n",
    );
    for atom in 0..2_000 {
        text.push_str(&format!(
            "C{atom} C {:.6} {:.6} {:.6} 1.0\n",
            (atom % 100) as f64 / 100.0,
            (atom % 97) as f64 / 97.0,
            (atom % 89) as f64 / 89.0
        ));
    }
    let structure = read_cif_stream(Cursor::new(text.as_bytes())).unwrap();
    assert_eq!(structure.sites().len(), 2_000);
    assert_eq!(structure.sites()[1_999].label(), "C1999");
    assert!(structure.is_p1());
}

#[test]
fn cif_stream_rejects_explicit_bonds_and_early_cartesian_loops() {
    let bonds = "data_b\n_cell_length_a 10\n_cell_length_b 10\n_cell_length_c 10\n_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\nloop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\nC1 C 0 0 0\nloop_\n_xtal_bond_atom_1\n_xtal_bond_atom_2\n0 0\n";
    let error = read_cif_stream(Cursor::new(bonds.as_bytes()))
        .expect_err("explicit bonds are not streamed");
    assert!(error.message().contains("explicit-bond"));

    let cartesian = "data_c\nloop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_Cartn_x\n_atom_site_Cartn_y\n_atom_site_Cartn_z\nC1 C 0 0 0\n";
    let error = read_cif_stream(Cursor::new(cartesian.as_bytes()))
        .expect_err("a Cartesian loop needs the cell first");
    assert!(error.message().contains("Cartesian"));
}

/// `BufRead` wrapper that counts the bytes handed to the reader.
struct CountingReader<R> {
    inner: R,
    consumed: usize,
}

impl<R> CountingReader<R> {
    fn new(inner: R) -> Self {
        Self { inner, consumed: 0 }
    }

    fn consumed(&self) -> usize {
        self.consumed
    }
}

impl<R: Read> Read for CountingReader<R> {
    fn read(&mut self, buffer: &mut [u8]) -> std::io::Result<usize> {
        let read = self.inner.read(buffer)?;
        self.consumed += read;
        Ok(read)
    }
}

impl<R: BufRead> BufRead for CountingReader<R> {
    fn fill_buf(&mut self) -> std::io::Result<&[u8]> {
        self.inner.fill_buf()
    }

    fn consume(&mut self, amount: usize) {
        self.consumed += amount;
        self.inner.consume(amount);
    }
}

#[test]
fn stream_frames_are_independent_of_the_file_length() {
    let text = xyz_text();
    let mut stream = XyzStream::new(Cursor::new(text.as_bytes()));
    let first: StreamFrame = stream.next_frame().unwrap().unwrap();
    assert_eq!(first.len(), ATOMS);
}
