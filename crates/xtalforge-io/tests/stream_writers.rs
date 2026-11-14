//! Integration tests for the streaming writers and extended-XYZ columns
//! (Task 17.2.2).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::cell::RefCell;
use std::collections::BTreeMap;
use std::io::{self, Cursor, Write};
use std::rc::Rc;

use serde_json::Value;
use xtalforge_core::{Error, Lattice, Site, StructureBuilder};
use xtalforge_io::{
    LammpsStream, LammpsWriter, StreamFrame, XyzStream, XyzWriter, read_cif, read_cif_stream,
    write_cif, write_cif_stream, write_lammps_frame, write_xyz_frame,
};

fn lattice() -> Lattice {
    Lattice::new([[12.0, 0.0, 0.0], [1.5, 11.0, 0.0], [0.5, 1.0, 10.0]]).unwrap()
}

fn lattice_text() -> String {
    lattice()
        .matrix()
        .iter()
        .flatten()
        .map(|value| format!("{value:.10}"))
        .collect::<Vec<_>>()
        .join(" ")
}

fn frac(index: usize, frame: u64) -> [f64; 3] {
    [
        (index as f64 * 0.125 + frame as f64 * 0.0625).rem_euclid(1.0),
        (index as f64 * 0.25 + frame as f64 * 0.125).rem_euclid(1.0),
        (index as f64 * 0.5 + frame as f64 * 0.25).rem_euclid(1.0),
    ]
}

fn frame(step: u64) -> StreamFrame {
    StreamFrame {
        elements: vec!["C".into(), "O".into(), "H".into(), "Zn".into()],
        frac: (0..4).map(|index| frac(index, step)).collect(),
        lattice: Some(lattice().matrix()),
        energy: Some(-12.5 - step as f64),
        step: Some(step),
        info: [("source".to_owned(), "unit test".to_owned())]
            .into_iter()
            .collect(),
        forces: None,
        columns: BTreeMap::new(),
    }
}

fn text_of(mut write: impl FnMut(&mut Vec<u8>)) -> String {
    let mut buffer = Vec::new();
    write(&mut buffer);
    String::from_utf8(buffer).unwrap()
}

#[test]
fn xyz_round_trips_frames_through_the_stream_reader() {
    let frames = vec![frame(0), frame(3)];
    let mut writer = XyzWriter::new(Vec::new());
    assert_eq!(writer.write_frames(frames.clone()).unwrap(), 2);
    let text = String::from_utf8(writer.into_inner()).unwrap();

    assert!(text.contains("Lattice=\"12.0000000000"));
    assert!(text.contains("energy=-12.5000000000"));
    assert!(text.contains("step=3"));
    assert!(text.contains("source=\"unit test\""));

    let mut stream = XyzStream::new(Cursor::new(text.as_bytes()));
    for expected in &frames {
        let actual = stream.next_frame().unwrap().expect("one frame");
        assert_eq!(actual.elements, expected.elements);
        assert_eq!(actual.step, expected.step);
        assert_eq!(actual.energy, expected.energy);
        assert_eq!(actual.lattice, expected.lattice);
        assert_eq!(
            actual.info.get("source").map(String::as_str),
            Some("unit test")
        );
        assert!(actual.forces.is_none());
        assert!(actual.columns.is_empty());
        for (actual, expected) in actual.frac.iter().zip(expected.frac.iter()) {
            for axis in 0..3 {
                assert!((actual[axis] - expected[axis]).abs() < 1e-9);
            }
        }
    }
    assert!(stream.next_frame().unwrap().is_none());
}

#[test]
fn xyz_round_trips_forces_and_custom_columns() {
    let mut expected = frame(4);
    expected.forces = Some(vec![
        [0.25, -0.5, 1.0],
        [0.0, 0.0, 0.0],
        [-1.0, 2.0, -3.0],
        [0.125, -0.25, 0.5],
    ]);
    expected.columns.insert(
        "charge".into(),
        vec![1.0.into(), 2.0.into(), 3.0.into(), 4.0.into()],
    );
    expected.columns.insert(
        "serial".into(),
        vec![7_i64.into(), 8_i64.into(), 9_i64.into(), 10_i64.into()],
    );
    expected.columns.insert(
        "fixed".into(),
        vec![
            Value::Bool(true),
            Value::Bool(false),
            Value::Bool(true),
            Value::Bool(true),
        ],
    );
    expected.columns.insert(
        "kind".into(),
        vec![
            "sp2 hybrid".into(),
            "sp3".into(),
            "terminal".into(),
            "metal".into(),
        ],
    );
    expected.columns.insert(
        "shift".into(),
        vec![
            Value::Array(vec![1.5.into(), 2.5.into()]),
            Value::Array(vec![0.0.into(), 0.0.into()]),
            Value::Array(vec![Value::from(-1.0), Value::from(3.0)]),
            Value::Array(vec![Value::from(4.0), Value::from(-2.0)]),
        ],
    );

    let text = text_of(|buffer| write_xyz_frame(buffer, &expected).unwrap());
    assert!(text.contains(
        "Properties=species:S:1:pos:R:3:forces:R:3:charge:R:1:fixed:L:1:kind:S:1:serial:I:1:shift:R:2"
    ));
    assert!(text.contains("\"sp2 hybrid\""));

    let actual = XyzStream::new(Cursor::new(text.as_bytes()))
        .next_frame()
        .unwrap()
        .expect("one frame");
    assert_eq!(actual.elements, expected.elements);
    assert_eq!(actual.forces, expected.forces);
    assert_eq!(actual.columns, expected.columns);
    assert_eq!(actual.step, expected.step);
    assert_eq!(actual.energy, expected.energy);
}

#[test]
fn xyz_reader_accepts_separate_fx_fy_fz_columns() {
    let text = format!(
        "1\nLattice=\"{}\" Properties=species:S:1:pos:R:3:fx:R:1:fy:R:1:fz:R:1 energy=-1.0 step=2\n\
         C 1.0 2.0 3.0 0.25 -0.5 1.0\n",
        lattice_text()
    );
    let frame = XyzStream::new(Cursor::new(text.as_bytes()))
        .next_frame()
        .unwrap()
        .expect("one frame");
    assert_eq!(frame.forces, Some(vec![[0.25, -0.5, 1.0]]));
    assert!(frame.columns.is_empty());
    assert_eq!(frame.step, Some(2));
    assert_eq!(frame.energy, Some(-1.0));
}

#[test]
fn xyz_reader_rejects_declarations_that_do_not_fit_the_rows() {
    for text in [
        "1\nProperties=species:S:1 energy=0\nC 1 2 3\n",
        "1\nProperties=species:S:1:pos:R:3:forces:R:3\nC 1 2 3 0.1 0.2\n",
        "1\nProperties=species:S:1:pos:R:3:bogus:Q:1\nC 1 2 3 0\n",
        "1\nProperties=species:S:1:pos:R:3:forces:R:2\nC 1 2 3 0.1 0.2\n",
    ] {
        let error = XyzStream::new(Cursor::new(text.as_bytes()))
            .next_frame()
            .expect_err("malformed columns are an error");
        assert_eq!(error.code(), "invalid_input");
    }
}

#[test]
fn xyz_writer_rejects_malformed_frames() {
    let good = frame(0);

    let mut short = good.clone();
    short.frac.pop();
    assert!(write_xyz_frame(&mut Vec::new(), &short).is_err());

    let mut non_finite = good.clone();
    non_finite.frac[0][0] = f64::NAN;
    assert!(write_xyz_frame(&mut Vec::new(), &non_finite).is_err());

    let mut non_finite_energy = good.clone();
    non_finite_energy.energy = Some(f64::INFINITY);
    assert!(write_xyz_frame(&mut Vec::new(), &non_finite_energy).is_err());

    let mut unknown = good.clone();
    unknown.elements[0] = "Qq".into();
    assert!(write_xyz_frame(&mut Vec::new(), &unknown).is_err());

    let mut short_forces = good.clone();
    short_forces.forces = Some(vec![[0.0, 0.0, 0.0]]);
    assert!(write_xyz_frame(&mut Vec::new(), &short_forces).is_err());

    let mut non_finite_force = good.clone();
    non_finite_force.forces = Some(vec![[f64::NAN, 0.0, 0.0]; 4]);
    assert!(write_xyz_frame(&mut Vec::new(), &non_finite_force).is_err());

    let mut short_column = good.clone();
    short_column
        .columns
        .insert("charge".into(), vec![1.0.into()]);
    assert!(write_xyz_frame(&mut Vec::new(), &short_column).is_err());

    let mut mixed_column = good.clone();
    mixed_column.columns.insert(
        "mixed".into(),
        vec![1.0.into(), "two".into(), 3.0.into(), 4.0.into()],
    );
    assert!(write_xyz_frame(&mut Vec::new(), &mixed_column).is_err());

    let mut reserved = good.clone();
    reserved.columns.insert(
        "forces".into(),
        vec![1.0.into(), 2.0.into(), 3.0.into(), 4.0.into()],
    );
    assert!(write_xyz_frame(&mut Vec::new(), &reserved).is_err());

    let empty = StreamFrame {
        lattice: Some(lattice().matrix()),
        ..StreamFrame::default()
    };
    assert!(write_xyz_frame(&mut Vec::new(), &empty).is_err());
}

#[test]
fn lammps_round_trips_and_is_deterministic() {
    let frames = vec![frame(0), frame(2)];
    let text = {
        let mut writer = LammpsWriter::new(Vec::new());
        assert_eq!(writer.write_frames(frames.clone()).unwrap(), 2);
        String::from_utf8(writer.into_inner()).unwrap()
    };
    let again = {
        let mut writer = LammpsWriter::new(Vec::new());
        writer.write_frames(frames.clone()).unwrap();
        String::from_utf8(writer.into_inner()).unwrap()
    };
    assert_eq!(text, again);
    assert!(text.contains("ITEM: BOX BOUNDS xy xz yz pp pp pp"));
    assert!(text.contains("ITEM: ATOMS id element x y z"));

    let mut stream = LammpsStream::new(Cursor::new(text.as_bytes()));
    for expected in &frames {
        let actual = stream.next_frame().unwrap().expect("one frame");
        assert_eq!(actual.step, expected.step);
        assert_eq!(actual.elements, expected.elements);
        let actual_lattice = Lattice::new(actual.lattice.expect("cell")).unwrap();
        for (actual, expected) in actual_lattice
            .parameters()
            .iter()
            .zip(lattice().parameters().iter())
        {
            assert!(
                (actual - expected).abs() < 1e-9,
                "cell parameter {actual} differs from {expected}"
            );
        }
        for (actual, expected) in actual.frac.iter().zip(expected.frac.iter()) {
            for axis in 0..3 {
                assert!((actual[axis] - expected[axis]).abs() < 1e-9);
            }
        }
    }
    assert!(stream.next_frame().unwrap().is_none());
}

#[test]
fn xyz_writer_accepts_plain_frames_without_a_cell() {
    let mut plain = frame(0);
    plain.lattice = None;
    plain.energy = None;
    plain.step = None;
    plain.info.clear();
    let text = text_of(|buffer| write_xyz_frame(buffer, &plain).unwrap());
    assert!(!text.contains("Lattice="));
    let actual = XyzStream::new(Cursor::new(text.as_bytes()))
        .next_frame()
        .unwrap()
        .expect("one frame");
    assert!(actual.lattice.is_none());
    assert_eq!(actual.frac, plain.frac);
}

#[test]
fn lammps_round_trips_a_cell_with_negative_tilt() {
    let cell = Lattice::new([[10.0, 0.0, 0.0], [-2.0, 9.0, 0.0], [1.0, -1.0, 8.0]]).unwrap();
    let mut input = frame(0);
    input.lattice = Some(cell.matrix());
    let text = text_of(|buffer| write_lammps_frame(buffer, &input).unwrap());
    let actual = LammpsStream::new(Cursor::new(text.as_bytes()))
        .next_frame()
        .unwrap()
        .expect("one frame");
    let back = Lattice::new(actual.lattice.expect("cell")).unwrap();
    for (actual, expected) in back.parameters().iter().zip(cell.parameters().iter()) {
        assert!(
            (actual - expected).abs() < 1e-9,
            "cell parameter {actual} differs from {expected}"
        );
    }
}

#[test]
fn lammps_writer_needs_a_lattice_and_a_known_element() {
    let mut plain = frame(0);
    plain.lattice = None;
    let error = write_lammps_frame(&mut Vec::new(), &plain).expect_err("no cell");
    assert!(error.message().contains("lattice"));

    let mut unknown = frame(0);
    unknown.elements[1] = "Qq".into();
    assert!(write_lammps_frame(&mut Vec::new(), &unknown).is_err());
}

#[test]
fn writers_propagate_source_and_io_errors() {
    let error = XyzWriter::new(Vec::new())
        .write_results(vec![Err(Error::InvalidInput("source failed".into()))])
        .expect_err("source error");
    assert_eq!(error.code(), "invalid_input");

    let error = write_xyz_frame(&mut FailingWriter, &frame(0)).expect_err("sink error");
    assert!(error.message().contains("could not write"));
    let error = write_lammps_frame(&mut FailingWriter, &frame(0)).expect_err("sink error");
    assert!(error.message().contains("could not write"));
}

#[test]
fn cif_stream_writer_matches_the_string_writer_and_round_trips() {
    let structure = StructureBuilder::new("stream", Lattice::cubic(4.0).unwrap())
        .site(
            Site::new("C1", "C", [0.1, 0.2, 0.3])
                .unwrap()
                .with_occupancy(0.5)
                .unwrap()
                .with_u_iso(0.01)
                .unwrap()
                .with_charge(-1.0)
                .unwrap(),
        )
        .site(Site::new("O1", "O", [0.4, 0.5, 0.6]).unwrap())
        .build()
        .unwrap();
    let expected = write_cif(&structure).unwrap();
    let mut buffer = Vec::new();
    write_cif_stream(&mut buffer, &structure).unwrap();
    assert_eq!(&buffer, expected.as_bytes());

    let text = String::from_utf8(buffer).unwrap();
    let back = read_cif(&text).unwrap();
    assert_eq!(back.sites().len(), 2);
    assert_eq!(back.sites()[0].occupancy(), 0.5);
    assert_eq!(back.sites()[0].charge(), Some(-1.0));
    assert_eq!(back.sites()[0].u_iso(), Some(0.01));

    let streamed = read_cif_stream(Cursor::new(text.as_bytes())).unwrap();
    assert_eq!(streamed.sites().len(), 2);
    assert_eq!(streamed.hall(), structure.hall());
    assert_eq!(streamed.lattice().parameters()[0], 4.0);

    let error = write_cif_stream(&mut FailingWriter, &structure).expect_err("sink error");
    assert!(error.message().contains("could not write the CIF"));
}

#[test]
fn cif_stream_writer_matches_the_corpus_reader() {
    let structure = read_cif(include_str!("../../../docs/public/samples/Si.cif")).unwrap();
    let mut buffer = Vec::new();
    write_cif_stream(&mut buffer, &structure).unwrap();
    assert_eq!(buffer, write_cif(&structure).unwrap().as_bytes());
    let back = read_cif_stream(Cursor::new(buffer.as_slice())).unwrap();
    assert_eq!(back.sites().len(), structure.sites().len());
    assert_eq!(
        back.lattice().parameters(),
        structure.lattice().parameters()
    );
}

#[test]
fn xyz_writer_consumes_lazy_frames_one_at_a_time() {
    // The million-atom claim is structural: `write_frames` consumes this
    // `FromFn` iterator one frame at a time, so peak memory is one frame and
    // not the whole trajectory.
    let sink = SharedSink::default();
    let mut writer = XyzWriter::new(sink.clone());
    let mut produced = 0_u64;
    let frames = std::iter::from_fn(|| {
        if produced == 4 {
            return None;
        }
        if produced > 0 {
            let written = String::from_utf8(sink.0.borrow().clone()).unwrap();
            assert_eq!(
                written.matches("Properties=").count(),
                produced as usize,
                "the writer asked for the next frame before writing the previous one"
            );
        }
        produced += 1;
        Some(frame(produced - 1))
    });
    assert_eq!(writer.write_frames(frames).unwrap(), 4);
}

#[test]
fn writer_accepts_a_reader_pipeline_without_collecting_frames() {
    let source = text_of(|buffer| {
        let mut writer = XyzWriter::new(buffer);
        writer.write_frames(vec![frame(0), frame(1)]).unwrap();
    });
    let mut stream = XyzStream::new(Cursor::new(source.as_bytes()));
    let mut writer = XyzWriter::new(Vec::new());
    let count = writer
        .write_results(std::iter::from_fn(|| stream.next_frame().transpose()))
        .unwrap();
    assert_eq!(count, 2);
    let text = String::from_utf8(writer.into_inner()).unwrap();
    assert_eq!(text.matches("Lattice=").count(), 2);
}

/// Sink that always fails, for the writer error paths.
struct FailingWriter;

impl Write for FailingWriter {
    fn write(&mut self, _buffer: &[u8]) -> io::Result<usize> {
        Err(io::Error::other("disk full"))
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

/// Sink shared with a lazy frame iterator.
#[derive(Clone, Default)]
struct SharedSink(Rc<RefCell<Vec<u8>>>);

impl Write for SharedSink {
    fn write(&mut self, buffer: &[u8]) -> io::Result<usize> {
        self.0.borrow_mut().extend_from_slice(buffer);
        Ok(buffer.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
