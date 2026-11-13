//! Vertex normals and OBJ/glTF/GLB export (Task 14.4.4).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use serde_json::Value;
use xtalforge_core::{Lattice, Site, StructureBuilder};
use xtalforge_view::{
    Triangle, ViewSettings, build_indexed_mesh, build_mesh, build_scene, vertex_normals, write_glb,
    write_gltf, write_obj, write_stl,
};

fn flat_triangle(first: [f64; 3], second: [f64; 3], third: [f64; 3]) -> Triangle {
    Triangle {
        normal: [0.0, 0.0, 1.0],
        vertices: [first, second, third],
    }
}

/// Two coplanar triangles that share an edge: a unit square.
fn quad() -> Vec<Triangle> {
    vec![
        flat_triangle([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]),
        flat_triangle([0.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]),
    ]
}

fn diamond() -> xtalforge_core::Structure {
    StructureBuilder::new("diamond", Lattice::cubic(3.566_986).unwrap())
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("C2", "C", [0.25, 0.25, 0.25]).unwrap())
        .build()
        .unwrap()
}

fn decode_character(character: u8) -> u8 {
    match character {
        b'A'..=b'Z' => character - b'A',
        b'a'..=b'z' => character - b'a' + 26,
        b'0'..=b'9' => character - b'0' + 52,
        b'+' => 62,
        b'/' => 63,
        b'=' => 0,
        _ => panic!("invalid base64 character"),
    }
}

fn base64_decode(text: &str) -> Vec<u8> {
    let mut out = Vec::new();
    for chunk in text.as_bytes().chunks(4) {
        let first = decode_character(chunk[0]);
        let second = decode_character(chunk[1]);
        let third = decode_character(chunk[2]);
        let fourth = decode_character(chunk[3]);
        out.push((first << 2) | (second >> 4));
        if chunk[2] != b'=' {
            out.push((second << 4) | (third >> 2));
        }
        if chunk[3] != b'=' {
            out.push((third << 6) | fourth);
        }
    }
    out
}

fn read_u32(bytes: &[u8], offset: usize) -> u32 {
    let mut array = [0u8; 4];
    array.copy_from_slice(&bytes[offset..offset + 4]);
    u32::from_le_bytes(array)
}

fn read_f32(bytes: &[u8], offset: usize) -> f32 {
    let mut array = [0u8; 4];
    array.copy_from_slice(&bytes[offset..offset + 4]);
    f32::from_le_bytes(array)
}

fn length(vector: [f64; 3]) -> f64 {
    (vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]).sqrt()
}

#[test]
fn vertex_normals_are_welded_and_unit_length() {
    let triangles = quad();
    let normals = vertex_normals(&triangles);
    assert_eq!(normals.len(), 2);
    for corners in &normals {
        for normal in corners {
            assert!((length(*normal) - 1.0).abs() < 1e-12);
            assert!(normal[0].abs() < 1e-12);
            assert!(normal[1].abs() < 1e-12);
            assert!((normal[2] - 1.0).abs() < 1e-12);
        }
    }
    let mesh = build_indexed_mesh(&triangles);
    assert_eq!(mesh.vertices.len(), 4);
    assert_eq!(mesh.indices.len(), 2);
    assert_eq!(mesh.triangle_count(), 2);
    assert!(!mesh.is_empty());
    assert_eq!(mesh.indices[0][0], mesh.indices[1][0]);
    for vertex in &mesh.vertices {
        assert!((length(vertex.normal) - 1.0).abs() < 1e-12);
    }
    let empty = build_indexed_mesh(&[]);
    assert!(empty.is_empty());
    assert_eq!(empty.triangle_count(), 0);
}

#[test]
fn obj_parses_back_with_welded_vertices_and_faces() {
    let text = write_obj(&quad(), "my mesh");
    assert!(text.starts_with("# xtalforge mesh\n"));
    assert!(text.contains("o my_mesh\n"));
    assert!(text.ends_with('\n'));
    let mut vertices = 0usize;
    let mut normals = 0usize;
    let mut faces = 0usize;
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("v ") {
            let values: Vec<f64> = rest
                .split(' ')
                .map(|value| value.parse().unwrap())
                .collect();
            assert_eq!(values.len(), 3);
            vertices += 1;
        } else if let Some(rest) = line.strip_prefix("vn ") {
            let values: Vec<f64> = rest
                .split(' ')
                .map(|value| value.parse().unwrap())
                .collect();
            assert_eq!(values.len(), 3);
            assert!((length([values[0], values[1], values[2]]) - 1.0).abs() < 1e-6);
            normals += 1;
        } else if let Some(rest) = line.strip_prefix("f ") {
            let corners: Vec<&str> = rest.split(' ').collect();
            assert_eq!(corners.len(), 3);
            for corner in corners {
                let (vertex, normal) = corner.split_once("//").unwrap();
                let vertex: usize = vertex.parse().unwrap();
                let normal: usize = normal.parse().unwrap();
                assert!(vertex >= 1 && vertex <= vertices);
                assert!(normal >= 1 && normal <= normals);
            }
            faces += 1;
        }
    }
    assert_eq!(vertices, 4);
    assert_eq!(normals, 4);
    assert_eq!(faces, 2);
}

#[test]
fn obj_handles_empty_names_and_empty_meshes() {
    let text = write_obj(&[], "");
    assert!(text.contains("o xtalforge\n"));
    assert!(!text.contains("\nv "));
    assert!(!text.contains("\nf "));
}

#[test]
fn gltf_parses_back_with_an_embedded_buffer() {
    let text = write_gltf(&quad(), "quad").unwrap();
    let document: Value = serde_json::from_str(&text).unwrap();
    assert_eq!(document["asset"]["version"], "2.0");
    assert!(
        document["asset"]["generator"]
            .as_str()
            .unwrap()
            .starts_with("xtalforge ")
    );
    assert_eq!(document["meshes"][0]["primitives"][0]["mode"], 4);
    assert_eq!(document["meshes"][0]["primitives"][0]["indices"], 2);
    let accessors = document["accessors"].as_array().unwrap();
    assert_eq!(accessors[0]["count"], 4);
    assert_eq!(accessors[0]["type"], "VEC3");
    assert_eq!(accessors[0]["componentType"], 5126);
    assert_eq!(accessors[0]["min"], serde_json::json!([0.0, 0.0, 0.0]));
    assert_eq!(accessors[0]["max"], serde_json::json!([1.0, 1.0, 0.0]));
    assert_eq!(accessors[2]["count"], 6);
    assert_eq!(accessors[2]["componentType"], 5125);
    let views = document["bufferViews"].as_array().unwrap();
    assert_eq!(views[0]["byteLength"], 4 * 3 * 4);
    assert_eq!(views[1]["byteOffset"], 4 * 3 * 4);
    assert_eq!(views[2]["byteOffset"], 4 * 3 * 4 * 2);
    assert_eq!(views[2]["byteLength"], 6 * 4);
    let uri = document["buffers"][0]["uri"].as_str().unwrap();
    let encoded = uri
        .strip_prefix("data:application/octet-stream;base64,")
        .unwrap();
    let bytes = base64_decode(encoded);
    assert_eq!(
        bytes.len(),
        document["buffers"][0]["byteLength"].as_u64().unwrap() as usize
    );
    let positions = document["bufferViews"][0]["byteOffset"].as_u64().unwrap() as usize;
    assert_eq!(read_f32(&bytes, positions), 0.0);
    assert_eq!(read_f32(&bytes, positions + 4 * 3), 1.0);
    let normals = document["bufferViews"][1]["byteOffset"].as_u64().unwrap() as usize;
    for corner in 0..4 {
        assert_eq!(read_f32(&bytes, normals + corner * 12 + 8), 1.0);
    }
    let indices = document["bufferViews"][2]["byteOffset"].as_u64().unwrap() as usize;
    for corner in 0..6 {
        assert!(read_u32(&bytes, indices + corner * 4) < 4);
    }
}

#[test]
fn glb_parses_back_with_padded_chunks() {
    let bytes = write_glb(&quad(), "quad").unwrap();
    assert_eq!(&bytes[0..4], b"glTF");
    assert_eq!(read_u32(&bytes, 4), 2);
    assert_eq!(read_u32(&bytes, 8) as usize, bytes.len());
    assert_eq!(bytes.len() % 4, 0);
    let json_length = read_u32(&bytes, 12) as usize;
    assert_eq!(&bytes[16..20], b"JSON");
    assert_eq!(json_length % 4, 0);
    let json_end = 20 + json_length;
    let document: Value = serde_json::from_slice(&bytes[20..json_end]).unwrap();
    assert!(document["buffers"][0].get("uri").is_none());
    let binary_length = read_u32(&bytes, json_end) as usize;
    assert_eq!(&bytes[json_end + 4..json_end + 8], b"BIN\0");
    assert_eq!(binary_length % 4, 0);
    let binary_start = json_end + 8;
    assert_eq!(binary_start + binary_length, bytes.len());
    let declared = document["buffers"][0]["byteLength"].as_u64().unwrap() as usize;
    assert!(declared <= binary_length);
    assert!(binary_length - declared < 4);
    assert_eq!(read_f32(&bytes, binary_start), 0.0);
    let normals = binary_start + 4 * 3 * 4;
    assert_eq!(read_f32(&bytes, normals + 8), 1.0);
}

#[test]
fn names_are_escaped_in_gltf_json() {
    let name = "quad \"quoted\" \\ slash <angle>";
    let text = write_gltf(&quad(), name).unwrap();
    let document: Value = serde_json::from_str(&text).unwrap();
    assert_eq!(document["nodes"][0]["name"].as_str().unwrap(), name);
    assert_eq!(document["meshes"][0]["name"].as_str().unwrap(), name);
    assert_eq!(document["scenes"][0]["name"].as_str().unwrap(), name);
}

#[test]
fn empty_and_non_finite_meshes_are_rejected() {
    assert!(write_gltf(&[], "empty").is_err());
    assert!(write_glb(&[], "empty").is_err());
    let broken = vec![flat_triangle(
        [f64::NAN, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
    )];
    assert!(write_gltf(&broken, "broken").is_err());
    assert!(write_glb(&broken, "broken").is_err());
    let infinite = vec![flat_triangle(
        [f64::INFINITY, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
    )];
    assert!(write_gltf(&infinite, "infinite").is_err());
}

#[test]
fn a_scene_mesh_exports_consistently_to_every_format() {
    let scene = build_scene(&diamond(), &ViewSettings::default()).unwrap();
    let triangles = build_mesh(&scene);
    let mesh = build_indexed_mesh(&triangles);
    assert_eq!(mesh.indices.len(), triangles.len());
    assert!(mesh.vertices.len() < triangles.len() * 3);
    for vertex in &mesh.vertices {
        assert!((length(vertex.normal) - 1.0).abs() < 1e-9);
    }
    let obj = write_obj(&triangles, "diamond");
    assert_eq!(obj.matches("\nf ").count(), triangles.len());
    let stl = write_stl(&triangles, "diamond");
    assert_eq!(stl.matches("facet normal").count(), triangles.len());
    let gltf = write_gltf(&triangles, "diamond").unwrap();
    let document: Value = serde_json::from_str(&gltf).unwrap();
    assert_eq!(
        document["accessors"][0]["count"].as_u64().unwrap() as usize,
        mesh.vertices.len()
    );
    assert_eq!(
        document["accessors"][2]["count"].as_u64().unwrap() as usize,
        triangles.len() * 3
    );
    assert_eq!(write_gltf(&triangles, "diamond").unwrap(), gltf);
}
