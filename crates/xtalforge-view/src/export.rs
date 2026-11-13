//! Wavefront OBJ, glTF 2.0 and GLB export for the printable mesh (Task
//! 14.4.4).
//!
//! The exporters reuse [`build_indexed_mesh`], so they share vertices and
//! carry one averaged normal per vertex. Positions and normals are stored as
//! `f32` in glTF, which is the only float width the format allows; OBJ keeps
//! six decimals. The glTF output has one scene, one node, one mesh and one
//! primitive, with a single material-less triangle list: colours, materials
//! and multiple LODs are out of scope.

use serde_json::{Value, json};
use xtalforge_core::Error;

use crate::mesh::Triangle;
use crate::normals::build_indexed_mesh;

/// Component type of a `float` accessor in glTF 2.0.
const COMPONENT_FLOAT: u32 = 5126;
/// Component type of an `unsigned int` accessor in glTF 2.0.
const COMPONENT_UNSIGNED_INT: u32 = 5125;
/// `ARRAY_BUFFER` buffer-view target.
const TARGET_ARRAY_BUFFER: u32 = 34962;
/// `ELEMENT_ARRAY_BUFFER` buffer-view target.
const TARGET_ELEMENT_ARRAY_BUFFER: u32 = 34963;
/// `TRIANGLES` primitive mode.
const MODE_TRIANGLES: u32 = 4;

/// Writes a Wavefront OBJ with welded vertices and one vertex normal each.
///
/// The name is used as the object name after replacing whitespace and control
/// characters with underscores; an empty name becomes `xtalforge`. Faces use
/// the `v//vn` form, and the output ends with a newline even when the mesh is
/// empty.
#[must_use]
pub fn write_obj(triangles: &[Triangle], name: &str) -> String {
    let mesh = build_indexed_mesh(triangles);
    let mut out = String::new();
    out.push_str("# xtalforge mesh\n");
    out.push_str(&format!("o {}\n", object_name(name)));
    out.push_str(&format!(
        "# {} vertices, {} triangles\n",
        mesh.vertices.len(),
        mesh.triangle_count()
    ));
    for vertex in &mesh.vertices {
        out.push_str(&format!(
            "v {:.6} {:.6} {:.6}\n",
            vertex.position[0], vertex.position[1], vertex.position[2]
        ));
    }
    for vertex in &mesh.vertices {
        out.push_str(&format!(
            "vn {:.6} {:.6} {:.6}\n",
            vertex.normal[0], vertex.normal[1], vertex.normal[2]
        ));
    }
    for corners in &mesh.indices {
        let [first, second, third] = *corners;
        out.push_str(&format!(
            "f {}//{} {}//{} {}//{}\n",
            first + 1,
            first + 1,
            second + 1,
            second + 1,
            third + 1,
            third + 1
        ));
    }
    out
}

/// Writes glTF 2.0 JSON with the binary buffer embedded as a base64 data URI.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an empty mesh or for a coordinate that
/// is not a finite `f32`.
pub fn write_gltf(triangles: &[Triangle], name: &str) -> Result<String, Error> {
    let buffers = mesh_buffers(triangles)?;
    let bytes = buffer_bytes(&buffers);
    let document = gltf_document(&buffers, name, Some(data_uri(&bytes)));
    serde_json::to_string(&document).map_err(|error| {
        Error::InvalidInput(format!("could not serialize the glTF document: {error}"))
    })
}

/// Writes a binary glTF (`.glb`) with a JSON chunk and a `BIN` chunk.
///
/// Both chunks are padded to a four-byte boundary as the specification
/// requires; the JSON chunk is padded with spaces and the binary chunk with
/// zeros.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an empty or non-finite mesh, or when a
/// buffer would exceed the 32-bit limits of the container.
pub fn write_glb(triangles: &[Triangle], name: &str) -> Result<Vec<u8>, Error> {
    let buffers = mesh_buffers(triangles)?;
    let document = gltf_document(&buffers, name, None);
    let mut json_chunk = serde_json::to_vec(&document).map_err(|error| {
        Error::InvalidInput(format!("could not serialize the glTF document: {error}"))
    })?;
    while json_chunk.len() % 4 != 0 {
        json_chunk.push(b' ');
    }
    let mut binary_chunk = buffer_bytes(&buffers);
    while !binary_chunk.len().is_multiple_of(4) {
        binary_chunk.push(0);
    }
    let total = 12 + 8 + json_chunk.len() + 8 + binary_chunk.len();
    let total_length = to_u32(total, "the GLB file")?;
    let json_length = to_u32(json_chunk.len(), "the glTF JSON chunk")?;
    let binary_length = to_u32(binary_chunk.len(), "the glTF binary chunk")?;
    let mut out = Vec::with_capacity(total);
    out.extend_from_slice(b"glTF");
    out.extend_from_slice(&2u32.to_le_bytes());
    out.extend_from_slice(&total_length.to_le_bytes());
    out.extend_from_slice(&json_length.to_le_bytes());
    out.extend_from_slice(b"JSON");
    out.extend_from_slice(&json_chunk);
    out.extend_from_slice(&binary_length.to_le_bytes());
    out.extend_from_slice(b"BIN\0");
    out.extend_from_slice(&binary_chunk);
    Ok(out)
}

/// Interleaved-free buffers: positions, then normals, then indices.
struct MeshBuffers {
    positions: Vec<f32>,
    normals: Vec<f32>,
    indices: Vec<u32>,
    vertices: usize,
    minimum: [f64; 3],
    maximum: [f64; 3],
}

impl MeshBuffers {
    fn position_bytes(&self) -> usize {
        self.positions.len() * 4
    }

    fn normal_bytes(&self) -> usize {
        self.normals.len() * 4
    }

    fn index_bytes(&self) -> usize {
        self.indices.len() * 4
    }
}

fn mesh_buffers(triangles: &[Triangle]) -> Result<MeshBuffers, Error> {
    let mesh = build_indexed_mesh(triangles);
    if mesh.indices.is_empty() {
        return Err(Error::InvalidInput(
            "a glTF mesh needs at least one triangle".to_owned(),
        ));
    }
    let mut buffers = MeshBuffers {
        positions: Vec::with_capacity(mesh.vertices.len() * 3),
        normals: Vec::with_capacity(mesh.vertices.len() * 3),
        indices: Vec::with_capacity(mesh.indices.len() * 3),
        vertices: mesh.vertices.len(),
        minimum: [f64::INFINITY; 3],
        maximum: [f64::NEG_INFINITY; 3],
    };
    for vertex in &mesh.vertices {
        let narrowed = vertex.position.map(|value| value as f32);
        if !vertex.position.iter().all(|value| value.is_finite())
            || !narrowed.iter().all(|value| value.is_finite())
        {
            return Err(Error::InvalidInput(
                "mesh coordinates must be finite f32 values".to_owned(),
            ));
        }
        for (slot, value) in buffers.minimum.iter_mut().zip(narrowed.iter()) {
            *slot = (*slot).min(f64::from(*value));
        }
        for (slot, value) in buffers.maximum.iter_mut().zip(narrowed.iter()) {
            *slot = (*slot).max(f64::from(*value));
        }
        buffers.positions.extend(narrowed);
        let normal = vertex.normal.map(|value| value as f32);
        if !vertex.normal.iter().all(|value| value.is_finite())
            || !normal.iter().all(|value| value.is_finite())
        {
            return Err(Error::InvalidInput(
                "mesh normals must be finite f32 values".to_owned(),
            ));
        }
        buffers.normals.extend(normal);
    }
    for corners in &mesh.indices {
        for index in corners {
            let narrowed = u32::try_from(*index).map_err(|_| {
                Error::InvalidInput("the mesh has more vertices than glTF can address".to_owned())
            })?;
            buffers.indices.push(narrowed);
        }
    }
    Ok(buffers)
}

fn gltf_document(buffers: &MeshBuffers, name: &str, uri: Option<String>) -> Value {
    let position_bytes = buffers.position_bytes();
    let normal_bytes = buffers.normal_bytes();
    let index_bytes = buffers.index_bytes();
    let buffer_length = position_bytes + normal_bytes + index_bytes;
    let buffer = match uri {
        Some(uri) => json!({"byteLength": buffer_length, "uri": uri}),
        None => json!({"byteLength": buffer_length}),
    };
    json!({
        "asset": {
            "version": "2.0",
            "generator": format!("xtalforge {}", xtalforge_core::version()),
        },
        "scene": 0,
        "scenes": [{"name": name, "nodes": [0]}],
        "nodes": [{"name": name, "mesh": 0}],
        "meshes": [{
            "name": name,
            "primitives": [{
                "attributes": {"POSITION": 0, "NORMAL": 1},
                "indices": 2,
                "mode": MODE_TRIANGLES,
            }],
        }],
        "accessors": [
            {
                "bufferView": 0,
                "byteOffset": 0,
                "componentType": COMPONENT_FLOAT,
                "count": buffers.vertices,
                "type": "VEC3",
                "min": buffers.minimum,
                "max": buffers.maximum,
            },
            {
                "bufferView": 1,
                "byteOffset": 0,
                "componentType": COMPONENT_FLOAT,
                "count": buffers.vertices,
                "type": "VEC3",
            },
            {
                "bufferView": 2,
                "byteOffset": 0,
                "componentType": COMPONENT_UNSIGNED_INT,
                "count": buffers.indices.len(),
                "type": "SCALAR",
            },
        ],
        "bufferViews": [
            {
                "buffer": 0,
                "byteOffset": 0,
                "byteLength": position_bytes,
                "target": TARGET_ARRAY_BUFFER,
            },
            {
                "buffer": 0,
                "byteOffset": position_bytes,
                "byteLength": normal_bytes,
                "target": TARGET_ARRAY_BUFFER,
            },
            {
                "buffer": 0,
                "byteOffset": position_bytes + normal_bytes,
                "byteLength": index_bytes,
                "target": TARGET_ELEMENT_ARRAY_BUFFER,
            },
        ],
        "buffers": [buffer],
    })
}

fn buffer_bytes(buffers: &MeshBuffers) -> Vec<u8> {
    let mut bytes = Vec::with_capacity(
        (buffers.positions.len() + buffers.normals.len() + buffers.indices.len()) * 4,
    );
    for value in &buffers.positions {
        bytes.extend_from_slice(&value.to_le_bytes());
    }
    for value in &buffers.normals {
        bytes.extend_from_slice(&value.to_le_bytes());
    }
    for value in &buffers.indices {
        bytes.extend_from_slice(&value.to_le_bytes());
    }
    bytes
}

fn data_uri(bytes: &[u8]) -> String {
    format!(
        "data:application/octet-stream;base64,{}",
        base64_encode(bytes)
    )
}

fn object_name(name: &str) -> String {
    let cleaned: String = name
        .chars()
        .map(|character| {
            if character.is_whitespace() || character.is_control() {
                '_'
            } else {
                character
            }
        })
        .collect();
    if cleaned.is_empty() {
        "xtalforge".to_owned()
    } else {
        cleaned
    }
}

fn to_u32(value: usize, what: &str) -> Result<u32, Error> {
    u32::try_from(value)
        .map_err(|_| Error::InvalidInput(format!("{what} exceeds the glTF 32-bit limit")))
}

fn base64_encode(bytes: &[u8]) -> String {
    let mut out = String::with_capacity(bytes.len().div_ceil(3) * 4);
    for chunk in bytes.chunks(3) {
        let first = chunk.first().copied().unwrap_or(0);
        let second = chunk.get(1).copied();
        let third = chunk.get(2).copied();
        out.push(base64_character(first >> 2));
        out.push(base64_character(
            ((first & 0b11) << 4) | (second.unwrap_or(0) >> 4),
        ));
        match second {
            Some(second) => out.push(base64_character(
                ((second & 0b1111) << 2) | (third.unwrap_or(0) >> 6),
            )),
            None => out.push('='),
        }
        match third {
            Some(third) => out.push(base64_character(third & 0b11_1111)),
            None => out.push('='),
        }
    }
    out
}

fn base64_character(value: u8) -> char {
    match value {
        0..=25 => char::from(b'A' + value),
        26..=51 => char::from(b'a' + (value - 26)),
        52..=61 => char::from(b'0' + (value - 52)),
        62 => '+',
        _ => '/',
    }
}
