//! Minimal PNG encoder (Task 8.1.6).
//!
//! The agent needs a raster picture without a native image library. This
//! encoder writes a truecolour (8-bit RGB, no alpha) PNG: signature, `IHDR`,
//! one zlib-compressed `IDAT` and `IEND`. `flate2` supplies the deflate
//! stream; the CRC is computed here.

use std::io::Write;
use std::sync::OnceLock;

use flate2::Compression;
use flate2::write::ZlibEncoder;

/// A truecolour image with 8 bits per channel.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RgbImage {
    width: u32,
    height: u32,
    pixels: Vec<u8>,
}

impl RgbImage {
    /// Creates a black image.
    ///
    /// # Panics
    ///
    /// Panics only when `width * height * 3` overflows `usize`, which needs a
    /// request above the documented image limit; the renderer checks its own
    /// dimensions first.
    #[must_use]
    pub fn new(width: u32, height: u32) -> Self {
        let len = (width as usize)
            .saturating_mul(height as usize)
            .saturating_mul(3);
        Self {
            width,
            height,
            pixels: vec![0; len],
        }
    }

    /// Creates an image filled with one colour.
    #[must_use]
    pub fn filled(width: u32, height: u32, color: [u8; 3]) -> Self {
        let mut image = Self::new(width, height);
        for pixel in image.pixels.chunks_mut(3) {
            if let [red, green, blue] = pixel {
                *red = color[0];
                *green = color[1];
                *blue = color[2];
            }
        }
        image
    }

    /// Width in pixels.
    #[must_use]
    pub fn width(&self) -> u32 {
        self.width
    }

    /// Height in pixels.
    #[must_use]
    pub fn height(&self) -> u32 {
        self.height
    }

    /// Raw RGB bytes, row major.
    #[must_use]
    pub fn pixels(&self) -> &[u8] {
        &self.pixels
    }

    /// Sets one pixel; out-of-range coordinates are ignored.
    pub fn set(&mut self, x: u32, y: u32, color: [u8; 3]) {
        if x >= self.width || y >= self.height {
            return;
        }
        let offset = (y as usize)
            .saturating_mul(self.width as usize)
            .saturating_add(x as usize)
            .saturating_mul(3);
        if let Some(pixel) = self.pixels.get_mut(offset..offset.saturating_add(3)) {
            pixel.copy_from_slice(&color);
        }
    }

    /// Encodes the image as PNG bytes.
    ///
    /// # Errors
    ///
    /// Returns [`std::io::Error`] when the deflate stream cannot be written;
    /// writing to an in-memory buffer does not fail in practice.
    pub fn to_png(&self) -> std::io::Result<Vec<u8>> {
        let mut raw = Vec::with_capacity(self.pixels.len() + self.height as usize);
        for row in self.pixels.chunks(self.width as usize * 3) {
            raw.push(0);
            raw.extend_from_slice(row);
        }
        let mut encoder = ZlibEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(&raw)?;
        let compressed = encoder.finish()?;

        let mut png = Vec::with_capacity(compressed.len() + 64);
        png.extend_from_slice(&[137, 80, 78, 71, 13, 10, 26, 10]);
        let mut ihdr = Vec::with_capacity(13);
        ihdr.extend_from_slice(&self.width.to_be_bytes());
        ihdr.extend_from_slice(&self.height.to_be_bytes());
        ihdr.extend_from_slice(&[8, 2, 0, 0, 0]);
        write_chunk(&mut png, b"IHDR", &ihdr);
        write_chunk(&mut png, b"IDAT", &compressed);
        write_chunk(&mut png, b"IEND", &[]);
        Ok(png)
    }
}

fn write_chunk(out: &mut Vec<u8>, kind: &[u8; 4], data: &[u8]) {
    out.extend_from_slice(&(data.len() as u32).to_be_bytes());
    out.extend_from_slice(kind);
    out.extend_from_slice(data);
    let mut crc = Crc32::new();
    crc.update(kind);
    crc.update(data);
    out.extend_from_slice(&crc.finish().to_be_bytes());
}

struct Crc32 {
    value: u32,
}

impl Crc32 {
    fn new() -> Self {
        Self { value: 0xFFFF_FFFF }
    }

    fn update(&mut self, bytes: &[u8]) {
        for byte in bytes {
            let index = ((self.value ^ u32::from(*byte)) & 0xFF) as usize;
            let entry = crc_table().get(index).copied().unwrap_or(0);
            self.value = entry ^ self.value.wrapping_shr(8);
        }
    }

    fn finish(&self) -> u32 {
        self.value ^ 0xFFFF_FFFF
    }
}

/// CRC-32 table (IEEE 802.3), built once on first use.
fn crc_table() -> &'static [u32; 256] {
    static TABLE: OnceLock<[u32; 256]> = OnceLock::new();
    TABLE.get_or_init(|| {
        let mut table = [0_u32; 256];
        for (index, slot) in table.iter_mut().enumerate() {
            let mut value = index as u32;
            for _ in 0..8 {
                value = if value & 1 == 1 {
                    0xEDB8_8320 ^ (value >> 1)
                } else {
                    value >> 1
                };
            }
            *slot = value;
        }
        table
    })
}
