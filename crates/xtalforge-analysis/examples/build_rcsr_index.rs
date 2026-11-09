//! Builds the embedded RCSR fingerprint index (Task 5.6.4).
//!
//! Input is the decompressed Crystal-Builder `rcsr-2019-06-01.json`
//! (`gzip -dc rcsr-2019-06-01.json.gz > /tmp/rcsr.json`); output is the
//! compact binary in `crates/xtalforge-analysis/data/rcsr_index.bin`.
//!
//! Every keyed catalogue entry is rebuilt from its canonical key and its
//! fingerprint is compared with the row before the entry is written. The
//! build fails on any mismatch, so the index cannot silently drift from the
//! data it describes.

use std::collections::BTreeMap;
use std::path::PathBuf;

use xtalforge_analysis::{NetFingerprint, QuotientGraph, fingerprint};

const MAGIC: &[u8; 8] = b"XFRCS01\n";

type Failure = Box<dyn std::error::Error>;

fn main() -> Result<(), Failure> {
    let mut arguments = std::env::args().skip(1);
    let input = arguments
        .next()
        .map_or_else(|| PathBuf::from("/tmp/rcsr.json"), PathBuf::from);
    let output = arguments.next().map_or_else(
        || PathBuf::from("crates/xtalforge-analysis/data/rcsr_index.bin"),
        PathBuf::from,
    );
    let text = std::fs::read_to_string(&input)?;
    let payload: serde_json::Value = serde_json::from_str(&text)?;
    let rows = payload
        .get("nets")
        .and_then(serde_json::Value::as_array)
        .ok_or("the RCSR JSON has no nets array")?;
    let mut names = BTreeMap::new();
    let mut entries = Vec::with_capacity(rows.len());
    let mut keyed = 0usize;
    for row in rows {
        let name = row[0].as_str().ok_or("a net has no name")?.to_owned();
        if names.insert(name.clone(), ()).is_some() {
            return Err(format!("duplicate RCSR name {name}").into());
        }
        let periodicity = row[1].as_u64().ok_or("a net has no periodicity")? as usize;
        let mut sequences: Vec<Vec<usize>> = row[2]
            .as_array()
            .ok_or("a net has no sequences")?
            .iter()
            .map(|sequence| {
                sequence
                    .as_array()
                    .ok_or("a sequence is not an array")
                    .map(|values| {
                        values
                            .iter()
                            .map(|value| value.as_u64().unwrap_or(0) as usize)
                            .collect()
                    })
            })
            .collect::<Result<_, _>>()?;
        sequences.sort();
        let mut symbols: Vec<String> = row
            .get(3)
            .and_then(serde_json::Value::as_array)
            .ok_or("a net has no symbols")?
            .iter()
            .map(|value| value.as_str().unwrap_or("").to_owned())
            .collect();
        symbols.sort();
        let expected = NetFingerprint {
            periodicity,
            sequences,
            symbols,
        };
        let key = row.get(7).and_then(|value| value.as_str()).unwrap_or("");
        let actual = if key.is_empty() {
            expected.clone()
        } else {
            keyed += 1;
            let graph = QuotientGraph::from_canonical_key(key)?;
            let found = fingerprint(&graph)?;
            if found != expected {
                return Err(format!("{name}: fingerprint does not match its canonical key").into());
            }
            found
        };
        entries.push((name, actual));
    }
    // Find a salt that gives every distinct fingerprint token a distinct
    // 32-bit hash. Shared tokens (ambiguous nets) are allowed to share.
    let mut salt = 0u32;
    loop {
        let mut seen: BTreeMap<u32, String> = BTreeMap::new();
        let mut ok = true;
        for (_name, fingerprint) in &entries {
            let token = fingerprint.token();
            let hash = fingerprint.hash_with_salt(salt);
            if let Some(other) = seen.get(&hash) {
                if *other != token {
                    ok = false;
                    break;
                }
            } else {
                seen.insert(hash, token.clone());
            }
        }
        if ok {
            break;
        }
        salt = salt
            .checked_add(1)
            .ok_or("no hash salt avoids a collision")?;
    }
    let mut bytes = Vec::new();
    bytes.extend_from_slice(MAGIC);
    let count = u32::try_from(entries.len())?;
    bytes.extend_from_slice(&count.to_le_bytes());
    bytes.extend_from_slice(&salt.to_le_bytes());
    let mut hashes: BTreeMap<u32, String> = BTreeMap::new();
    for (name, fingerprint) in &entries {
        let hash = fingerprint.hash_with_salt(salt);
        hashes.insert(hash, name.clone());
        let name_bytes = name.as_bytes();
        let name_len = u8::try_from(name_bytes.len())?;
        bytes.push(name_len);
        bytes.extend_from_slice(name_bytes);
        bytes.push(u8::try_from(fingerprint.periodicity)?);
        let coordination = fingerprint.coordination();
        bytes.push(u8::try_from(coordination.len())?);
        for value in &coordination {
            bytes.push(u8::try_from(*value)?);
        }
        bytes.extend_from_slice(&hash.to_le_bytes());
    }
    std::fs::write(&output, &bytes)?;
    println!(
        "wrote {} entries ({} keyed, {} unkeyed, {} shared hashes, salt {}) to {} ({} bytes)",
        entries.len(),
        keyed,
        entries.len() - keyed,
        entries.len() - hashes.len(),
        salt,
        output.display(),
        bytes.len()
    );
    Ok(())
}
