//! Change tracking (Task 8.1.4).
//!
//! The document keeps one [`CommandRecord`] per edit with its
//! [`Origin`](xtalforge_commands::Origin). This module groups consecutive
//! edits made through the agent into bursts, so a user can review the
//! assistant's work and undo it as a unit.
//!
//! Ported from Crystal-Builder `xtal/agent/changes.py` (MIT).

use serde::{Deserialize, Serialize};
use xtalforge_commands::{CommandRecord, Origin};

/// A consecutive run of edits made through the agent.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ChangeBurst {
    /// Position of the first edit in the document history, zero based.
    pub start: usize,
    /// Number of edits in the burst.
    pub len: usize,
    /// One description per edit, in application order.
    pub descriptions: Vec<String>,
    /// Number of `undo` calls that remove the whole burst when it is the
    /// most recent history entry; `None` when later user edits exist.
    pub undo_steps: Option<usize>,
}

/// Groups the agent-origin runs of a command log.
///
/// A run is consecutive records whose origin is [`Origin::Agent`]; a user or
/// script edit between two agent edits starts a new burst. The log only holds
/// entries that were not trimmed from the undo stack.
#[must_use]
pub fn changes_for(history: &[CommandRecord]) -> Vec<ChangeBurst> {
    let mut bursts = Vec::new();
    let mut index = 0;
    while let Some(record) = history.get(index) {
        if record.origin != Origin::Agent {
            index += 1;
            continue;
        }
        let start = index;
        while let Some(record) = history.get(index) {
            if record.origin != Origin::Agent {
                break;
            }
            index += 1;
        }
        let Some(run) = history.get(start..index) else {
            break;
        };
        let len = index - start;
        let descriptions = run
            .iter()
            .map(|record| record.description.clone())
            .collect();
        let at_end = index == history.len();
        bursts.push(ChangeBurst {
            start,
            len,
            descriptions,
            undo_steps: if at_end { Some(len) } else { None },
        });
    }
    bursts
}

/// Number of edits in the log made through the agent.
#[must_use]
pub fn agent_edit_count(history: &[CommandRecord]) -> usize {
    history
        .iter()
        .filter(|record| record.origin == Origin::Agent)
        .count()
}
