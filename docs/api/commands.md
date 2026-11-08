# Commands and document

`xtalforge-commands` is the only mutation path. A `Command` returns the new
`Structure` from `apply`; `CommandStack` owns the undo/redo history, the clean
flag and a bounded log. `Document` owns the structure, the stack, the selection,
the atom groups and the origin of each edit.

```rust
use xtalforge::commands::{CommandStack, Origin, atoms::MoveSites};

let mut stack = CommandStack::new();
let edited = stack.push(Box::new(MoveSites::new(vec![0], [0.1, 0.0, 0.0])), &structure, Origin::User)?;
```

## Stack

`CommandStack::push(command, &structure, origin)` applies and records;
`undo(&mut structure)` and `redo(&mut structure)` move through the history.
`can_undo`, `can_redo`, `undo_depth`, `redo_depth`, `is_dirty`, `mark_clean`,
`max_depth`/`set_max_depth` (default 200) and `log()` expose the state.
`MacroCommand` groups commands; `ReplaceStructure`, `SnapshotEdit` and
`StructureOperation` cover whole-structure and function-based edits.

## Document

`Document::new(structure)` starts clean. `apply(command, origin)`, `undo`,
`redo`, `dirty`, `mark_clean`, `structure()`, `selection()` and the atom-group
methods are the public state. `Origin` is `User`, `Script` or `Agent`.

## Catalogue

- `atoms`: `AddSites`, `DeleteSites`, `SetElement`, `MoveSites`, `SetOccupancy`.
- `bonds`: `AddBond`, `RemoveBond`, `SetBondRules`.
- `cell`: `SetLattice` (`KeepFractional`/`KeepCartesian`), `Supercell`,
  `ShiftOrigin`, `WrapIntoCell`.
- `symmetry`: `ReduceToP1`, `FindSymmetry`, `SetSpaceGroup`, `Standardize`.
- `groups`: `AtomGroup`, `AddAtomGroup`, `RemoveAtomGroup`, `SetGroupColor`,
  `SetGroupHidden`.
- `clipboard`: `Fragment`, `copy_selection`, `PasteFragment`.

Every command is tested for apply and undo; invalid indices and values return
`Error::InvalidInput`.
