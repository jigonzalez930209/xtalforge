<!--
Pull request template. The checks mirror AGENTS.md sections 2 and 3.
-->

## What changed

<!-- One paragraph: the capability or fix, and the roadmap subtask. -->

- Subtasks:
- Scope:

## Checklist

- [ ] Roadmap checkboxes updated in `docs/roadmap/`
- [ ] Parity matrix updated in `docs/api/parity.md`
- [ ] API page added or updated for every public operation
- [ ] Algorithm page added or updated when behaviour changed
- [ ] `CHANGELOG.md` `Unreleased` entry added
- [ ] `scripts/check-rules.sh` passes (500 lines, no warning suppressions)
- [ ] `cargo fmt --check` passes
- [ ] `cargo clippy --workspace --all-targets --all-features -- -D warnings` passes
- [ ] `cargo test --workspace --all-features` passes
- [ ] Python and WASM tests pass when the surface is affected
- [ ] `scripts/check-parity.sh` passes
- [ ] Commits are directory-scoped with `feat(<scope>):` or `ref(<scope>):`
      and the completed tasks are listed from hardest to easiest

## Notes for the reviewer

<!-- Deviations, tolerances, known gaps. -->
