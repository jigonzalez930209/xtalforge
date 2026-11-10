# Phase 0 — Remaining work

> **Update 2026-11-10.** Actions were enabled on 2026-11-08, only the tag-driven
> `release.yml` is kept, and the workflow has published `v0.0.1`, `v0.0.2`,
> `v0.0.3`, `v1.0.0-rc.1` and `v1.0.0` (crates.io, PyPI, npm, GitHub release and
> Pages).
> The names are published. Tasks 0.8.1-0.8.5 (CI matrix) stay open by
> maintainer decision because the runner cost was not wanted; the checks run
> locally (`AGENTS.md` section 5.4). The sections below are the original
> 2026-11-08 plan, kept for the record.

Status on 2026-11-08: every Phase 0 subtask is implemented locally except the
items below. The full local check list passes; the open items are external
settings or approvals, not missing code.

## 0.6.3 Release workflow

`.github/workflows/release.yml` is written: tag trigger, version/changelog
checks, Rust/Python/WASM/docs builds, crates.io and PyPI publishing, npm
provenance, GitHub release and Pages deployment. It has **not run** because
GitHub Actions are disabled at repository level by maintainer instruction.

To close it:

1. Enable Actions: `gh api --method PUT repos/jigonzalez930209/xtalforge/actions/permissions -F enabled=true`.
2. Run `workflow_dispatch` with `dry_run: true`.
3. Fix any job that fails; repeat until the dry run is green.
4. Before `v0.0.1`, run it with publishing disabled one final time.

## 0.6.5 Name reservation

`xtalforge` and `xtalforge-wasm` were free on 2026-11-08 (crates.io search
empty, PyPI 404, npm 404). No placeholder was published; the roadmap requires
the maintainer's explicit agreement before publishing a `0.0.0` placeholder.
`xtalforge-core`, `xtalforge-io` and `xtalforge-symmetry` are also free and will
be reserved with the same release.

To close it: approve the placeholders, add the registry tokens, publish
`0.0.0` placeholders, and record the dates in `AGENTS.md` section 6.

## 0.7.5 Documentation deployment

The docs deployment is the `docs` job of `.github/workflows/release.yml`. It
runs on a version tag, waits for the verification and WASM jobs, downloads the
already-built web package (no second Rust build), runs VitePress and deploys to
GitHub Pages. It has not run because Actions are disabled. A manual
`Docs` workflow remains for recovery deployments; it does not run on `main`.
Versioned documentation per MINOR release starts at `v0.1.0` and is not needed
for Phase 0.

To close it: enable Actions, run the workflow on `main` and on a test tag, and
verify the deployed URL from the repository Pages settings.

## 0.8.1–0.8.5 Continuous integration

The CI, rules and manual docs workflows were removed by maintainer decision on
2026-11-08 to avoid runner cost. The local check list in `AGENTS.md` section
5.4 is the gate. The only workflow left is the tag-driven `release.yml`. **No
CI job runs.**

To close it:

1. Enable Actions.
2. Mark `Rules / Repository rules`, `Format and lints`, the Rust matrix, Python,
   WASM and Docs jobs as required status checks in branch protection.
3. Run the matrix on a pull request and fix every failure.
4. Record the measured runtime; the target is under 15 minutes.

## Deferred surface details (tracked, not blocking Phase 0 exit)

- **0.5.4 Python interop.** ASE `Atoms` and pymatgen `Structure` helpers need
  the public `Structure` type; they start in Phase 1 Task 1.17.2.
- **0.5.5 threaded WASM.** The opt-in `wasm-bindgen-rayon` build is documented
  in `docs/api/bindings.md`; it is enabled only when a threaded use case
  exists, because it requires cross-origin isolation.

## Closing checklist

- [x] Actions enabled by the maintainer (2026-11-08).
- [x] Release workflow green (superseded by four real tag releases).
- [ ] CI matrix green on Linux, macOS and Windows (declined: no CI workflow).
- [x] Docs deployment green (GitHub Pages, every tag).
- [x] Name placeholders approved and published (2026-11-09).
- [x] `AGENTS.md` section 6 updated with the reservation dates.
