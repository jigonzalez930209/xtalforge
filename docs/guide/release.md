# Release and deployment

The release is driven by `.github/workflows/release.yml`. A tag `vX.Y.Z` builds the CLI for Linux, Windows and macOS, publishes the
GitHub release, then publishes PyPI, npm and crates.io from `ubuntu-24.04`,
then builds and deploys the documentation. A manual `workflow_dispatch` only
runs the CLI builds and the version check. The separate `Docs` workflow is
manual-only and exists as a recovery path; it does not run on pushes to
`main`.

## 1. Enable Actions and Pages

1. **Settings → Actions → General**: allow GitHub Actions.
2. **Settings → Pages → Build and deployment**: source **GitHub Actions**,
   not “Deploy from a branch”. The site is VitePress; the legacy Jekyll build
   of `main` is not the site.

## 2. Secrets

Add them in **Settings → Secrets and variables → Actions → New repository
secret**. The workflow reads the names below.

| Secret | Registry | How to create it | Required |
|---|---|---|---|
| `CARGO_REGISTRY_TOKEN` | crates.io | crates.io → Account Settings → API Tokens → New Token; scopes `publish-new`, `publish-update`, `yank` | yes |
| `NPM_TOKEN` | npm | npm → Access Tokens → Generate New Token; granular token with read/write on `xtalforge-wasm` | yes |
| `PYPI_API_TOKEN` | PyPI | PyPI → Account settings → API tokens; scope the `xtalforge` project (use an account token for the first upload) | yes |
| `GITHUB_TOKEN` | GitHub | Provided automatically by Actions | automatic |

The workflow uploads the wheel and sdist with `twine`, using
`PYPI_API_TOKEN`. If you prefer PyPI Trusted Publishing, replace the twine
step with `pypa/gh-action-pypi-publish` and remove the token secret.

### Optional signing secrets

The current installers are unsigned. If code signing is added, these are the
conventional names:

| Secret | Use |
|---|---|
| `APPLE_CERTIFICATE` | base64 Developer ID certificate |
| `APPLE_CERTIFICATE_PASSWORD` | certificate password |
| `APPLE_SIGNING_IDENTITY` | Developer ID Application identity |
| `APPLE_ID` / `APPLE_TEAM_ID` / `APPLE_APP_SPECIFIC_PASSWORD` | macOS notarization |
| `WINDOWS_CERTIFICATE` / `WINDOWS_CERTIFICATE_PASSWORD` | Windows Authenticode |

They are not consumed by the workflow yet; unsigned macOS and Windows builds
will trigger the usual Gatekeeper/SmartScreen warnings.

## 3. What the workflow publishes

- **crates.io**: `xtalforge-core`, `xtalforge-io`, `xtalforge-symmetry`,
  `xtalforge-view`, `xtalforge-commands`, `xtalforge-ff`,
  `xtalforge-analysis`, `xtalforge-build`, `xtalforge-mof`,
  `xtalforge-modules`, the `xtalforge` facade and `xtalforge-cli`, in dependency
  order. Existing versions are skipped so a rerun is idempotent.
- **PyPI**: an abi3 manylinux wheel and an sdist built on `ubuntu-24.04`,
  uploaded with `twine`.
- **npm**: `xtalforge-wasm`, built for the bundler target on `ubuntu-24.04`
  with provenance.
- **GitHub release**:
  - Linux x86_64: `.deb`, `.rpm`, `.AppImage`, `.tar.gz`
  - Windows x86_64: `.zip`
  - macOS universal (x86_64 + arm64): `.tar.gz`, `.dmg`
  - `SHA256SUMS.txt`
- **GitHub Pages**: the VitePress site, built on `ubuntu-24.04` after the
  package publications, using the WASM web build from the publish job.

## 4. Caching and runtime

Every Rust job uses `Swatinem/rust-cache` with a job-specific `shared-key`, so
the next release reuses the dependency and target caches instead of compiling
from zero. The CLI jobs for Linux, Windows and macOS run in parallel; the macOS
job builds a universal binary. Package publication happens in one
`ubuntu-24.04` job (PyPI, npm, crates.io), and the docs job reuses the
`wasm-web` artifact instead of rebuilding Rust. No test suite runs in the
release workflow.

## 5. Release checklist

1. Set the version in `Cargo.toml` and run `scripts/sync-version.sh`.
2. Move the `Unreleased` changelog entries under `## [X.Y.Z] - YYYY-MM-DD`.
3. Run the full local CI gate from `AGENTS.md` section 5.4 on the exact commit
   to tag. There is no GitHub CI; this gate is mandatory.
4. Commit, tag and push:

   ```bash
   git tag -a vX.Y.Z -m "xtalforge X.Y.Z"
   git push origin vX.Y.Z
   ```

5. Watch the `Release` workflow; it does not run tests. The GitHub release
   and Pages deployment are the last jobs.
6. If a publish step fails after some registries accepted the version, fix the
   cause and re-run the workflow. crates.io and PyPI versions cannot be
   overwritten; the workflow skips crates whose exact version already exists.
7. The first crates.io release creates many new crates. crates.io rate-limits
   new crate creation; the workflow retries with backoff and may take a while.
