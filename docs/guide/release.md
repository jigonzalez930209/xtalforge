# Release and deployment

The release is driven by `.github/workflows/release.yml`. A tag `vX.Y.Z` runs
the verification, builds and publishes the crates, wheels, npm package, CLI
packages and documentation, then creates the GitHub release. A manual
`workflow_dispatch` runs the same builds without publishing. The separate
`Docs` workflow is manual-only and exists as a recovery path; it does not run
on pushes to `main`.

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
| `PYPI_API_TOKEN` | PyPI | PyPI → Account settings → API tokens; scope the `xtalforge` project after the first upload | no |
| `GITHUB_TOKEN` | GitHub | Provided automatically by Actions | automatic |

PyPI uses **Trusted Publishing** in the workflow, so no PyPI secret is needed:
create a pending publisher in PyPI → Publishing with owner
`jigonzalez930209`, repository `xtalforge`, workflow `release.yml` and (leave
the environment empty). If the project cannot use Trusted Publishing yet, add
`PYPI_API_TOKEN` and pass it as `password` to
`pypa/gh-action-pypi-publish`.

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
- **PyPI**: abi3 wheels for Linux, macOS x86_64/arm64 and Windows, uploaded
  with attestations. The sdist is not built at this stage.
- **npm**: `xtalforge-wasm`, built for the bundler target with provenance.
- **GitHub release**:
  - Linux x86_64: `.deb`, `.rpm`, `.AppImage`, `.tar.gz`
  - Windows x86_64: `.zip`
  - macOS x86_64 and arm64: `.tar.gz`, `.dmg`
  - `SHA256SUMS.txt`
- **GitHub Pages**: the VitePress site, deployed only on a version tag after
  the verification and WASM jobs pass.

## 4. Caching and runtime

Every Rust job uses `Swatinem/rust-cache` with a job-specific `shared-key`, so
the next release reuses the dependency and target caches instead of compiling
from zero. The docs job does **not** rebuild Rust: it downloads the `wasm-web`
artifact produced by the `wasm` job, runs `npm ci` and builds VitePress. The
WASM job and the verification run in parallel; wheels and the three CLI
platforms also run in parallel after verification.

## 5. Release checklist

1. Set the version in `Cargo.toml` and run `scripts/sync-version.sh`.
2. Move the `Unreleased` changelog entries under `## [X.Y.Z] - YYYY-MM-DD`.
3. Run the local release checks from `AGENTS.md` section 5.4.
4. Commit, tag and push:

   ```bash
   git tag -a vX.Y.Z -m "xtalforge X.Y.Z"
   git push origin vX.Y.Z
   ```

5. Watch the `Release` workflow; the GitHub release and Pages deployment are
   the last jobs.
6. If a publish step fails after some registries accepted the version, fix the
   cause and re-run the workflow. crates.io and PyPI versions cannot be
   overwritten; the workflow skips crates whose exact version already exists.
7. The first crates.io release creates many new crates. crates.io rate-limits
   new crate creation; the workflow retries with backoff and may take a while.
