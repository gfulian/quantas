# Quantas release procedure

A release is permitted only from a clean tagged commit.  Numerical reliability and
reproducibility take precedence over publication speed.

## 1. Prepare

- Confirm the version in `src/quantas/_version.py`, `CITATION.cff`, and the changelog.
- Review the frozen `quantas.api` namespace snapshots and registry capabilities; every public-symbol change must be deliberate.
- Close or explicitly defer every release-blocking item in `ROADMAP.md`.
- Regenerate and verify the examples manifest.
- Confirm that public source archives contain no internal phase reports, transient
  validation products, credentials, local paths, or unpublished datasets.

## 2. Validate

From a fresh source checkout or unpacked archive, run:

```bash
bash tools/validate_release.sh
```

The script recreates ``.venv-validation``, installs ``.[dev]``, stops at the
first failed command, records the complete environment and output under
``validation/``, and executes the manifest, Ruff, mypy, compile, staged pytest,
architecture, Sphinx, build, Twine, distribution, and repository checks.

Also inspect the scientific validation matrix and compare all release reference
results against the approved baseline.

## 3. Test distribution

Create a GitHub Actions manual release run targeting the protected `testpypi`
environment.  Install the uploaded version from TestPyPI in a clean environment and
repeat the CLI/import smoke tests.  TestPyPI and PyPI are separate repositories and
must have separate Trusted Publisher configurations.

## 4. Publish

- Create a signed `vX.Y.Z` tag on the validated commit.
- Publish a GitHub release from that tag.
- Allow the protected `pypi` environment to publish through OpenID Connect Trusted
  Publishing; no long-lived PyPI API token belongs in repository secrets.
- Verify metadata, files, installation, `quantas --version`, and the rendered project
  page after publication.

## 5. Record

- Archive the exact source commit, distribution hashes, CI run, environment matrix,
  and validation summary in the separate internal development record.
- Begin the next changelog section immediately.
