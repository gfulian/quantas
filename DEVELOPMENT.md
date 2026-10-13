# Quantas development workflow

Quantas uses a small number of long-lived branches and short-lived branches for
individual pieces of work. The intention is to keep the repository easy to
understand while still allowing scientific features to be developed and
validated independently.

## Long-lived branches

The `master` branch is the current public baseline. It should always represent
a working and tested version of Quantas. Repository visitors who clone Quantas
without selecting another branch receive this version.

Ongoing development is integrated through `development`. New scientific
features, interfaces, fixes, and larger documentation changes normally begin
from this branch and return to it after review and validation.

The historical Quantas 0.9.1 line is retained as `legacy/0.9.1`. It is frozen
except for exceptional historical or archival reasons and is not used as a
starting point for Quantas 2 development.

During the transition to this branch model, `dev/refactor` is retained
temporarily because existing teaching and laboratory material still refers to
it. New development should use `development` instead. Once those external
dependencies have been updated, `dev/refactor` can be retired.

## Feature and fix branches

Work that should not be developed directly on `development` can use a
short-lived branch.

Typical names are:

```text
feature/eosfit-pvt
feature/full-qha
feature/phonopy-interface

fix/qha-volume-units
fix/eos-report-format
```

These names are descriptive rather than prescriptive. The important point is
that the branch has one clear purpose and can be removed once its work has been
integrated.

A new feature normally starts from an up-to-date development branch:

```bash
git switch development
git pull --ff-only
git switch -c feature/example
```

When worktrees are being used, the same branch can instead be created directly
as a new worktree.

## Local worktrees

The main Quantas development workspace uses Git worktrees so that several
branches can remain checked out at the same time.

A typical local layout is:

```text
worktrees/
  - quantas-master/
  - quantas-development/
  - quantas-dev-refactor/
  - quantas-legacy-0.9.1/
  - features/
  - fixes/

Feature worktrees can then be created below features/, for example:

```bash
git -C ./worktrees/quantas-development worktree add `
    "/path/to/quantas/worktrees/features/quantas-new-feature" `
    -b feature/new-feature development
```

The exact local directory structure is not part of the public Quantas
repository. It is simply a convenient way to keep active branches isolated
from one another.

## Integrating development

A feature should return to development only after its scientific and
technical behaviour has been checked.

For small changes this may require only focused tests. Larger changes should
also exercise the relevant module, API, CLI, persistence, documentation, and
scientific regression coverage.

The complete repository test runner is:

```bash
python tools/run_tests.py all -- -q
```

Static and architectural checks used during development include:

```bash
python -m ruff check src tests tools docs/tools
python -m mypy src/quantas
python tools/check_architecture.py
```

Documentation changes should also be checked with the project's Sphinx build,
particularly when warnings are treated as errors.

## Promoting development to master

`development` is promoted to `master` when it represents a coherent,
well-tested public checkpoint.

This is deliberately different from merging every individual feature directly
into master. The integration branch may contain several related changes that
need to be validated together before they become the public baseline.

Before promotion, the complete validation gate should be run from a clean
checkout or equivalent clean source tree. Changes to versions, release notes,
and tags are handled as part of the release process rather than during ordinary
feature development.

## Tags and releases

Branches describe where development is happening. Tags identify immutable
versions of the software.

Release tags use the package version with a leading v, for example:

```text
v0.9.1
v2.0.0b13
v2.0.0rc1
v2.0.0
v2.0.1
```

A tag should refer only to a commit that has passed the appropriate release
validation. Detailed maintainer instructions for creating and publishing a
release are kept in the development documentation rather than in this file.

## Hot fixes

If a serious problem is discovered in the current public `master` while
development already contains unrelated future work, the correction can be
made on a short-lived `hotfix/...` branch created from `master`.

After validation, the fix is merged into both `master` and `development` so
that the correction is not lost from future versions.

Routine fixes that do not require this separation should simply follow the
normal development workflow.

## Keeping the repository clean

Finished feature and fix branches should normally be removed after their work
has been integrated. Their commits remain part of Git history, so keeping an
obsolete branch indefinitely does not preserve additional scientific
information.

Likewise, worktrees are working directories rather than archives. Once the
corresponding branch has been merged and the worktree is clean, the worktree
can be removed.

The aim is for the long-lived remote branch list to remain small:

```text
master
development
legacy/0.9.1
```

with only genuinely active feature, fix, or temporary compatibility branches
alongside them.