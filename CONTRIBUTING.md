# Contributing to Quantas

Contributions to Quantas are welcome, particularly when they improve the
scientific reliability of the package, extend its interoperability, add useful
validation data, or make the documentation easier to use.

Quantas is first a scientific software: changes to numerical methods, physical
conventions, units, or data interpretation therefore require a little more care
than ordinary code maintenance. The aim is not simply to keep the test suite
green, but to make sure that a change remains scientifically traceable and that
the same calculation gives the same result regardless of whether it is reached
through the Python API or the command-line interface.

## Preparing a development environment

Quantas currently supports Python 3.10 through 3.13.

A normal development checkout can be installed in an isolated environment with:

```bash
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
```

The editable installation is convenient during development because changes to
the source tree are immediately available to the active environment.

For substantial changes, please run the complete test suite before considering
the work finished:

```bash
python tools/run_tests.py all -- -q
```

The repository also uses Ruff, mypy, architecture checks, documentation builds,
and distribution checks as part of its validation process. The exact release
procedure is documented separately because release validation is intentionally
more demanding than routine development.

## Working on Quantas

New development normally starts from the `development` branch. Changes that are
large enough to benefit from isolation should be developed on short-lived
branches such as `feature/...` or `fix/...` and merged back into `development`
when they are ready.

The `master` branch represents the current tested public baseline and should not
be used as the ordinary integration branch for unfinished work.

The historical Quantas 0.9.1 code is retained separately for reference. It is
useful when checking legacy behaviour or historical naming, but it is not the
architectural model for Quantas 2.

The detailed branch and worktree conventions used by the project are described
in [DEVELOPMENT.md](DEVELOPMENT.md).

## Architecture

Quantas is designed so that scientific calculations do not depend on a
particular frontend.

`quantas.core` contains reusable numerical, physical, chemical, geometrical,
and event-related building blocks. It must remain independent of Click, Rich,
Dash, Matplotlib, Plotly, terminal streams, progress bars, and GUI objects.

`quantas.models` contains shared calculators, readers, exporters, and passive
data contracts. Active or stateful objects use normal classes, while
dataclasses are reserved for data that are genuinely passive.

Code-specific parsing belongs under `quantas.interfaces`. Scientific workflows
such as Elasticity, SEISMIC, HA, QHA, Thermoelasticity, and EOS live under
quantas.modules.

Rendering is kept separate from calculation. The command-line interface may use
Click and Rich, but scientific code must not depend on either of them. The
public Python API and the CLI ultimately call the same scientific
implementations.

## Scientific changes

Changes to formulas or numerical methods should begin by understanding and
capturing the behaviour that is already present.
Where possible, a numerical change should be supported by characterization
tests and compared with an analytical result, the historical Quantas
implementation, a published reference, or an independent implementation.

Formulas, physical conventions, units, array shapes, and accepted tolerances
should not change incidentally as part of a structural cleanup.
Real scientific calculations and native HDF5 values use double precision
(`float64`), while complex quantities use `complex128`. Display precision is a
presentation choice and must not alter stored scientific data.

A numerical modification should normally be kept separate from a large
refactoring. This makes both the scientific reasoning and the review history
much easier to follow.
When a scientific behaviour changes intentionally, the corresponding tests,
documentation, provenance, and public-facing interfaces should be updated
together.

## Public API and data formats

Quantas 2 is still in beta. Its public API is deliberately controlled by
contract tests, but limited changes are still possible before the first stable
2.0 release.

Changes to `quantas.api`, native HDF5 schemas, structured input formats, or
other application-facing contracts should therefore be deliberate and
documented. Compatibility requirements will become stricter once the first
stable Quantas 2 release is published.

## Documentation and validation

New public functionality should include enough documentation for another user
to understand what the feature does, how to use it, and which scientific
assumptions it makes.

Tests should be chosen for the behaviour being changed rather than simply for
code coverage. Analytical limits, regression data, invalid-input cases,
persistence round trips, and comparisons with independent references are all
valuable when they are scientifically relevant.

Bug reports and validation datasets are particularly useful when they include
enough information to reproduce the calculation.

## Pull requests

A pull request should describe the problem being addressed, the scientific or
technical reasoning behind the chosen solution, and the validation that has
been performed.

Small, focused changes are easier to review and safer to integrate than a
mixture of numerical, architectural, formatting, and documentation changes in
one large patch.

If a change affects scientific results, please make that explicit.