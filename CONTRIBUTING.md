# Contributing to Quantas

Quantas is a scientific library first. Numerical reliability and stable public
contracts take priority over cleanup, optimization, or frontend convenience.

## Development setup

Create an isolated Python 3.10 or newer environment and install the complete
development dependency set:

```bash
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
```

Run the deterministic validation matrix with:

```bash
python tools/run_tests.py all -- -q
ruff check src tests tools
mypy
python -m compileall -q src tests tools
```

Build and validate distributable archives with:

```bash
python -m build
python -m twine check dist/*
python tools/check_distribution.py dist
```

## Architectural rules

- `quantas.core` contains reusable numerical physics and must not import Click,
  Rich, Matplotlib, Dash, or terminal objects.
- `quantas.models` contains passive contracts and shared active base classes.
- `quantas.interfaces` contains code-specific parsers.
- `quantas.modules` contains complete scientific workflows.
- `quantas.renderers` owns frontend-specific presentation.
- `quantas.cli` parses options and delegates; it does not own scientific logic.
- A separate Quantas GUI package should depend on the organized
  `quantas.api` namespaces, registry capabilities, and neutral
  result/report/plot contracts, never on CLI or internal module facades.

## Scientific changes

Before changing a formula or numerical method:

1. capture current behaviour with characterization tests;
2. preserve units, shapes, conventions, and tolerances;
3. compare against analytical, legacy, published, or external references;
4. keep calculations in `float64` or `complex128`;
5. update theory, API, CLI, persistence, and validation documentation together.

Do not combine a numerical change with a large structural refactor.

## Public contracts

Quantas is still in pre-release refactoring. The organized `quantas.api` surface is
frozen by contract tests, but no public deprecation cycle is active yet. Introduce a formal
deprecation policy only when public releases begin. Until then, changes to
application-facing facades must be deliberate, documented, and covered by API
contract tests.
