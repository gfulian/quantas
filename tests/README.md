# Quantas test suites

Tests are organized by scientific domain and classified on a second axis by execution
concern.  Unit tests use minimal synthetic data when this isolates one formula or
failure mode; `tests/examples/` uses the curated real inputs under `examples/` for
reader, workflow, persistence, and scientific-regression coverage.

## Direct pytest usage

```bash
pytest --suite core -v
pytest --suite qha -v
pytest --suite examples -v
pytest --suite scientific-regression -v
pytest --suite qha-plotting -v
pytest --suite cli -v
pytest --suite all -v
```

The short module names select Python-library/workflow tests only.  Add `-all`, `-cli`,
or `-plotting` where that module exposes the corresponding suite.  The `examples`
suite is intentionally isolated because it reads large external-code outputs and real
scientific datasets.

List every preset with:

```bash
pytest --list-suites
```

## Deterministic complete run

```bash
python tools/run_tests.py all -- -q
```

The runner executes fifteen isolated processes:

```text
core -> elasticity -> seismic -> HA -> QHA -> EOS -> thermoelasticity
     -> CLI -> six domain plotting suites -> real examples
```

Isolation prevents numerical-thread, terminal, parser, and graphical renderer state
from leaking between stages.  Automatic loading of unrelated pytest plugins is
disabled; each process receives an isolated Matplotlib configuration and a per-stage
timeout.  Use `--stage-timeout 0` to disable the timeout or `--allow-plugin-autoload`
when an external plugin is deliberately required.
