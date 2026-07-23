Inputs, options, and results
============================

The public API of Quantas uses a common vocabulary whenever it is scientifically
appropriate.

- **Input** objects carry normalized scientific data.
- **Options** objects carry numerical and workflow configuration.
- **Result** objects collect metadata, workflow outputs, warnings, and
  references to generated reports or plots.

The canonical public entry points are exposed from :mod:`quantas.api` and often
follow the pattern below:

- ``read_input()``
- ``normalize_input()``
- ``run()``
- ``get_result()``
- ``read_result()``
- ``write_result()``
- ``build_report()``
- ``build_plots()``

Some workflows, especially EOS, legitimately expose a more specialized
vocabulary because the scientific problem is organized around datasets,
fit requests, diagnostics, and archives.

For the shared passive contracts used across modules, see :mod:`quantas.api.common`.
