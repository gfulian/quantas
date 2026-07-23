Reporting, warnings, and events
================================

Quantas distinguishes between scientific results and frontend presentation.

Reports and plots
-----------------

Workflows can build frontend-neutral report tables and plot specifications.
These contracts can then be rendered by supported frontends without changing
any numerical result.

Events and progress
-------------------

Long-running workflows emit structured events.  These events can represent:

- informational messages;
- progress updates;
- warnings;
- structured result notifications.

Operational progress is meant for observers and frontends.  It is not a part of
scientific result persistence.

Warnings
--------

Scientifically meaningful warnings, extrapolation notices, and other important
workflow diagnostics are stored with results when appropriate.
