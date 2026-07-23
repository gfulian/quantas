HDF5 and provenance
===================

HDF5 is the native result format of Quantas.

A result archive typically stores:

- program and version;
- module and method identifiers;
- schema version;
- creation metadata;
- normalized input;
- workflow options;
- results and uncertainties;
- warnings and scientifically relevant events;
- numerical precision metadata;
- citation keys or related provenance information.

Readers identify Quantas HDF5 content from metadata rather than from the file
name alone.

Workflow-specific HDF5 payloads are described in :doc:`../../formats/index`.
