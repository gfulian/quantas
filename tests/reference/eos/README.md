# EOS analytical reference data

These JSON files contain frozen, implementation-independent values used by the
EOS scientific regression tests. They were generated from explicit analytical
formulae outside the production implementation and retained under stable,
scientific names rather than refactoring-phase names.

- `eos_vt_formula_reference.json`: volume-temperature model values and thermal
  expansion coefficients.
- `eos_pvt_formula_reference.json`: coupled pressure-volume-temperature values
  and coupling properties.

Generation scripts and internal development records are kept in the separate
maintainer development archive. Changes to these values require a documented
scientific revalidation, not a routine test update.
