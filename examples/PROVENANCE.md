# Example provenance

The examples are retained as scientific inputs, not as anonymous test blobs.  Each
file should remain traceable to a calculation, publication, or documented tutorial.
Where a file header contains a citation, that header is authoritative and must not be
removed during formatting or fixture reduction.

## Elasticity

- `elasticity/calcite.dat`: calcite B3LYP elastic tensor used to demonstrate
  second-order elasticity analysis.
- `elasticity/hydroxylapatite.dat`: hydroxylapatite elastic tensor and density used
  throughout the elasticity and seismic validation suites.
- `seismic/tutorial_api.py`: executable public-API tutorial that reuses the
  hydroxylapatite tensor without duplicating its scientific data.

## Equations of state

- `eos/PV_quartz.dat`: experimental quartz P-V dataset distributed as an EOS fitting
  example.
- `eos/PV_topaz.dat`: topaz volume and axial compression data from Gatta et al.
  (2014), with the complete citation retained in the file.
- `eos/VT_rutile.dat`: rutile thermal-expansion data derived from the EosFit7 example
  distribution; the source and literature selection are recorded in the header.
- `eos/T_triclinic.dat`: triclinic plagioclase thermal data citing Tribaudino et al.
  (2011).
- `eos/PVT_NaF.dat`: NaF P-V-T data from Clough et al. (2025), including molar-volume
  units.  This file validates the historical EOS convention in which `VSCALE` can
  carry an absolute unit.
- `eos/specs/*.spec`: tutorial batch plans authored for Quantas 2.0 from the
  curated datasets above.  They add no observational data and record only model,
  solver, acceptance, and presentation choices.
- `eos/scripts/*.py`: executable public-API examples derived from the same
  datasets and batch plans; they are documentation assets rather than independent
  scientific sources.

## HA and QHA

- `qha/crystal-qha/`: MgO/periclase B3LYP CRYSTAL QHA output and the normalized
  Quantas YAML generated from it.
- `qha/crystal-phonons/`: dolomite PBE0 CRYSTAL phonon calculations over seven
  volumes, their file list, and the normalized Quantas YAML.

The large CRYSTAL outputs are used only by the isolated examples/scientific-regression
suite.  Unit tests should use reduced fixtures unless the complete output is required
to validate parser behavior.

## Thermoelasticity

- `thermoelasticity/crystal_outputs/`: real dolomite PBE0 and MgO B3LYP CRYSTAL
  second-order elastic calculations over volume/pressure series.
- `thermoelasticity/dol_pbe0_thermoelastic.yaml`: normalized quasi-static
  thermoelastic input generated from the dolomite series.
- `thermoelasticity/prem_custom_temperature.*`: an illustrative PREM-pressure profile
  combined with a user-defined temperature table.  The temperature table is not a
  geological reference model and must be replaced for quantitative interpretation.
- `thermoelasticity/dolomite_continental_profile.yaml`: illustrative composed
  profile used by the thermoelasticity tutorial. Pressure is evaluated from PREM;
  temperature is generated from an explicitly parameterized layered conductive
  continental model. It is a reproducible teaching path, not a universal regional
  geotherm.
- `thermoelasticity/tutorial_api.py`: executable public-API reproduction of the
  complete dolomite workflow. It combines the curated QHA and static SOEC inputs and
  generates only derived results outside the source example tree.


## Interoperability

- `interoperability/workflow_cli.sh` and `interoperability/workflow_api.py`:
  executable documentation examples that combine the existing curated dolomite
  QHA and thermoelastic datasets. They introduce no new observational or
  first-principles data; all generated HDF5, report, and state files are derived
  products written outside the source example tree.

## Input/output utilities

- `io/inspect_hdf5.py`: documentation example that traverses an existing native
  Quantas HDF5 file and reports groups, datasets, attributes, units, dtypes,
  shapes, and storage layout without loading complete large arrays.
- `io/extract_hdf5.py`: documentation example that copies selected numerical
  datasets to a compressed NumPy archive and writes a JSON manifest preserving
  source paths and relevant metadata. It creates no independent scientific data.

## Maintenance rules

1. Preserve original scientific numbers and header citations.
2. Normalize distributed text to UTF-8/LF without changing numeric content.
3. Regenerate `manifest.json` and `MANIFEST.sha256` after an intentional change.
4. Add or update a scientific-regression assertion when replacing a real dataset.
5. Do not duplicate large files under `tests/data`; tests should reference this
   directory through the repository root.
