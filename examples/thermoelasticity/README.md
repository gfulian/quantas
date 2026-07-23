# Thermoelasticity examples

This directory contains curated inputs and executable examples for the Quantas
quasi-static thermoelastic workflow.

## Dolomite tutorial

- `dol_pbe0_thermoelastic.yaml`: normalized static SOEC-versus-volume input.
- `dolomite_continental_profile.yaml`: illustrative layered continental
  pressure-temperature profile used by the public tutorial.
- `tutorial_api.py`: end-to-end public-API workflow covering QHA preparation,
  cold finite-strain calibration, Point, Grid, and Profile analyses, tabular
  export, interoperability, and selected figures.
- `crystal_outputs/`: original CRYSTAL calculations retained for provenance and
  input-generation validation.

Run the complete API tutorial from the repository root with:

```bash
python examples/thermoelasticity/tutorial_api.py --output-dir thermoelastic_tutorial
```

The profile specification is pedagogical rather than a universal geological
model.  Its pressure branch follows PREM while its temperature branch is a
layered conductive continental construction with explicitly stated parameters.
Replace those parameters for quantitative interpretation of a specific region.

## Additional examples

- `prem_custom_temperature.*`: PREM pressure combined with a user-supplied
  temperature table.
- `scripts/`: focused API examples for input generation, fitting, depth
  profiles, and plotting foundations.
