# EOS examples

The EOS examples cover volumetric and axial P--V fitting, V--T thermal
expansion, and coupled P--V--T fitting.  File headers retain source citations,
column definitions, and unit conventions.

## Tutorial material

The `specs/` directory contains three strict `QUANTAS EOS SPEC 1` batch plans:

- `quartz_pv_tutorial.spec`: solver and P--V model comparison;
- `rutile_vt_tutorial.spec`: solver, V--T model, and axial-expansion comparison;
- `naf_pvt_tutorial.spec`: comparison of four P--V--T couplings, including MGD.

The `scripts/` directory contains public-API equivalents for running the batch
plans and representative P--V, V--T, and P--V--T fits.  These scripts import
only `quantas.api` and are used by the documentation as executable examples.

Run a batch from the command line, for example, with:

```console
quantas eos run PV_quartz.dat --spec specs/quartz_pv_tutorial.spec \
    --output quartz_pv.hdf5 --report quartz_pv.log --force
```

See the EOS tutorial section of the manual for interpretation, diagnostics,
post-fit calculation, plotting, and exercises.
