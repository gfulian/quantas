#!/usr/bin/env bash
set -euo pipefail

output_dir="${1:-interoperability_cli}"
mkdir -p "$output_dir"

quantas qha run examples/qha/crystal-phonons/dol_pbe0.yaml \
  --scheme freq \
  --minimization poly \
  --thermal-expansion mixed_derivative \
  --temperature 300 1800 100 \
  --pressure 0 10 2 \
  --energy-degree 3 \
  --frequency-degree 3 \
  --output "$output_dir/dolomite_qha.hdf5" \
  --report "$output_dir/dolomite_qha.log" \
  --no-progress \
  --force

quantas thermoelasticity run \
  examples/thermoelasticity/dol_pbe0_thermoelastic.yaml \
  "$output_dir/dolomite_qha.hdf5" \
  --reference-eos BM3 \
  --finite-strain-order 3 \
  --adiabatic auto \
  --validation standard \
  --output "$output_dir/dolomite_fit.hdf5" \
  --report "$output_dir/dolomite_fit.log" \
  --no-progress \
  --force

quantas thermoelasticity analysis point \
  "$output_dir/dolomite_fit.hdf5" 5 800 \
  --tensor-condition adiabatic \
  --extrapolation fail \
  --output "$output_dir/dolomite_5GPa_800K.dat" \
  --force

quantas elasticity run "$output_dir/dolomite_5GPa_800K.dat" \
  --output "$output_dir/dolomite_state_elasticity.hdf5" \
  --report "$output_dir/dolomite_state_elasticity.log" \
  --no-progress \
  --force

quantas seismic run "$output_dir/dolomite_5GPa_800K.dat" \
  --level phase \
  --ntheta 31 \
  --nphi 61 \
  --output "$output_dir/dolomite_state_seismic.hdf5" \
  --report "$output_dir/dolomite_state_seismic.log" \
  --no-progress \
  --force
