# QUANTAS EOS SPEC 1
[metadata]
title = NaF P-V-T coupling comparison

[input]

[batch]
failure_policy = continue

[defaults]
accept = no
covariance_scaling = inflate-only
max_iterations = 10000

[defaults.pvt]
pv_model = BM3
solver = effective-variance
inner_max_iterations = 10000

[presentation]
detail = extended
show_uncertainties = yes
max_data_rows = 8

[job linear-coupling]
domain = pvt
targets = volume
pv_model = BM3
vt_model = berman:quadratic
coupling = linear
solver = effective-variance
initial.dK0_dT = -0.02
bound.dK0_dT = none : 0.0

[job anderson-gruneisen]
domain = pvt
targets = volume
pv_model = BM3
vt_model = berman:quadratic
coupling = anderson-gruneisen
solver = effective-variance
initial.delta = 4.0
bound.delta = 0.0 : 20.0

[job hp-einstein]
domain = pvt
targets = volume
pv_model = BM3
coupling = thermal-pressure
thermal_pressure_model = holland-powell-einstein
solver = effective-variance
initial.theta_e = 500.0
bound.theta_e = 1.0 : 5000.0

[job mgd-full]
domain = pvt
targets = volume
pv_model = BM4
coupling = thermal-pressure
thermal_pressure_model = mgd
volume_basis = molar-formula-unit
formula = NaF
solver = effective-variance
fix.temperature_ref = 295.0
initial.theta_d0 = 459.0
bound.theta_d0 = 1.0 : 5000.0
initial.gamma0 = 1.5
bound.gamma0 = 0.01 : 10.0
initial.q = 1.0
bound.q = -10.0 : 10.0
accept = yes
