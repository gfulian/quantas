# QUANTAS EOS SPEC 1
[metadata]
title = Quartz P-V solver and model comparison

[input]

[batch]
failure_policy = continue

[defaults]
accept = no
covariance_scaling = inflate-only
max_iterations = 10000

[defaults.pv]
model = BM3

[presentation]
detail = extended
show_uncertainties = yes
max_data_rows = 8

[job bm3-ols]
domain = pv
targets = volume
solver = ols

[job bm3-wls]
domain = pv
targets = volume
solver = wls

[job bm3-effective-variance]
domain = pv
targets = volume
solver = effective-variance
inner_max_iterations = 10000
accept = yes

[job bm2-effective-variance]
domain = pv
targets = volume
model = BM2
solver = effective-variance
inner_max_iterations = 10000

[job pt3-effective-variance]
domain = pv
targets = volume
model = PT3
solver = effective-variance
inner_max_iterations = 10000

[job v3-effective-variance]
domain = pv
targets = volume
model = V3
solver = effective-variance
inner_max_iterations = 10000
