# QUANTAS EOS SPEC 1
[metadata]
title = Rutile V-T solver and model comparison

[input]

[batch]
failure_policy = continue

[defaults]
accept = no
covariance_scaling = inflate-only
max_iterations = 10000

[defaults.vt]
model = berman:quadratic

[presentation]
detail = extended
show_uncertainties = yes
max_data_rows = 8

[job berman-ols]
domain = vt
targets = volume
solver = ols

[job berman-wls]
domain = vt
targets = volume
solver = wls

[job berman-effective-variance]
domain = vt
targets = volume
solver = effective-variance
inner_max_iterations = 10000
accept = yes

[job fei-effective-variance]
domain = vt
targets = volume
model = fei:inverse-square
solver = effective-variance
inner_max_iterations = 10000

[job salje-effective-variance]
domain = vt
targets = volume
model = salje
solver = effective-variance
inner_max_iterations = 10000

[job axial-a]
domain = vt
targets = a
model = berman:quadratic
solver = effective-variance
inner_max_iterations = 10000

[job axial-c]
domain = vt
targets = c
model = berman:quadratic
solver = effective-variance
inner_max_iterations = 10000
