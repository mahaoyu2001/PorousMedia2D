# Stress-tail validation

The comparison includes direct FEM, GND with analytical moments, GND with
each realization's FEM mean/variance/kurtosis, and full compensated Campbell
inversion. The weak-field series accelerates the full exponent; it is not a
global fourth-cumulant truncation.

```matlab
addpath('output/Round2_Tail_Validation/Programs','-begin')
run_tail_validation('plots')        % bundled CSV -> figures
run_tail_validation('export')       % compact caches -> tables and figures
run_tail_validation('base')         % base inversion or valid cache
run_tail_validation('convergence')  % all six numerical configurations
```

Small result records in `Results/Cache` include probabilities, quantiles and
diagnostics. No raw fields or solver databases are included. Configuration
matching prevents reuse of a different numerical setup. To force fresh
inversion, copy the source to a separate work directory without the caches,
or move the relevant configuration files aside. The final reported setup is
`SpatialVeryFine`: radial 4096, angular 8192, Fourier grid 2^20, window 16.

`run_tail_validation('fem')` and `validate_tail_sources` require exported
raw fields in `output/FieldValue`, or the external input override. Source
byte/time checks must be repeated for independently regenerated fields.
Accuracy depends on component and threshold. Full inversion is exact only
for the chosen Poisson superposition, not the interacting FEM medium.
