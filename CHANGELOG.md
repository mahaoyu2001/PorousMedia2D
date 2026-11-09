# Changelog

## 2.0.0 - 2026-11-09

- Align analytical moments with the second-revision matrix-conditioned
  Poisson model: variance `f*K_i2/(1-f)^2`, skewness
  `K_i3/(sqrt(f)*K_i2^(3/2))`, kurtosis `3+K_i4/(f*K_i2^2)`.
- Retain isolated-pore coefficients and background `Sigma/(1-f)`.
  The dilute identification is `rho*A0 approximately f`.
- Evaluate coefficients numerically without Symbolic Math Toolbox;
  include the optional symbolic derivation of Eq. (17).
- Add corrected moment, PDF/CDF, JS/KL, exceedance and quantile programs,
  processed FEM inputs, compact metrics, and standalone FIG templates.
- Add full compensated Campbell inversion, FEM-moment GND comparisons,
  numerical sensitivity results, and hard-core pair-prefactor verification.
- Use relative paths and accept field tables with or without text headers.
- Preserve the original Abaqus workflow, coordinates and MIT license.
  This source update does not perform new FEM solves.

## Earlier public source

Initial public commit `a0cb34a` and MIT addition `1bc6179` are retained.
Older analytical moment phase factors are superseded by version 2.0.0.
