# PorousMedia2D

Version **2.0.0** accompanies the second revision of manuscript
**AMMOD-D-26-03080**, *Higher-order statistical characterization of stress
heterogeneity induced by random pores*. Versioned source:
https://github.com/mahaoyu2001/PorousMedia2D/releases/tag/v2.0.0.

The package includes Abaqus scripts, fixed pore-center inputs, material
parameters, corrected MATLAB programs and compact processed FEM results.
Large CAE/ODB files and raw fields are excluded and can be regenerated.

## Analytical convention

The isolated-pore coefficients `K_in` are unchanged. With the dilute
identification `eta=rho*A0 approximately f`, the matrix-conditioned Poisson
disturbance model uses:

```text
mean_i     = Sigma_i/(1-f)
variance_i = f*K_i2/(1-f)^2
skewness_i = K_i3/(sqrt(f)*K_i2^(3/2))
kurtosis_i = 3 + K_i4/(f*K_i2^2)
```

The background `Sigma/(1-f)` remains a mean-field approximation. `K_i3=0`
does not establish exact FEM symmetry. Accuracy depends on porosity, component
and statistic. The Gaussian predictions use the analytical mean and variance
shared with GND; Fig. 2(c)'s descriptive reference instead uses FEM moments.

## Layout

```text
MaterialData/                       matrix parameters
ModelData/                          fixed pore-center inputs
python script/Modeling/             Abaqus geometry, mesh, loading, jobs
python script/PostProcessing/       ODB -> S11 S22 S12 EVOL export
output/cal_stress_statistics.m      corrected analytical moments
output/Round2_Minimal_Correction/    moment/distribution/tail workflows
output/Round2_Tail_Validation/       Campbell and four-way tail comparisons
output/Pair_Correlation_Prefactor/   hard-core prefactor FFT check
output/Symbolic_Derivation/          optional Eq. (17) symbolic derivation
output/Style_Templates/             original standalone figure settings
tests/                             release verification
```

## Software and quick reproduction

- Abaqus/CAE 2022 and Standard for FEM. Scripts retain Python 2.7-compatible
  syntax and must use the Abaqus interpreter.
- MATLAB R2023a for the tested post-processing workflow. Core coefficients,
  moments, GND and Campbell functions need no Symbolic Math Toolbox; only
  the optional symbolic derivation requires it.

Run from the repository root in MATLAB:

```matlab
addpath('tests')
test_release
addpath('output/Round2_Minimal_Correction/Programs','-begin')
run_round2('check')
run_round2('tables')
run_round2('single')
run_round2('distributions')
addpath('output/Round2_Tail_Validation/Programs','-begin')
run_tail_validation('plots')
```

These commands use bundled processed statistics and compact field metrics;
they do not run Abaqus. Workflow READMEs explain raw-field processing and
Campbell inversion. Paths default to this repository. Optionally set
`POROUSMEDIA_OUTPUT_DIR` to an external FEM output directory containing the
input workbooks and `FieldValue` folder.

## FEM reproduction

Five fixed realizations a-e are supplied at target geometric porosities
0.5%, 1%, 2%, 3%, 5%, 7.5%, and 10%. Hard-core inputs exclude overlap.
Unconstrained (`r0`) inputs permit overlap/coalescence and use more pores
to compensate expected coverage loss. Fixed-count numerical realizations
are not identical to an infinite stationary Poisson process. Supplementary
Material 3 reports realized-porosity differences at common target values.

Edit **USER SETTINGS** in `python script/Modeling/Run2D_Geom_Batch.py`.
The default is five hard-core 5% realizations; select base name
`2DVoid_Vol5_r0` for the overlap-allowed case. From the repository root:

```powershell
abaqus cae noGUI="python script\Modeling\Run2D_Geom_Batch.py"
```

The runner builds geometry/mesh checkpoints and input files. Set
`SUBMIT_JOBS=True` explicitly to solve. Settings are a unit-square RVE,
pore radius 0.005, nominal CPE3 mesh size 0.0004, small-deformation plane
strain and uniform normal tractions. Inputs include periodic pore copies;
the finite RVE solve uses the stated traction conditions. Pores may intersect
the boundary. A coarser installation test changes the resulting statistics.

Set matching model names in the field exporter, then run:

```powershell
abaqus python "python script\PostProcessing\ExtractLE_Vol.py"
```

Each `output/FieldValue/<job>_S2.txt` contains `S11 S22 S12 EVOL`, with
integration-area weight EVOL. Corrected readers accept its optional header.
The two original comparison entry scripts in `output` remain available and
call corrected core moments. Their `auto` mode uses caches without raw files.

## Data, citation and license

Measured-statistics workbooks, per-realization metrics and CSV comparisons
support reproduction without large ODB files. Reference statistics for
additional loading/15-20% studies do not imply their raw solver outputs are
included. SD is sample scatter across five realizations, not a confidence
interval. Relative errors use the FEM mean as reference.

Cite tag **v2.0.0** to identify this source version. `CITATION.cff` supplies
metadata. This release does not claim a separately registered DOI.
License: MIT.
