# PorousMedia2D

This repository contains the minimal Abaqus–MATLAB workflow used to reproduce
representative two-dimensional porous-media examples and the pore-exclusion
sensitivity analysis reported in the accompanying work.

The public package intentionally contains source code, fixed pore-centre input
data, material parameters, analytical coefficients, and compact processed
statistics. Large generated files such as `.cae`, `.odb`, and element-level
field tables are excluded and can be regenerated with the workflow below.

## Reproducible examples

1. **Hard-core FEM model**: circular pore centres have a minimum separation of
   `2a`, where `a` is the pore radius.
2. **Poisson (`r_min = 0`) FEM model**: no minimum centre separation is imposed
   and circular pores may overlap.

Five fixed realizations (`a`–`e`) are provided for porosities 0.5%, 1%, 2%, 3%,
5%, 7.5%, and 10%. The files are fixed numerical inputs rather than random
seeds, so every user starts from the same geometry.

## Repository layout

```text
PorousMedia2D/
├── MaterialData/                 material properties (`key value` format)
├── ModelData/                    pore-centre coordinates (`x y` per row)
├── python script/
│   ├── Modeling/                 geometry, mesh, loading, and job creation
│   └── PostProcessing/           ODB-to-text field export
└── output/                       MATLAB statistics and plotting scripts
```

## Software

- Abaqus/CAE 2022 (Abaqus Python and the Standard solver)
- MATLAB R2020b or newer
- MATLAB Symbolic Math Toolbox for the analytical prediction

The Abaqus scripts retain Python 2.7-compatible syntax. They must be run with
the Abaqus interpreter, not a system Python installation.

## 1. Build, mesh, load, and solve the RVE models

Open `python script/Modeling/Run2D_Geom_Batch.py` and edit only its **USER
SETTINGS** block. The default reproduces the five hard-core 5% realizations:

```python
MODEL_NAME_LIST = ['2DVoid_Vol5']
LABEL_LIST = ['a', 'b', 'c', 'd', 'e']
```

For the overlap-allowed comparison, change the base name to:

```python
MODEL_NAME_LIST = ['2DVoid_Vol5_r0']
```

From the repository root, run:

```powershell
abaqus cae noGUI="python script\Modeling\Run2D_Geom_Batch.py"
```

The CAE window is not required. By default, the script creates geometry and
mesh checkpoints, writes Abaqus input files, and saves one final CAE database
containing all requested labels. Set `SUBMIT_JOBS = True` to launch the solver.

The reported model settings are defined in one place in the runner:

- RVE size: `1.0 x 1.0`
- pore radius: `0.005`
- CPE3 mesh size: `0.0004`
- vertical resultant: `1.0 x RVE width`
- plane-strain, small-deformation static analysis

The reported mesh is computationally expensive. A larger mesh size can be used
for a quick installation test, but its statistics should not be compared with
the bundled reference results.

If execution stops after geometry or meshing, set `START_STAGE` to `mesh` or
`load`, respectively. The runner then opens the corresponding checkpoint
instead of repeating the completed stages.

## 2. Export element fields

After the jobs finish, set the same model names in
`python script/PostProcessing/ExtractLE_Vol.py`, then run:

```powershell
abaqus python "python script\PostProcessing\ExtractLE_Vol.py"
```

For each ODB, the script writes
`output/FieldValue/<job-name>_S2.txt`. Its columns are:

```text
S11  S22  S12  EVOL
```

`EVOL` is the element area used as the statistical weight.

## 3. Reproduce the figures

Run either MATLAB entry script:

- `output/Deformation_Statistic_Nmoment_Geom_r0_Comparison_ErrorBand.m`
  compares analytical predictions, hard-core FEM, and Poisson FEM over all
  supplied porosities.
- `output/Deformation_Statistic_Nmoment_MinDistance.m` isolates the 5% case and
  plots variance and kurtosis versus normalized minimum centre distance.

Both scripts expose `DATA_SOURCE = 'auto' | 'field' | 'cache'` near the top.
`auto` recomputes statistics when the complete raw field set exists; otherwise
it reads the bundled processed Excel workbook. `field` also refreshes that
workbook. Figures are written to `output/pic/` as PNG and MATLAB FIG files.

## Input naming convention

`ModelData/<base-name>_<label>.txt` must match the names configured in the
Abaqus runner. For example:

```text
2DVoid_Vol5_a.txt       hard-core, 5%, realization a
2DVoid_Vol5_r0_a.txt    r_min = 0, 5%, realization a
```

All quantities use one internally consistent unit system. With the supplied
material file, stress is expressed in the same unit as its elastic modulus.

## Generated files

The `.gitignore` excludes solver databases, checkpoints, logs, extracted raw
field tables, and generated figures. The compact processed workbooks in
`output/` are intentionally versioned so the published plots can be recreated
without distributing very large ODB files.

## Citation and license

Please cite the associated paper and the archived repository record. Before
public release, add the final paper citation, Zenodo DOI, authors, and the
license selected by the repository owner.
