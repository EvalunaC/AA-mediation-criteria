# Quantifying Mediated Treatment Effects to Strengthen Accelerated Approval Evidence in Oncology

Reproduction code for the revised manuscript (PO-26-00422). All results are
based on 500 simulation replications per scenario (and per confounding cell)
plus a real-data case study of study 20050203 (panitumumab + FOLFOX), which
ships with the `MultiStateModels` package.

## Requirements

- R >= 4.5 with packages: `rjags`, `survival`, `dplyr`, `tidyr`, `ggplot2`,
  `patchwork` (or `gridExtra`)
- JAGS 4.x
- The `MultiStateModels` package (branch `multimed`), cloned into the repo
  root and installed:

  ```sh
  git clone -b multimed https://github.com/EvalunaC/MultiStateModels.git
  R CMD INSTALL MultiStateModels
  ```

  The clone must remain present: `patches/` rebuilds two package functions
  from its source file, and the JAGS model file is read from it. Set the
  `MSM_DIR` environment variable if the clone lives elsewhere.

## Layout

```
simulation/   generates the raw model fits (JAGS; compute-intensive)
patches/      two fixes applied on load: robust MCMC initialization for
              FitJagsMstate; corrected + collapsed P01 integration and
              robust quadrature for logRRCalc
analysis/     reproduces every table and figure from the raw fits (fast)
data/         raw fits land here (not tracked; see below)
output/       tables/ and figures/ (not tracked)
```

## Regenerating the raw data

From `simulation/` (each script is resumable; one `.rds` per replication;
roughly 1-4 minutes of CPU per fit):

```sh
bash generate_maturity.sh          # 6 scenarios x 500 reps x 4 interim cuts -> data/maturity     (3,000 files)
bash generate_sample_size.sh       # n in {100,200,300} x 6 x 500, month 24  -> data/sample_size  (9,000 files)
bash generate_confounding.sh       # 6 scenarios x 4 u-levels x 500 reps     -> data/confounding  (12,000 files)
bash generate_case.sh              # 2 cohorts x 6 interim cuts              -> data/case_interim (12 files)
```

Set `RSCRIPT=/path/to/Rscript` to choose the R build; pass core counts as
arguments (see each script). Seeds are fixed per (scenario, cell,
replication), so regeneration is deterministic. The raw fits used for the
manuscript are also available from the authors on request.

## Reproducing the results

```sh
cd analysis
Rscript run_all.R
```

Reads `data/`, writes every manuscript and supplementary table and figure to
`output/`, and verifies key manuscript numbers with hard assertions
(`stopifnot`) along the way. No MCMC is run at this stage; the full pipeline
takes minutes.

The decision criteria and the calibrated thresholds
(Delta* = 0.825/0.835/0.830/0.775 at months 12/18/24/30) are defined once, in
`analysis/00_criteria.R`; every analysis decides through those functions.
