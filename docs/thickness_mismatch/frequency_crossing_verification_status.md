# Frequency Crossing Verification Status

Last updated: 2026-06-02.

This document records the diagnostic-only status of the recent positive-gap
checks for possible real frequency crossings in the Euler-Bernoulli analytic
thickness-mismatch model. It is a documentation/status snapshot based on
already generated outputs; no FEM, Gmsh, CalculiX, 3D workflow, article
figure, article text, old determinant, or old solver is part of this status
step.

## Scope

The checks were designed to separate four different questions:

- whether adjacent eigenvalue branches have a true crossing or a positive
  gap;
- whether a close approach is a true crossing or an avoided crossing;
- whether sorted-root gap analysis agrees with, or differs from, descendant
  label tracking;
- whether the conclusion is specific to the diagnostic Euler-Bernoulli
  analytic thickness-mismatch model.

The primary classification is based on adjacent sorted-root gaps. Descendant
labels are tracked separately as diagnostic metadata because branch labels can
be sensitive near symmetric or nearly degenerate parameter values.

## Eta Equals Zero, Positive Mu

Source outputs:

- `results/eta0_eps0p0025_mu_positive_gap_verification_summary.csv`
- `results/eta0_eps0p0025_mu_positive_gap_verification_details.csv`
- `results/eta0_eps0p0025_mu_positive_gap_verification_report.md`
- `results/eta0_eps0p0025_mu_positive_gap_scaling.png`
- `results/eta0_eps0p0025_mu_positive_gap_min_beta.png`
- `results/eta0_eps0p0025_mu_positive_gap_lambda_beta_selected.png`

Parameters and tested values:

- `eta = 0`
- `epsilon = 0.0025`
- `beta` scanned over the audit range
- adjacent sorted pairs `1-2` through `6-7`
- tested positive `mu` values:
  `1e-6`, `2e-6`, `5e-6`, `1e-5`, `2e-5`, `5e-5`,
  `1e-4`, `2e-4`, `5e-4`, `1e-3`, `2e-3`, `5e-3`,
  `1e-2`, `0.02`, `0.05`, `0.1`, `0.2`, `0.5`, `0.9`

Status:

- all `19/19` positive-mu cases were classified as `resolved_positive_gap`;
- no true crossing was found for tested `mu > 0`;
- smallest gap: `0.000131412306811`;
- smallest-gap point: `mu = 1e-6`, `beta = 13.3082838569 deg`,
  sorted pair `5-6`;
- all resolved points scaling estimate:
  `gap ~= 1.55525 * mu^0.737954`;
- small-mu subset, `mu <= 1e-2`, scaling estimate:
  `gap ~= 1.09713 * mu^0.706164`.

Interpretation:

- `mu = 0` is a special symmetric/degenerate case;
- in the tested `mu > 0` cases, the degeneracy unfolds into resolved positive
  gaps, interpreted as avoided crossings;
- the earlier coarse-scan rearrangement signal near `mu = 0.001..0.002` was
  rejected as a true eigenvalue crossing by the strict positive-gap
  verification.

## Eta Nonzero In Tested Ranges

Source outputs:

- `results/eta_parameter_positive_gap_verification_summary.csv`
- `results/eta_parameter_positive_gap_verification_details.csv`
- `results/eta_parameter_positive_gap_verification_report.md`
- `results/eta_parameter_positive_gap_min_gap_vs_eta.png`
- `results/eta_parameter_positive_gap_min_gap_vs_beta_eta_scan.png`
- `results/eta_parameter_positive_gap_min_gap_vs_mu_eta_scan.png`
- `results/eta_parameter_positive_gap_classification_map.png`
- `results/eta_parameter_positive_gap_lambda_eta_selected.png`
- `results/eta_parameter_positive_gap_lambda_beta_selected.png`
- `results/eta_parameter_positive_gap_lambda_mu_selected.png`

Status:

- checked `47` cases total: `8` eta-scans, `15` beta-scans, and `24`
  mu-scans;
- all tested cases were classified as `resolved_positive_gap`;
- true crossings found: none;
- unresolved possible crossings: none;
- smallest eta-nonzero gap: `4.07469969166e-4`;
- smallest-gap point: `beta = 5 deg`, `eta = 0.5`,
  `mu = 0.158647263871`, sorted pair `1-2`.

Special checkpoints:

- `beta = 15 deg`, `eta = 0.5`, sorted pair `4-5`:
  at `mu = 0.379`, gap `0.0126731097004`;
- `beta = 15 deg`, `eta = 0.5`, sorted pair `5-6`:
  at `mu = 0.716`, gap `0.0287236957149`;
- local refinement found positive minima near both checkpoints;
- `beta = 5 deg`, `eta = 0.5`, article false-crossing checkpoint:
  descendant `1-2` gap `0.00040816999683`, sign change `no`,
  sorted positions `1/2`.

Artifact handling:

- the raw default sign scan with step `0.01` was found to miss close roots
  near `mu = 0.15864`;
- strict local root repairs and refined local `mu` windows were added in the
  diagnostic audit;
- the final run recorded `150` strict local grid repairs;
- after local refinement, tracking had no unresolved cases, and MAC tracking
  was stable in the special windows.

## Safe Statement

For the tested parameter cases, no true eigenvalue crossings were found for
eta != 0 or for eta = 0 with mu > 0. Close approaches were resolved as
positive gaps and are interpreted as avoided crossings. The exactly symmetric
cases eta=0 and/or mu=0 require separate treatment because branch labels can
be symmetry-sensitive.

## Unsafe Statement

This diagnostic audit does not prove analytically that crossings are
impossible for all continuous parameter values. It supports the numerical
statement only for the tested ranges and local refinements.

## Article Implication

- Article figures should not show apparent crossings caused by coarse tracking
  or missed close roots.
- If an apparent crossing appears, it should be checked by positive-gap
  verification before being interpreted physically.
- Article text should avoid absolute statements such as "there are no
  crossings for all parameters".
- Prefer wording such as "within the tested parameter range and numerical
  tolerance".

## Protected Areas

This status document does not change, regenerate, or validate through:

- `paper_thickness_mismatch_article/main.tex`
- article figures
- `paper_dorofeev_style/`
- `src/my_project/analytic/formulas.py`
- old determinants or old solvers
- FEM models
- Gmsh or CalculiX workflows
- baseline results
