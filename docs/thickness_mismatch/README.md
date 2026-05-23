# Thickness-Mismatch Diagnostic Model

This folder documents the initial diagnostic-only extension for two coupled
circular rods with different radii. It is separate from the current article
workflow and does not replace the baseline equal-radius determinant.

## Parameters Held Fixed

- The total length is fixed: `l1 + l2 = 2 l`.
- The length-ratio parameter remains
  `mu = (l2 - l1) / (l1 + l2)`, so `l1 = l (1 - mu)` and
  `l2 = l (1 + mu)`.
- `r0` is the common radius in the equal-radius limit `eta = 0`.
- The total mass of the two-rod system is held fixed relative to the
  equal-radius system with radius `r0`.
- The base thickness parameter is `epsilon = r0 / (2 l)`.

## Radius-Ratio Parameter

The radius mismatch is parameterized by

```text
eta = (r2 - r1) / (r1 + r2)
```

The mass-preserving factors are

```text
denom = sqrt(1 + 2 mu eta + eta^2)
tau1 = (1 - eta) / denom
tau2 = (1 + eta) / denom
r1 = r0 tau1
r2 = r0 tau2
```

For circular sections,

```text
S1 = S0 tau1^2
S2 = S0 tau2^2
J1 = J0 tau1^4
J2 = J0 tau2^4
```

The mass conservation identity is

```text
(1 - mu) tau1^2 + (1 + mu) tau2^2 = 2
```

The local thickness parameters are

```text
epsilon1 = epsilon tau1
epsilon2 = epsilon tau2
```

## Determinant Convention

The diagnostic determinant uses the same unknown ordering as the baseline
analytic model:

```text
Z = (A1, B1, A2, B2, P1, P2)
```

The bending arguments are

```text
Lambda1 = Lambda (1 - mu) / sqrt(tau1)
Lambda2 = Lambda (1 + mu) / sqrt(tau2)
```

The axial arguments are kept in terms of the base thickness parameter:

```text
theta1 = epsilon Lambda^2 (1 - mu)
theta2 = epsilon Lambda^2 (1 + mu)
```

The row structure and signs follow the baseline determinant in
`src/my_project/analytic/formulas.py`. Radius mismatch enters through the
relative weights

```text
rotation:      tau_i^(-1/2)
moment:        tau_i^3
transv. shear: tau_i^(5/2)
axial force:   tau_i^2
```

At `eta = 0`, `tau1 = tau2 = 1`, and the new determinant is expected to reduce
directly to the baseline equal-radius determinant.

## Symmetry Expectation

Exchanging the two rods corresponds to simultaneous sign reversal

```text
mu -> -mu
eta -> -eta
```

The sorted spectrum should therefore be invariant under this paired
transformation.

## Diagnostic Code

- `src/my_project/analytic/formulas_thickness_mismatch.py` implements
  `det_eta(Lambda, beta, mu, epsilon, eta)` and local root scanning helpers.
- `scripts/analysis/check_thickness_mismatch_eta_zero_limit.py` checks mass
  conservation, the `eta = 0` determinant/root limit, and swap symmetry.
- `scripts/analysis/plot_lambda_eta_thickness_mismatch.py` writes initial
  diagnostic plots and CSVs for `Lambda(eta)` and an optional `Lambda(mu)`
  eta-sweep.
- `scripts/analysis/track_lambda_eta_thickness_mismatch.py` starts branches at
  `eta=0` and continues them to positive and negative eta by unique nearest-root
  assignment. This separates branch-continuation diagnostics from the sorted
  root index shown in the first `Lambda(eta)` plot.
- `scripts/analysis/track_lambda_mu_thickness_mismatch_eta_sweep.py` starts
  branches at `mu=0` for fixed `eta=-0.1, 0, 0.1` and continues them to
  `mu=0.9` by the same nearest-root assignment. For `mu>0`, `eta=0.1` means
  that the longer second rod is thicker, while `eta=-0.1` means that the
  shorter first rod is thicker. Edit the user-parameter block at the top of
  this script to set one `BETA_DEG`, `EPSILON`, `ETA_VALUES`, the mu range, and
  branch/root scan settings. By default this script saves only a PNG; CSV and
  Markdown report files are optional debug outputs controlled by `SAVE_CSV` and
  `SAVE_REPORT`.
- `scripts/analysis/plot_thickness_mismatch_branch5_shapes.py` reconstructs
  full analytic deformed shapes for tracked branch 5 in the thickness-mismatch
  model. It tracks branch 5 from `mu=0` separately for each configured eta,
  reconstructs the determinant null-vector shape, and writes diagnostic PNG
  panels only.

Sorted root curves can be misleading near veering or close branch interactions:
the sorted index can change its branch identity even when the continued branch
varies smoothly. The tracking CSV records the nearest sorted index at each eta
or mu so these switches are visible.

## Diagnostic Outputs

The first diagnostic run writes:

- `results/thickness_mismatch_eta_zero_roots_check.csv`
- `results/thickness_mismatch_swap_symmetry_check.csv`
- `results/thickness_mismatch_lambda_eta_beta15_eps0p0025.png`
- `results/thickness_mismatch_lambda_eta_beta15_eps0p0025.csv`
- `results/thickness_mismatch_lambda_mu_beta15_eps0p0025_eta_sweep.png`
- `results/thickness_mismatch_lambda_mu_beta15_eps0p0025_eta_sweep.csv`
- `results/thickness_mismatch_lambda_eta_beta15_eps0p0025_tracked.csv`
- `results/thickness_mismatch_lambda_eta_beta15_eps0p0025_sorted_for_tracking.csv`
- `results/thickness_mismatch_lambda_eta_beta15_eps0p0025_tracked.png`
- `results/thickness_mismatch_lambda_eta_tracking_report.md`
- `results/thickness_mismatch_lambda_mu_beta{...}_eps{...}_eta_sweep_tracked.png`
- optional, only when enabled in the script:
  `results/thickness_mismatch_lambda_mu_beta{...}_eps{...}_eta_sweep_tracked.csv`
  and
  `results/thickness_mismatch_lambda_mu_beta{...}_eps{...}_eta_sweep_tracking_report.md`
- `results/thickness_mismatch_branch5_shapes_beta15_eps0p0025_mu0p0.png`
- `results/thickness_mismatch_branch5_shapes_beta15_eps0p0025_mu0p1.png`

These files are diagnostic artifacts only. They do not change the article
figures, the FEM model, or the baseline equal-radius determinant.

Run the eta-tracking diagnostic with:

```bash
python scripts/analysis/track_lambda_eta_thickness_mismatch.py
```

Run the fixed-eta `Lambda(mu)` tracking diagnostic with:

```bash
python scripts/analysis/track_lambda_mu_thickness_mismatch_eta_sweep.py
```

Run the tracked branch-5 shape comparison diagnostic with:

```bash
python scripts/analysis/plot_thickness_mismatch_branch5_shapes.py
```
