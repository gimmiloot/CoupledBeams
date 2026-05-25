# Thickness-Mismatch Diagnostic Model

This folder documents the initial diagnostic-only extension for two coupled
circular rods with different radii. It is separate from the current article
workflow and does not replace the baseline equal-radius determinant.

For project-wide rules on branch identity, thin-rod applicability, diagnostic
workflow, and consistency checks, see `../project_rules.md`.

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

## Branch and Mode Numbering Convention

For fixed-eta thickness-mismatch `Lambda(mu)` diagnostics, a branch means a
descendant of a mode shape, usually seeded at `mu=0`. Sorted positions are
diagnostic metadata only, and unresolved low-MAC or large-jump assignments do
not rename descendant branches. See `../project_rules.md` for the
project-wide rule.

## Thin-Rod Applicability Rule

For the mass-preserving eta parameterization, the diameter-to-length ratios are

```text
2*r1/l1 = 4*epsilon*tau1/(1 - mu)
2*r2/l2 = 4*epsilon*tau2/(1 + mu)
```

Diagnostic plots draw segments solid only when both ratios satisfy the thin-rod
diameter criterion, and dashed otherwise. See `../project_rules.md` for the
project-wide rule.

## Isolated-Rod Reference Curves

Isolated-rod reference curves are interpretation aids, not proof of branch
identity or modal exchange. Diagnostic scripts must state which boundary
condition family is used.

For the thickness-mismatch eta scaling, the isolated-rod references use the
same project normalization as the coupled system:

```text
Lambda_rod_1 = alpha_n*sqrt(tau1)/(1 - mu)
Lambda_rod_2 = alpha_n*sqrt(tau2)/(1 + mu)
```

The clamped-supported / clamped-pinned reference family uses the project
convention already present in the diagnostic scripts, with roots of
`tan(alpha)=tanh(alpha)`. The clamped-clamped / fixed-fixed reference family
uses roots of `cosh(alpha) cos(alpha)=1` and the same eta-normalized
`Lambda_rod_i` conversion above.

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

For the current script audit, preferred entry points, historical diagnostics,
and future wrapper/refactor TODOs, see
`../../scripts/analysis/thickness_mismatch/README.md`.

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
  `mu=0.9` by adjacent-step analytic shape-MAC assignment with a marked
  warning when a step contains a low-MAC assignment. Nearest-frequency
  assignment is kept only as diagnostic metadata. For
  `mu>0`, `eta=0.1` means that the longer second rod is thicker, while `eta=-0.1` means that the
  shorter first rod is thicker. Edit the user-parameter block at the top of
  this script to set one `BETA_DEG`, `EPSILON`, `ETA_VALUES`, the mu range, and
  branch/root scan settings. By default this script saves only a PNG; CSV and
  Markdown report files are optional debug outputs controlled by `SAVE_CSV` and
  `SAVE_REPORT`.
- `scripts/analysis/plot_lambda_mu_thickness_mismatch_beta15_eta_large_slenderness.py`
  is a large-eta `Lambda(mu)` diagnostic for `eta=-0.5, 0, 0.5` at
  `beta=15 deg` and `epsilon=0.0025`. It uses fixed-eta shape-MAC tracking with
  low-MAC warnings, stores warnings when nearest-frequency and MAC assignments
  differ, tracks the first 10 branches while plotting the first 7,
  and uses the common diameter-to-length helper, marking curve segments as dashed
  whenever the thin-rod criterion `2*r_i/l_i <= 0.1` is violated for either rod.
- `scripts/analysis/fem_check_thickness_mismatch_eta_p0p5_beta15.py` performs a
  diagnostic-only FEM check for `eta=0.5`, `beta=15 deg`, and `epsilon=0.0025`.
  It uses the same planar Euler-Bernoulli frame-element convention as the
  existing FEM code, but assembles a separate two-radius diagnostic model with
  arm-specific `S_i` and `J_i`; it writes CSV/PNG/Markdown outputs and branch
  5/6/7 shape-MAC diagnostics.
- `scripts/analysis/check_thickness_mismatch_branch_identity_eta_p0p5.py`
  performs a refined local descendant-branch audit for `eta=0.5`,
  `beta=15 deg`, and `epsilon=0.0025`. It seeds branches at `mu=0`, checks
  descendants 5--7 on the refined `mu=0.65..0.90` grid, and reports whether
  any apparent sorted-position change is supported or suspicious.
- `scripts/analysis/plot_thickness_mismatch_eta_p0p5_global_spectrum.py`
  builds an overview diagnostic for `eta=0.5`, `beta=15 deg`, and
  `epsilon=0.0025`, plotting the first eight sorted roots, the first eight
  descendant branches, and the accepted canonical sorted position of each
  descendant over the full `mu=0..0.9` range.
- `scripts/analysis/plot_lambda_mu_thickness_mismatch_beta15_eta_descendants.py`
  writes the current diagnostic `Lambda(mu)` eta-sweep for
  `eta=-0.5, 0, 0.5` at `beta=15 deg` and `epsilon=0.0025`, plotting the
  first six descendant branches seeded at `mu=0` with the diameter-based
  solid/dashed applicability split.
- `scripts/analysis/plot_thickness_mismatch_fem_comparison_eta_p0p5_beta15.py`
  overlays eta=0.5 analytic descendant branches with sorted FEM frequency
  markers on a sparse `mu` grid and reports same-index and nearest-FEM-mode
  frequency diagnostics without changing the FEM physical model.
- `scripts/analysis/plot_lambda_mu_thickness_mismatch_eta_p0p5_with_isolated_rods.py`
  plots eta=0.5 descendant branches together with the existing
  clamped-supported single-rod reference convention, adapted to the
  thickness-mismatch Lambda scaling as
  `Lambda=alpha*sqrt(tau_i)/(1 +/- mu)`.
- `scripts/analysis/plot_lambda_mu_thickness_mismatch_beta45_eta_p0p5_with_isolated_rods.py`
  repeats the same eta=0.5 isolated-rod diagnostic at `beta=45 deg` and
  reports a cautious nearest-reference-grid comparison against the beta-15
  case.
- `scripts/analysis/plot_lambda_mu_thickness_mismatch_beta45_eta_p0p5_with_isolated_rods_fixed_fixed.py`
  repeats the beta-45 eta=0.5 isolated-rod diagnostic with
  clamped-clamped / fixed-fixed single-rod references. It changes only the
  isolated-rod alpha roots, preserving the same eta-model Lambda scaling
  `Lambda=alpha*sqrt(tau_i)/(1 +/- mu)`.
- `scripts/analysis/plot_thickness_mismatch_branch5_shapes.py` reconstructs
  full analytic deformed shapes for tracked branch 5 in the thickness-mismatch
  model. It tracks branch 5 from `mu=0` separately for each configured eta,
  reconstructs the determinant null-vector shape, and writes diagnostic PNG
  panels only.

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
- `results/thickness_mismatch_lambda_mu_beta15_eps0p0025_eta_m0p5_0_p0p5_thickness_ratio_split.png`
- `results/thickness_mismatch_lambda_mu_beta15_eps0p0025_eta_m0p5_0_p0p5_mac_tracking_warnings.csv`
- `results/thickness_mismatch_lambda_mu_beta15_eps0p0025_eta_m0p5_0_p0p5_mac_tracking_report.md`
- `results/thickness_mismatch_fem_check_beta15_eps0p0025_eta_p0p5.csv`
- `results/thickness_mismatch_fem_check_beta15_eps0p0025_eta_p0p5.png`
- `results/thickness_mismatch_fem_check_beta15_eps0p0025_eta_p0p5_report.md`
- `results/thickness_mismatch_branch_identity_eta_p0p5_beta15_refined.csv`
- `results/thickness_mismatch_branch_identity_eta_p0p5_beta15_refined.png`
- `results/thickness_mismatch_branch_identity_eta_p0p5_beta15_refined_report.md`
- `results/thickness_mismatch_eta_p0p5_beta15_global_spectrum_8modes.csv`
- `results/thickness_mismatch_eta_p0p5_beta15_global_spectrum_8modes.png`
- `results/thickness_mismatch_eta_p0p5_beta15_global_spectrum_8modes_zoom_4_8.png`
- `results/thickness_mismatch_eta_p0p5_beta15_global_spectrum_8modes_report.md`
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

Run the large-eta `Lambda(mu)` diagnostic with the `2*r_i/l_i <= 0.1`
solid/dashed split and shape-MAC tracking with:

```bash
python scripts/analysis/plot_lambda_mu_thickness_mismatch_beta15_eta_large_slenderness.py
```

Run the eta=0.5 FEM frequency/MAC diagnostic with:

```bash
python scripts/analysis/fem_check_thickness_mismatch_eta_p0p5_beta15.py
```

Run the eta=0.5 refined branch-identity audit with:

```bash
python scripts/analysis/check_thickness_mismatch_branch_identity_eta_p0p5.py
```

Run the eta=0.5 global sorted-spectrum and descendant overview with:

```bash
python scripts/analysis/plot_thickness_mismatch_eta_p0p5_global_spectrum.py
```

Run the tracked branch-5 shape comparison diagnostic with:

```bash
python scripts/analysis/plot_thickness_mismatch_branch5_shapes.py
```
