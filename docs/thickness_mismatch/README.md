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

For the thickness-mismatch diagnostics, a branch is a descendant of a mode
shape, not the k-th entry of the sorted spectrum at each parameter value.

- `sorted mode k` or `k-th sorted frequency` means the k-th frequency in
  increasing order at one fixed parameter point.
- `descendant branch k` means the continuation of the k-th mode shape selected
  at the tracking start point, usually `mu=0` for fixed-eta `Lambda(mu)`
  diagnostics.
- In this diagnostic study, `branch k` means `descendant branch k` unless a
  script explicitly says it is plotting sorted roots.
- The current sorted position of a descendant may change, but that sorted
  position is diagnostic metadata. It is not the definition of branch identity.
- A change of sorted position must be supported by local mode-shape/MAC
  evidence. A jump over one or more neighboring sorted roots, such as a
  descendant recorded as moving from sorted position 5 to 7 in one step, is a
  suspicious assignment and requires a refined local mu-grid check before it is
  interpreted as a real change of sorted position.
- Low-MAC or low-margin assignments do not change the accepted descendant
  identity and do not change the accepted canonical sorted position. They must
  be recorded as unresolved candidate assignments; the previous canonical
  sorted position is retained until a refined local analysis provides a
  high-confidence continuation.

MAC diagnostics can indicate possible ambiguity in the branch group 5--7. Such
results should be reported as shape-assignment ambiguity until a refined local
grid confirms or rejects the sorted-position change.

## Thin-Rod Applicability Rule

For all diagnostic plots in the thickness-mismatch study, the applicability of
the thin-rod model is checked with the diameter-to-length ratio

```text
2*r_i/l_i <= 0.1,  i = 1, 2.
```

Here `2*r_i` is the circular-section diameter. Do not replace this rule with
the radius-to-length criterion `r_i/l_i <= 0.1`.

In the mass-preserving eta parameterization this becomes

```text
thickness_ratio_1 = 2*r1/l1 = 4*epsilon*tau1/(1 - mu)
thickness_ratio_2 = 2*r2/l2 = 4*epsilon*tau2/(1 + mu)
```

A plotted curve segment is valid and drawn solid only when both ratios are at
or below `0.1`. If either ratio exceeds `0.1`, the segment is drawn dashed and
the script or report must warn with the eta value, the affected mu range, the
rod number, and the maximum value of `2*r_i/l_i` on the computed grid.

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
  whenever the thin-rod criterion `2r_i/l_i <= 0.1` is violated for either rod.
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
- `scripts/analysis/plot_thickness_mismatch_branch5_shapes.py` reconstructs
  full analytic deformed shapes for tracked branch 5 in the thickness-mismatch
  model. It tracks branch 5 from `mu=0` separately for each configured eta,
  reconstructs the determinant null-vector shape, and writes diagnostic PNG
  panels only.

Sorted root curves can be misleading near veering or close branch interactions:
the sorted index can change its branch identity even when the continued branch
varies smoothly. Tracking CSV files record the current sorted position as
diagnostic metadata only; they do not redefine the descendant branch number.
Any large sorted-position jump must be treated as a suspicious assignment until
a refined local MAC audit supports it.

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

Run the large-eta `Lambda(mu)` diagnostic with the `2r_i/l_i <= 0.1`
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
