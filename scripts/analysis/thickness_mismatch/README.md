# Thickness-Mismatch Script Map

This directory is a navigation layer for the diagnostic-only
thickness-mismatch study. The runnable scripts still live in the flat
`scripts/analysis/` layout for compatibility. Moving them in this pass would
require broad import and output-path churn without changing any diagnostic
result.

Project-wide rules for descendant branch identity, sorted-position metadata,
low-MAC assignments, thin-rod applicability, and diagnostic/article separation
live in `../../../docs/project_rules.md`.

## Preferred Entry Points

| Task | Preferred command | Compatibility or historical scripts | Main outputs | Notes |
| --- | --- | --- | --- | --- |
| Eta-zero and swap checks | `python scripts/analysis/check_thickness_mismatch_eta_zero_limit.py` | none | `results/thickness_mismatch_eta_zero_roots_check.csv`, `results/thickness_mismatch_swap_symmetry_check.csv` | Sanity checks for the diagnostic eta extension. |
| Descendant `Lambda(mu)` eta sweep | `python scripts/analysis/plot_lambda_mu_thickness_mismatch_beta15_eta_descendants.py` | `plot_lambda_mu_thickness_mismatch_beta15_eta_large_slenderness.py`, `track_lambda_mu_thickness_mismatch_eta_sweep.py` | `results/thickness_mismatch_lambda_mu_beta15_eps0p0025_eta_m0p5_0_p0p5_descendants.png` and report | Current clean presentation-style plot; tracking warnings stay in reports/data, not as x-markers. |
| Eta=0.5 global spectrum overview | `python scripts/analysis/plot_thickness_mismatch_eta_p0p5_global_spectrum.py` | none | `results/thickness_mismatch_eta_p0p5_beta15_global_spectrum_8modes.*` | Shows sorted roots, descendant branches, canonical sorted positions, and unresolved candidate assignments. |
| Refined branch-identity audit | `python scripts/analysis/check_thickness_mismatch_branch_identity_eta_p0p5.py` | none | `results/thickness_mismatch_branch_identity_eta_p0p5_beta15_refined.*` | Local descendant audit for branches 5--7 at eta=0.5. |
| Eta=0.5 analytic/FEM overlay | `python scripts/analysis/plot_thickness_mismatch_fem_comparison_eta_p0p5_beta15.py` | `fem_check_thickness_mismatch_eta_p0p5_beta15.py` | `results/thickness_mismatch_fem_comparison_beta15_eps0p0025_eta_p0p5.*` | Lightweight plotting comparison. The older FEM check keeps the detailed CSV/MAC diagnostics. |
| Eta=0.5 isolated-rod references, beta=15 | `python scripts/analysis/plot_lambda_mu_thickness_mismatch_eta_p0p5_with_isolated_rods.py` | none | `results/thickness_mismatch_lambda_mu_beta15_eps0p0025_eta_p0p5_with_isolated_rods.*` | Uses the clamped-supported / clamped-pinned reference convention and eta-normalized `Lambda_rod`. |
| Eta=0.5 isolated-rod references, beta=45 | `python scripts/analysis/plot_lambda_mu_thickness_mismatch_beta45_eta_p0p5_with_isolated_rods.py` | none | `results/thickness_mismatch_lambda_mu_beta45_eps0p0025_eta_p0p5_with_isolated_rods.*` | Same reference convention as beta=15, different beta only. |
| Eta=0.5 fixed-fixed isolated-rod references, beta=45 | `python scripts/analysis/plot_lambda_mu_thickness_mismatch_beta45_eta_p0p5_with_isolated_rods_fixed_fixed.py` | none | `results/thickness_mismatch_lambda_mu_beta45_eps0p0025_eta_p0p5_with_isolated_rods_fixed_fixed.*` | Uses clamped-clamped / fixed-fixed reference roots. |
| Branch-shape diagnostics | `python scripts/analysis/plot_thickness_mismatch_branch5_shapes.py` | none | `results/thickness_mismatch_branch5_shapes_*.png` | Shape-only diagnostic for tracked branch 5. |

## Historical And One-Off Diagnostics

- `plot_lambda_eta_thickness_mismatch.py` is the initial sorted-root
  `Lambda(eta)` diagnostic.
- `track_lambda_eta_thickness_mismatch.py` is the initial eta-continuation
  diagnostic.
- `track_lambda_mu_thickness_mismatch_eta_sweep.py` remains useful as an
  editable single-beta tracking sandbox, but the current clean eta sweep plot
  is `plot_lambda_mu_thickness_mismatch_beta15_eta_descendants.py`.
- `plot_lambda_mu_thickness_mismatch_beta15_eta_large_slenderness.py` is kept
  as the earlier large-eta diagnostic with warning artifacts.
- `fem_check_thickness_mismatch_eta_p0p5_beta15.py` is the heavier
  analytic-vs-FEM diagnostic with detailed CSV/MAC output.

## Helper Sources

- `src/my_project/analytic/formulas_thickness_mismatch.py` implements the
  diagnostic determinant and eta model helpers.
- `scripts/lib/thickness_mismatch_mac_tracking.py` implements diagnostic
  descendant tracking and unresolved-assignment flags.
- `scripts/lib/thickness_mismatch_diagnostic_helpers.py` collects shared
  plotting, validity, and isolated-rod reference utilities.

The CS/CP reference convention uses roots of `tan(alpha)=tanh(alpha)`. The
CC/FF reference convention uses roots of `cosh(alpha) cos(alpha)=1`. Both use
the eta-model normalization
`Lambda_rod_i = alpha_n*sqrt(tau_i)/(1 +/- mu)`.

## Future Refactor TODO

When a broader import cleanup is explicitly requested, the flat scripts can be
converted into parameterized entry points under this directory:

- `plot_descendant_eta_sweep.py`
- `plot_with_isolated_rods.py`
- `compare_fem.py`
- `audit_branch_identity.py`
- `check_eta_zero_limit.py`

The old script names should then remain as compatibility wrappers that print a
short deprecation note and call the new entry point with the old parameters and
output paths.
