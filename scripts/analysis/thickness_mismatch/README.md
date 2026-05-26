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
| Single fixed-fixed rod EB/Timoshenko check | `python scripts/analysis/compare_single_rod_eb_timoshenko.py` | none | `results/single_rod_fixed_fixed_eb_timoshenko_frequencies.csv`, `results/single_rod_fixed_fixed_eb_timoshenko_lambda_vs_epsilon.png`, `results/single_rod_fixed_fixed_eb_timoshenko_report.md` | Preferred diagnostic for the one-rod Timoshenko check before coupled-beam use; parameterized by epsilon list and marked with the `2r/L <= 0.1` thin-rod criterion. |
| Single fixed-fixed rod 3D solid FEM workflow | `python scripts/analysis/solid_fem_single_rod_fixed_fixed.py` | optional `GMSH_EXE`, `CCX_EXE` | `docs/thickness_mismatch/solid_fem_audit.md`, `results/single_rod_fixed_fixed_3d_solid_fem_report.md`, `results/single_rod_fixed_fixed_3d_solid_fem_comparison.csv`, `results/single_rod_fixed_fixed_3d_solid_fem_mode_metrics.csv`, `results/single_rod_fixed_fixed_3d_solid_fem_bending_doublet_comparison.csv`, `results/solid_fem_single_rod/` | Diagnostic-only external solid-FEM workflow seed. Writes Gmsh/CalculiX input templates, builds fixed-end node sets for ccx, parses `.frd` displacement modes when available, classifies bending/axial/torsion-like modes, and compares classified bending doublets against the one-rod EB/Timoshenko references; external solvers remain optional. |
| Coupled equal rods multi-beta 3D solid FEM workflow | `python scripts/analysis/solid_fem_coupled_equal_rods.py` | `solid_fem_coupled_equal_rods_beta15.py` remains for reproducing the original beta=15 diagnostic outputs | `results/coupled_equal_rods_3d_solid_fem_report.md`, `results/coupled_equal_rods_3d_solid_fem_geometry_audit.csv`, `results/coupled_equal_rods_3d_solid_fem_sorted_comparison_by_beta.csv`, `results/coupled_equal_rods_3d_solid_fem_mode_metrics_by_beta.csv`, `results/coupled_equal_rods_3d_solid_fem_in_plane_comparison_by_beta.csv`, `results/solid_fem_coupled_equal_rods/` | Diagnostic-only beta sweep for equal rods at `beta=15,45,90 deg`, `mu=0`, `eta=0`; audits fused-joint overlap using `r/sin(beta)`, keeps coordinate-derived outer fixed node sets, parses `.frd` mode shapes, and compares classified in-plane solid modes with beta-specific EB/Timoshenko sorted references. |
| Coupled equal rods point-joint 3D solid FEM workflow | `python scripts/analysis/solid_fem_coupled_equal_rods_point_joint.py` | fused workflow remains the comparison baseline | `results/coupled_equal_rods_point_joint_3d_solid_fem_report.md`, `results/coupled_equal_rods_point_joint_3d_solid_fem_geometry_audit.csv`, `results/coupled_equal_rods_point_joint_3d_solid_fem_sorted_comparison_by_beta.csv`, `results/coupled_equal_rods_point_joint_3d_solid_fem_mode_metrics_by_beta.csv`, `results/coupled_equal_rods_point_joint_3d_solid_fem_in_plane_comparison_by_beta.csv`, `results/coupled_equal_rods_point_joint_3d_solid_fem_mac_matrix_eb.csv`, `results/coupled_equal_rods_point_joint_3d_solid_fem_mac_matrix_timoshenko.csv`, `results/coupled_equal_rods_point_joint_3d_solid_fem_mac_matches.csv`, `results/solid_fem_coupled_equal_rods_point_joint/` | Diagnostic-only alternative joint idealization: two separate cylinder meshes whose inner end-face node sets are coupled to one shared CalculiX rigid reference node after a local `*RIGID BODY` capability probe. Adds centerline MAC-like matching against EB/Timoshenko analytic shape reconstructions, which is preferred over simple classified in-plane ordering for mode-identity diagnostics. |
| Coupled equal rods beta=15 3D solid FEM workflow | `python scripts/analysis/solid_fem_coupled_equal_rods_beta15.py` | optional `GMSH_EXE`, `CCX_EXE` | `results/coupled_equal_rods_beta15_3d_solid_fem_report.md`, `results/coupled_equal_rods_beta15_3d_solid_fem_sorted_comparison.csv`, `results/coupled_equal_rods_beta15_3d_solid_fem_mode_metrics.csv`, `results/coupled_equal_rods_beta15_3d_solid_fem_in_plane_comparison.csv`, `results/solid_fem_coupled_equal_rods_beta15/` | Diagnostic-only 3D solid check for beta=15, mu=0, eta=0 equal rods. Generates a fused two-cylinder mesh, coordinate-derived outer fixed node sets, CalculiX modal runs, `.frd` shape classification, and classified in-plane FEM vs EB/Timoshenko sorted-frequency comparisons. |
| Coupled equal rods EB/Timoshenko check | `python scripts/analysis/compare_coupled_equal_rods_eb_timoshenko.py` | none | `results/coupled_equal_rods_beta15_eb_timoshenko_frequencies.csv`, `results/coupled_equal_rods_beta15_timoshenko_kappa_sensitivity.csv`, `results/coupled_equal_rods_beta15_eb_timoshenko_lambda_vs_epsilon.png`, `results/coupled_equal_rods_beta15_eb_timoshenko_report.md` | Preferred diagnostic for the first coupled-rod Timoshenko sanity check at `beta=15 deg`, `mu=0`, `eta=0`; compares sorted frequencies only, marks `4*epsilon <= 0.1`, records cut-off margins, and runs a small kappa sensitivity check. |
| Eta-zero and swap checks | `python scripts/analysis/check_thickness_mismatch_eta_zero_limit.py` | none | `results/thickness_mismatch_eta_zero_roots_check.csv`, `results/thickness_mismatch_swap_symmetry_check.csv` | Sanity checks for the diagnostic eta extension. |
| Descendant `Lambda(mu)` eta sweep | `python scripts/analysis/plot_lambda_mu_thickness_mismatch_beta15_eta_descendants.py` | `plot_lambda_mu_thickness_mismatch_beta15_eta_large_slenderness.py`, `track_lambda_mu_thickness_mismatch_eta_sweep.py` | `results/thickness_mismatch_lambda_mu_beta15_eps0p0025_eta_m0p5_0_p0p5_descendants.png` and report | Current clean presentation-style plot; tracking warnings stay in reports/data, not as x-markers. |
| Eta=0.5 global spectrum overview | `python scripts/analysis/plot_thickness_mismatch_eta_p0p5_global_spectrum.py` | none | `results/thickness_mismatch_eta_p0p5_beta15_global_spectrum_8modes.*` | Shows sorted roots, descendant branches, canonical sorted positions, and unresolved candidate assignments. |
| Refined branch-identity audit | `python scripts/analysis/check_thickness_mismatch_branch_identity_eta_p0p5.py` | none | `results/thickness_mismatch_branch_identity_eta_p0p5_beta15_refined.*` | Local descendant audit for branches 5--7 at eta=0.5. |
| Eta=0.5 analytic/FEM overlay | `python scripts/analysis/plot_thickness_mismatch_fem_comparison_eta_p0p5_beta15.py` | `fem_check_thickness_mismatch_eta_p0p5_beta15.py` | `results/thickness_mismatch_fem_comparison_beta15_eps0p0025_eta_p0p5.*` | Lightweight plotting comparison. The older FEM check keeps the detailed CSV/MAC diagnostics. |
| Eta=0.5 isolated-rod references, beta=15 | `python scripts/analysis/plot_lambda_mu_thickness_mismatch_eta_p0p5_with_isolated_rods.py` | none | `results/thickness_mismatch_lambda_mu_beta15_eps0p0025_eta_p0p5_with_isolated_rods.*` | Uses the clamped-supported / clamped-pinned reference convention and eta-normalized `Lambda_rod`. |
| Eta=0.5 isolated-rod references, beta=45 | `python scripts/analysis/plot_lambda_mu_thickness_mismatch_beta45_eta_p0p5_with_isolated_rods.py` | none | `results/thickness_mismatch_lambda_mu_beta45_eps0p0025_eta_p0p5_with_isolated_rods.*` | Same reference convention as beta=15, different beta only. |
| Eta=0.5 fixed-fixed isolated-rod references, beta=45 | `python scripts/analysis/plot_lambda_mu_thickness_mismatch_beta45_eta_p0p5_with_isolated_rods_fixed_fixed.py` | none | `results/thickness_mismatch_lambda_mu_beta45_eps0p0025_eta_p0p5_with_isolated_rods_fixed_fixed.*` | Uses clamped-clamped / fixed-fixed reference roots. |
| Branch-shape eta overlay | `python scripts/analysis/plot_thickness_mismatch_branch_shapes_vs_eta.py` | `plot_thickness_mismatch_branch5_shapes.py` | `results/thickness_mismatch_shapes_branch5_beta15_mu0p0_eta_m0p5_0_p0p5.png`, `results/thickness_mismatch_shapes_branch5_beta15_mu0p15_eta_m0p5_0_p0p5.png` | Parameterized descendant-shape overlay for several eta values at target mu values, using shape-MAC tracking from `mu=0`. |
| Branch-shape diagnostics | `python scripts/analysis/plot_thickness_mismatch_branch5_shapes.py` | none | `results/thickness_mismatch_branch5_shapes_*.png` | Shape-only diagnostic for tracked branch 5. |

## Current FEM Validation Status

The current diagnostic FEM validation path has three levels:

- `solid_fem_single_rod_fixed_fixed.py` is the one-rod fixed-fixed 3D solid FEM
  benchmark and currently gives the clearest support for the Timoshenko
  correction through classified bending-doublet comparisons.
- `solid_fem_coupled_equal_rods.py` is the fused-cylinder coupled-rod geometry
  audit across `beta=15,45,90 deg`; it is useful diagnostically but is not a
  direct validation of the ideal point-joint one-dimensional model.
- `solid_fem_coupled_equal_rods_point_joint.py` is the current best coupled-rod
  validation path, using separate solid cylinders with a CalculiX rigid
  reference joint and MAC-like centerline shape matching.

The status, limitations, and article-ready roadmap are tracked in
`../../../docs/thickness_mismatch/fem_validation_status.md`.

## Timoshenko Diagnostic Notes

The single-rod and coupled equal-rod EB/Timoshenko diagnostics now record the
Timoshenko cut-off frequency

```text
Omega_c = sqrt(kappa*G*A/(rho*I))
```

For the coupled equal-rod project normalization `Omega = epsilon*Lambda^2`,
with `r = 2*epsilon` and `L_segment = 1`, the report also records

```text
Lambda_c = (kappa/(2*(1 + nu)))**0.25/epsilon
```

The current diagnostic Timoshenko extension modifies only the flexural part:
- bending moment: `M = E I psi'`
- shear force: `Q = kappa G A (w' - psi)`

The axial part remains classical:
- `N = E A u'`

This is acceptable for the current low flexural-frequency diagnostics, but
higher-frequency longitudinal refinements may require Mindlin-Herrmann /
related rod models.

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
- `plot_thickness_mismatch_branch5_shapes.py` is the older branch-5 shape
  diagnostic. Prefer `plot_thickness_mismatch_branch_shapes_vs_eta.py` when
  comparing several eta values on one shape plot or changing branch/mu/eta
  parameters.

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
