# Thickness-Mismatch Script Map

This directory is a navigation layer for the diagnostic-only
thickness-mismatch study. The runnable scripts still live in the flat
`scripts/analysis/` layout for compatibility. Moving them in this pass would
require broad import and output-path churn without changing any diagnostic
result.

Project-wide rules for descendant branch identity, sorted-position metadata,
low-MAC assignments, thin-rod applicability, and diagnostic/article separation
live in `../../../docs/project_rules.md`.

The current positive-gap status for frequency-crossing checks is summarized in
`../../../docs/thickness_mismatch/frequency_crossing_verification_status.md`.

The out-of-plane EB + Saint-Venant torsion subsystem is documented, covered by
lightweight determinant sanity tests, and now has a diagnostic 1D EB+torsion
FEM validation audit listed below.

## Directory Scaffold

New thickness-mismatch diagnostic scripts should be added under this directory
by diagnostic role:

```text
maps/
  eta-mu-beta frequency maps, heatmaps, sensitivity maps

shapes/
  mode-shape plotting, sorted-mode shapes, descendant-branch shapes

audits/
  local checks, sorted-vs-descendant audits, gap checks

lambda_mu/
  Lambda(mu) diagnostic plots, single-rod spectra overlays

postprocess/
  CSV-only summaries and global trend analysis
```

Existing flat scripts in `scripts/analysis/` remain supported for
compatibility. Future migration should use thin wrappers at the old paths, not
breaking moves.

Diagnostic scripts must explicitly state whether they use sorted frequencies
or descendant branches. Sorted frequencies are used for spectral maps and gap
metrics; descendant branches are used for tracking modal character.

## Smoke Mode Convention

Heavy diagnostic entry points should support `--smoke` or document an existing
`--quick` mode as the smoke mode. Smoke mode means a tiny grid, very fast
runtime, output under `results/_smoke/` or `results/smoke/`, and no overwrite
of normal diagnostic outputs unless a future script explicitly implements and
requires `--force`.

The current eta-mu-beta map, sorted `Lambda(mu)` reference plot, descendant
mode-shape beta scan, and fixed sorted-mode beta scan all expose `--smoke`.

## Mode-Shape Refactor TODO

`plot_mode_shapes_eta_beta_scan.py` and
`plot_mode_shapes_eta_beta_scan_sorted_modes.py` should eventually become one
parameterized entry point with `selection = sorted | descendant`. Shared
functionality to extract later includes beta grid construction, root solving,
coefficient extraction, geometry reconstruction, sign normalization,
grid/contact-sheet plotting, and CSV summary writing.

Example future descendant command:

```bash
python scripts/analysis/thickness_mismatch/shapes/plot_beta_scan.py \
  --eta 0.1 \
  --mu 0.3 \
  --epsilon 0.0025 \
  --beta-start 0 \
  --beta-end 11 \
  --selection descendant \
  --branch 2
```

Example future sorted-frequency command:

```bash
python scripts/analysis/thickness_mismatch/shapes/plot_beta_scan.py \
  --eta 0.1 \
  --mu 0.3 \
  --epsilon 0.0025 \
  --beta-start 0 \
  --beta-end 15 \
  --selection sorted \
  --sorted-modes 2 3
```

## Preferred Entry Points

| Task | Preferred command | Compatibility or historical scripts | Main outputs | Notes |
| --- | --- | --- | --- | --- |
| Out-of-plane EB+torsion analytic vs 1D FEM validation | `python scripts/analysis/thickness_mismatch/audits/compare_out_of_plane_analytic_vs_1d_fem.py` | none | `results/out_of_plane_fem_validation/out_of_plane_analytic_vs_1d_fem_comparison.csv`, `results/out_of_plane_fem_validation/out_of_plane_analytic_vs_1d_fem_convergence.csv`, `results/out_of_plane_fem_validation/out_of_plane_analytic_vs_1d_fem_report.md` | Diagnostic-only sorted-frequency validation for the out-of-plane determinant against an independent one-dimensional Euler--Bernoulli plus Saint-Venant torsion FEM. It checks determinant signs, joint coupling, and Lambda scaling within the same continuum model; it does not implement or run 3D FEM and does not use descendant tracking. |
| Single fixed-fixed rod EB/Timoshenko check | `python scripts/analysis/compare_single_rod_eb_timoshenko.py` | none | `results/single_rod_fixed_fixed_eb_timoshenko_frequencies.csv`, `results/single_rod_fixed_fixed_eb_timoshenko_lambda_vs_epsilon.png`, `results/single_rod_fixed_fixed_eb_timoshenko_report.md` | Preferred diagnostic for the one-rod Timoshenko check before coupled-beam use; parameterized by epsilon list and marked with the `2r/L <= 0.1` thin-rod criterion. |
| Single fixed-fixed rod 3D solid FEM workflow | `python scripts/analysis/solid_fem_single_rod_fixed_fixed.py` | optional `GMSH_EXE`, `CCX_EXE` | `docs/thickness_mismatch/solid_fem_audit.md`, `results/single_rod_fixed_fixed_3d_solid_fem_report.md`, `results/single_rod_fixed_fixed_3d_solid_fem_comparison.csv`, `results/single_rod_fixed_fixed_3d_solid_fem_mode_metrics.csv`, `results/single_rod_fixed_fixed_3d_solid_fem_bending_doublet_comparison.csv`, `results/solid_fem_single_rod/` | Diagnostic-only external solid-FEM workflow seed. Writes Gmsh/CalculiX input templates, builds fixed-end node sets for ccx, parses `.frd` displacement modes when available, classifies bending/axial/torsion-like modes, and compares classified bending doublets against the one-rod EB/Timoshenko references; external solvers remain optional. |
| Coupled equal rods multi-beta 3D solid FEM workflow | `python scripts/analysis/solid_fem_coupled_equal_rods.py` | `solid_fem_coupled_equal_rods_beta15.py` remains for reproducing the original beta=15 diagnostic outputs | `results/coupled_equal_rods_3d_solid_fem_report.md`, `results/coupled_equal_rods_3d_solid_fem_geometry_audit.csv`, `results/coupled_equal_rods_3d_solid_fem_sorted_comparison_by_beta.csv`, `results/coupled_equal_rods_3d_solid_fem_mode_metrics_by_beta.csv`, `results/coupled_equal_rods_3d_solid_fem_in_plane_comparison_by_beta.csv`, `results/solid_fem_coupled_equal_rods/` | Diagnostic-only beta sweep for equal rods at `beta=15,45,90 deg`, `mu=0`, `eta=0`; audits fused-joint overlap using `r/sin(beta)`, keeps coordinate-derived outer fixed node sets, parses `.frd` mode shapes, and compares classified in-plane solid modes with beta-specific EB/Timoshenko sorted references. |
| Coupled equal rods point-joint 3D solid FEM workflow | `python scripts/analysis/solid_fem_coupled_equal_rods_point_joint.py` | fused workflow remains the comparison baseline | Baseline mode (`MESH_SIZE_FACTORS=[1.0]`): `results/coupled_equal_rods_point_joint_3d_solid_fem_report.md`, `results/coupled_equal_rods_point_joint_3d_solid_fem_geometry_audit.csv`, `results/coupled_equal_rods_point_joint_3d_solid_fem_sorted_comparison_by_beta.csv`, `results/coupled_equal_rods_point_joint_3d_solid_fem_mode_metrics_by_beta.csv`, `results/coupled_equal_rods_point_joint_3d_solid_fem_in_plane_comparison_by_beta.csv`, `results/coupled_equal_rods_point_joint_3d_solid_fem_mac_matrix_eb.csv`, `results/coupled_equal_rods_point_joint_3d_solid_fem_mac_matrix_timoshenko.csv`, `results/coupled_equal_rods_point_joint_3d_solid_fem_mac_matches.csv`, `results/solid_fem_coupled_equal_rods_point_joint/`; mesh-convergence mode: `results/coupled_equal_rods_point_joint_mesh_convergence_report.md`, `results/coupled_equal_rods_point_joint_mesh_convergence_frequencies.csv`, `results/coupled_equal_rods_point_joint_mesh_convergence_mac.csv`, `results/coupled_equal_rods_point_joint_mesh_convergence_summary.csv`, `results/solid_fem_coupled_equal_rods_point_joint_mesh_convergence/`; planar-constrained mode: `results/coupled_equal_rods_point_joint_planar_3d_solid_fem_report.md`, `results/coupled_equal_rods_point_joint_planar_3d_solid_fem_sorted_comparison.csv`, `results/coupled_equal_rods_point_joint_planar_3d_solid_fem_mode_metrics.csv`, `results/coupled_equal_rods_point_joint_planar_3d_solid_fem_mac_matrix_eb.csv`, `results/coupled_equal_rods_point_joint_planar_3d_solid_fem_mac_matrix_timoshenko.csv`, `results/coupled_equal_rods_point_joint_planar_3d_solid_fem_mac_matches.csv`, `results/solid_fem_coupled_equal_rods_point_joint_planar/` | Diagnostic-only alternative joint idealization: two separate cylinder meshes whose inner end-face node sets are coupled to one shared CalculiX rigid reference node after a local `*RIGID BODY` capability probe. Adds centerline MAC-like matching against EB/Timoshenko analytic shape reconstructions, which is preferred over simple classified in-plane ordering for mode-identity diagnostics. Optional in-script partial-run filters `RUN_BETA_DEG_VALUES` and `RUN_EPSILON_VALUES` restrict reruns; multiple `MESH_SIZE_FACTORS` route cases into separate `mesh_<factor>` convergence folders; planar-constrained mode is an auxiliary planar-subspace diagnostic and does not replace full 3D validation. |
| Equal-thickness coupled rods `Lambda(mu)` EB/Timoshenko/FEM diagnostic | `python scripts/analysis/plot_coupled_equal_rods_beta15_eps0p01_eb_timoshenko_fem_vs_mu.py` | earlier uncorrected outputs remain historical diagnostics | `results/coupled_equal_thickness_beta15_eps0p01_eb_timoshenko_3d_fem_vs_mu_corrected.png`, `results/coupled_equal_thickness_beta15_eps0p01_eb_timoshenko_3d_fem_vs_mu_corrected.csv`, `results/coupled_equal_thickness_beta15_eps0p01_eb_timoshenko_3d_fem_vs_mu_corrected_report.md`, `results/solid_fem_coupled_equal_rods_beta15_eps0p01_planar_mu/` | Corrected diagnostic-only `beta=15 deg`, `epsilon=0.01`, `eta=0` plot over `mu=0..0.9`; `mu` changes rod lengths, so this is equal-thickness rather than equal-length except at `mu=0`. Analytic EB and Timoshenko curves are descendant branches; sparse FEM points use planar point-joint 3D solid FEM and MAC-like branch assignment. EB thin-rod diameter validity and Timoshenko cut-off validity are reported separately; do not dash Timoshenko curves using the EB diameter criterion. |
| Equal-thickness full rigid-joint `Lambda(mu)` epsilon sweep | `python scripts/analysis/plot_rigid_joint_equal_thickness_beta15_lambda_mu_eps_sweep.py` | does not replace the earlier planar `epsilon=0.01` overlay | `results/rigid_joint_equal_thickness_beta15_eps0p01_lambda_mu_eb_timoshenko_fem.png`, `results/rigid_joint_equal_thickness_beta15_eps0p025_lambda_mu_eb_timoshenko_fem.png`, `results/rigid_joint_equal_thickness_beta15_eps0p05_lambda_mu_eb_timoshenko_fem.png`, `results/rigid_joint_equal_thickness_beta15_lambda_mu_eb_timoshenko_fem.csv`, `results/rigid_joint_equal_thickness_beta15_lambda_mu_eb_timoshenko_fem_report.md`, `results/solid_fem_rigid_joint_equal_thickness_beta15_lambda_mu/` | Current main diagnostic visualization for comparing EB, Timoshenko, and full rigid end-face 3D FEM along `mu` at `beta=15 deg`, `eta=0`, and `epsilon=0.01,0.025,0.05`. Analytic curves are descendant branches over `mu=0..0.9`; FEM uses a 16-point `mu` grid with no planar constraint, no fused volume, and no patch tuning. Main plots include only strong/moderate non-ambiguous MAC rows; weak, ambiguous, and duplicate matches remain in the CSV/report. |
| Variable-length Timoshenko limit audit | `python scripts/analysis/audit_variable_length_timoshenko_limits.py` | reuses `scripts/lib/variable_length_timoshenko.py`, copied from the corrected `Lambda(mu)` diagnostic formulas | `results/variable_length_timoshenko_limits_audit.csv`, `results/variable_length_timoshenko_limits_audit.md` | Diagnostic-only sorted-root audit for the variable-length Timoshenko implementation. It checks `mu=0` consistency against the older equal-rods Timoshenko diagnostic, the `epsilon -> 0` EB limit for several beta/mu values, and `beta=0` straight-rod mu-invariance. Current run passes all three checks for `eta=0`; that audit does not verify `eta != 0` tau scaling, descendant branch identity, FEM validation, or high-frequency cut-off behavior. |
| Equal-thickness `mu=0` FEM sanity audit | `python scripts/analysis/audit_coupled_equal_thickness_beta15_eps0p01_mu0_fem_sanity.py` | consumes the corrected plot CSV and generated planar FEM case | `results/coupled_equal_thickness_beta15_eps0p01_mu0_fem_sanity_audit.md`, `results/coupled_equal_thickness_beta15_eps0p01_mu0_fem_sanity_audit.csv` | Diagnostic audit of the `beta=15 deg`, `epsilon=0.01`, `eta=0`, `mu=0` discrepancy. Records the current variable-length Timoshenko formulas, checks `Lambda_FEM=sqrt(Omega/epsilon)`, audits point-joint planar constraints/node sets, compares sorted and MAC-order planar FEM against EB/Timoshenko, and uses the existing 1D EB frame FEM as a sorted-frequency sanity anchor. |
| Point-joint 3D FEM discrepancy isolation audit | `python scripts/analysis/audit_3d_fem_point_joint_mu0_eps0p01.py` | reuses point-joint Gmsh/CalculiX helpers and writes only isolated audit outputs | `results/3d_fem_point_joint_mu0_eps0p01_audit.md`, `results/3d_fem_point_joint_mu0_eps0p01_audit.csv`, `results/solid_fem_point_joint_mu0_eps0p01_audit/` | Diagnostic-only isolation pass for `beta=15 deg`, `mu=0`, `eta=0`, `epsilon=0.01`. It compares the current planar-constrained point-joint, the full point-joint without planar constraints, and a straight two-half-rod point-joint sanity case. Current result: the discrepancy is already present in full point-joint, while the straight joint sanity case has small strong-MAC errors; this points away from the planar constraint as the sole source and away from a generic rigid-reference-joint failure. |
| Point-joint coupling-model audit | `python scripts/analysis/audit_3d_fem_joint_coupling_models.py` | writes isolated Gmsh/CalculiX cases only under the audit output directory | `results/3d_fem_joint_coupling_models_audit.md`, `results/3d_fem_joint_coupling_models_audit.csv`, `results/solid_fem_joint_coupling_models_audit/` | Diagnostic-only comparison of alternative joint coupling idealizations at `epsilon=0.01`, `mu=0`, `eta=0`, and `beta=0,15,45,90 deg`. It tests full rigid end faces and central rigid patches with radii `0.25*r` and `0.5*r`; the `0.5*r` patch is the current best diagnostic candidate because it passes the beta=0 sanity check and reduces the beta=15 error, but no model is article-ready. |
| Point-joint constraint bracketing audit | `python scripts/analysis/audit_3d_fem_joint_constraint_bracketing.py` | writes isolated Gmsh/CalculiX cases only under the bracketing audit output directory | `results/3d_fem_joint_constraint_bracketing_audit.md`, `results/3d_fem_joint_constraint_bracketing_audit.csv`, `results/solid_fem_joint_constraint_bracketing_audit/` | Diagnostic-only bracketing study at `epsilon=0.01`, `mu=0`, `eta=0`, and `beta=0,15,45,90 deg`. It keeps the rigid-end-face and center-patch probes, adds translation-only central/average equation probes, a rotation-clamped reference-node bracket, and a two-reference-node constraint attempt. Current result: translation-only probes are too soft and fail beta=0 sanity, the rotation clamp reproduces rigid-end-face errors while remaining physically incompatible with the analytic free point-joint rotation, the two-reference-node deck is unsupported by the attempted CalculiX constraint combination, and `center_patch_0p5` remains diagnostic promising but not article-ready or patch-radius calibration. |
| Rigid end-face thickness trend audit | `python scripts/analysis/audit_3d_fem_rigid_joint_thickness_trend.py` | writes isolated full rigid end-face Gmsh/CalculiX cases only under the trend output directory | `results/3d_fem_rigid_joint_thickness_trend_report.md`, `results/3d_fem_rigid_joint_thickness_trend.csv`, `results/3d_fem_rigid_joint_thickness_trend_error_ratio.png`, `results/solid_fem_rigid_joint_thickness_trend/` | Diagnostic-only comparative trend study for the physically clear full rigid end-face point-joint at `beta=15 deg`, `mu=0`, `eta=0`, and `epsilon=0.005,0.01,0.025,0.05`. It uses no planar constraint, no fused volume, and no patch tuning; MAC-matched strong/moderate rows show the mean Timoshenko/EB error ratio decreasing from about `1.005` at `epsilon=0.005` to about `0.621` at `epsilon=0.05`, supporting the comparative thickness trend while leaving mesh convergence and visual mode review pending. |
| Rigid end-face equal-thickness frequency trend plot | `python scripts/analysis/plot_rigid_joint_equal_thickness_beta15_eps_trend.py` | consumes the existing rigid end-face trend CSV; does not rerun Gmsh/CalculiX | `results/rigid_joint_equal_thickness_beta15_frequency_comparison_eps_trend.png`, `results/rigid_joint_equal_thickness_beta15_frequency_comparison_eps_trend.csv`, `results/rigid_joint_equal_thickness_beta15_frequency_comparison_eps_trend_report.md` | Diagnostic-only visualization of the `beta=15 deg`, `mu=0`, `eta=0` rigid-joint thickness benchmark for `epsilon=0.01,0.025,0.05`. It plots EB, Timoshenko, and MAC-eligible full rigid end-face 3D FEM frequencies for the first six branches, records weak/ambiguous exclusions, and keeps EB diameter applicability separate from the Timoshenko cut-off without vertical applicability lines or Timoshenko style changes based on the EB criterion. |
| Coupled equal rods beta=15 3D solid FEM workflow | `python scripts/analysis/solid_fem_coupled_equal_rods_beta15.py` | optional `GMSH_EXE`, `CCX_EXE` | `results/coupled_equal_rods_beta15_3d_solid_fem_report.md`, `results/coupled_equal_rods_beta15_3d_solid_fem_sorted_comparison.csv`, `results/coupled_equal_rods_beta15_3d_solid_fem_mode_metrics.csv`, `results/coupled_equal_rods_beta15_3d_solid_fem_in_plane_comparison.csv`, `results/solid_fem_coupled_equal_rods_beta15/` | Diagnostic-only 3D solid check for beta=15, mu=0, eta=0 equal rods. Generates a fused two-cylinder mesh, coordinate-derived outer fixed node sets, CalculiX modal runs, `.frd` shape classification, and classified in-plane FEM vs EB/Timoshenko sorted-frequency comparisons. |
| Coupled equal rods EB/Timoshenko check | `python scripts/analysis/compare_coupled_equal_rods_eb_timoshenko.py` | none | `results/coupled_equal_rods_beta15_eb_timoshenko_frequencies.csv`, `results/coupled_equal_rods_beta15_timoshenko_kappa_sensitivity.csv`, `results/coupled_equal_rods_beta15_eb_timoshenko_lambda_vs_epsilon.png`, `results/coupled_equal_rods_beta15_eb_timoshenko_report.md` | Preferred diagnostic for the first coupled-rod Timoshenko sanity check at `beta=15 deg`, `mu=0`, `eta=0`; compares sorted frequencies only, marks `4*epsilon <= 0.1`, records cut-off margins, and runs a small kappa sensitivity check. |
| Variable-length thickness-mismatch Timoshenko tau audit | `python scripts/analysis/audit_variable_length_thickness_timoshenko_limits.py` | reuses `scripts/lib/variable_length_timoshenko.py` and the EB eta determinant in `src/my_project/analytic/formulas_thickness_mismatch.py` | `results/variable_length_thickness_timoshenko_limits_audit.csv`, `results/variable_length_thickness_timoshenko_limits_audit.md` | Diagnostic-only sorted-root audit for tau-aware eta support. It checks eta=0 regression against the legacy variable-length Timoshenko path, eta-to-zero continuity, `(mu, eta) -> (-mu, -eta)` swap symmetry, the `epsilon -> 0` EB thickness-mismatch limit, and beta=0 straight composite-rod behavior. Current run passes these gates for the tested first six sorted roots; descendant branch identity and energy plots remain separate. |
| Timoshenko energy partition audit | `python scripts/analysis/audit_timoshenko_energy_partition_beta15_eps0p01_eta0_eta0p5.py` | reuses `scripts/lib/variable_length_timoshenko.py` reconstruction and energy helpers | `results/timoshenko_energy_partition_beta15_eps0p01_eta0_eta0p5.csv`, `results/timoshenko_energy_partition_beta15_eps0p01_eta0_eta0p5_report.md`, diagnostic PNGs in `results/` | Diagnostic-only descendant audit at `beta=15 deg`, `epsilon=0.01`, `eta=0,0.5`, and `mu=0..0.9`. It partitions Timoshenko potential energy into bending/shear/axial terms per rod, records EB diameter and Timoshenko cutoff margins, and compares matching EB descendants. Use it only to explain why EB/Timoshenko frequencies can remain close when the EB diameter criterion fails locally; it does not change article figures, FEM validation, old determinants, old solvers, or baseline results. |
| Eta-zero and swap checks | `python scripts/analysis/check_thickness_mismatch_eta_zero_limit.py` | none | `results/thickness_mismatch_eta_zero_roots_check.csv`, `results/thickness_mismatch_swap_symmetry_check.csv` | Sanity checks for the diagnostic eta extension. |
| Descendant `Lambda(mu)` eta sweep | `python scripts/analysis/plot_lambda_mu_thickness_mismatch_beta15_eta_descendants.py --betas 7.5` | `plot_lambda_mu_thickness_mismatch_beta15_eta_large_slenderness.py`, `track_lambda_mu_thickness_mismatch_eta_sweep.py` | `results/thickness_mismatch_lambda_mu_beta{beta}_eps0p0025_eta_m0p5_0_p0p5_descendants.*` | Current clean presentation-style plot; no args keeps the historical beta=15 output path, `--betas` writes one PNG/report/tracking-CSV/warning-CSV set per angle, and tracking warnings stay in reports/data, not as x-markers. |
| Eta=0.5 global spectrum overview | `python scripts/analysis/plot_thickness_mismatch_eta_p0p5_global_spectrum.py` | none | `results/thickness_mismatch_eta_p0p5_beta15_global_spectrum_8modes.*` | Shows sorted roots, descendant branches, canonical sorted positions, and unresolved candidate assignments. |
| Lambda(beta) sorted/descendant identity audit | `python scripts/analysis/audit_lambda_beta_sorted_descendant_thickness_mismatch.py` | `formulas_thickness_mismatch.py`, `thickness_mismatch_mac_tracking.py` | `results/eta0p5_eps0p0025_mu0_lambda_beta_sorted_descendant_audit.*`, `results/eta0p5_eps0p0025_mu0_sorted5_identity_summary.csv`, plus eta=0 comparison outputs | Analytic-only EB beta sweep at `epsilon=0.0025`, `mu=0`, and `eta=0.5` with an eta=0 comparison. Descendants are seeded at `beta=0`, `mu=0`; sorted position remains metadata. The current audit finds sorted 5 at `beta=15 deg`, `eta=0.5` is descendant 5, while the eta=0 comparison has sorted 5 as descendant 6. |
| Eta scan first-six rearrangement audit | `python scripts/analysis/audit_eta_scan_first6_rearrangement.py` | `audit_lambda_beta_sorted_descendant_thickness_mismatch.py`, `formulas_thickness_mismatch.py`, `thickness_mismatch_mac_tracking.py` | `results/eta_scan_eps0p0025_mu0_first6_rearrangement.csv`, `results/eta_scan_eps0p0025_mu0_first6_rearrangement_summary.csv`, `results/eta_scan_eps0p0025_mu0_first6_rearrangement_report.md`, `results/eta_scan_eps0p0025_mu0_first6_*.png` | Analytic-only EB eta audit over `eta=-0.5..0.5`, `beta=0..90 deg`, `epsilon=0.0025`, `mu=0`. Rearrangement means any descendant among 1..6 changes sorted position over beta. The current grid finds no rearrangement for `|eta| >= 0.018`, clean rearrangement near zero for `0.002 <= |eta| <= 0.016`, and the exact `eta=0` equal-thickness case as tracking-unreliable with pairwise sorted-position exchanges. |
| Eta=0 mu scan first-six rearrangement audit | `python scripts/analysis/audit_mu_scan_eta0_first6_rearrangement.py` | `audit_lambda_beta_sorted_descendant_thickness_mismatch.py`, `formulas_thickness_mismatch.py`, `thickness_mismatch_mac_tracking.py` | `results/mu_scan_eta0_eps0p0025_first6_rearrangement.csv`, `results/mu_scan_eta0_eps0p0025_first6_rearrangement_summary.csv`, `results/mu_scan_eta0_eps0p0025_first6_rearrangement_report.md`, `results/mu_scan_eta0_eps0p0025_*.png` | Analytic-only EB mu audit over `mu=0..0.9`, `eta=0`, `beta=0..90 deg`, and `epsilon=0.0025`. Rearrangement means any descendant among 1..6 changes sorted position over beta after beta=0 seed degeneracy handling. The current grid finds true rearrangement only for `0.001 <= mu <= 0.002`, marks `mu=0` tracking-unreliable, and finds no rearrangement for `0.003 <= mu <= 0.9`. |
| Eta=0 fixed-mu `Lambda(beta)` rearrangement check | `python scripts/analysis/plot_lambda_beta_eta0_eps0p0025_mu_rearrangement_check.py` | `audit_mu_scan_eta0_first6_rearrangement.py`, `formulas_thickness_mismatch.py`, `thickness_mismatch_mac_tracking.py` | `results/lambda_beta_eta0_eps0p0025_mu0p001_first6.png`, `results/lambda_beta_eta0_eps0p0025_mu0p002_first6.png`, `results/lambda_beta_eta0_eps0p0025_mu0p003_first6.png`, `results/lambda_beta_eta0_eps0p0025_mu_rearrangement_check.csv`, `results/lambda_beta_eta0_eps0p0025_mu_rearrangement_check_report.md` | Analytic-only EB refined beta-grid plot for `mu=0.001,0.002,0.003`, `eta=0`, and `epsilon=0.0025`. It uses the mu-scan shape-MAC descendant tracker, beta step `0.1 deg`, warning markers, and spike checks. Current result: no reproducible first-six rearrangement, no real crossings, and no spike artifacts for these three fixed-mu checks; the earlier `0.25 deg` narrow mu-scan signal is grid/tracking-sensitive diagnostic evidence. |
| Eta=0 positive-mu adjacent-gap verification | `python scripts/analysis/audit_eta0_mu_positive_gap_rearrangement_verification.py` | `audit_mu_scan_eta0_first6_rearrangement.py`, `formulas_thickness_mismatch.py`, `thickness_mismatch_mac_tracking.py` | `results/eta0_eps0p0025_mu_positive_gap_verification_summary.csv`, `results/eta0_eps0p0025_mu_positive_gap_verification_details.csv`, `results/eta0_eps0p0025_mu_positive_gap_verification_report.md`, `results/eta0_eps0p0025_mu_positive_gap_scaling.png`, `results/eta0_eps0p0025_mu_positive_gap_min_beta.png`, `results/eta0_eps0p0025_mu_positive_gap_lambda_beta_selected.png` | Strict analytic-only EB audit for positive `mu=1e-6..0.9` at `eta=0`, separating adjacent sorted-root gaps from descendant-label tracking. It refines beta-local gap minima for sorted pairs 1-2 through 6-7 and compares default/strict root settings. Current result: all 19 tested positive `mu` rows have resolved positive sorted gaps; the smallest is about `1.314e-4` at `mu=1e-6`, `beta=13.308 deg`, pair 5-6. Small-mu descendant swaps are fragile tracking diagnostics, not accepted true crossings. |
| Eta-parameter positive-gap verification | `python scripts/analysis/audit_eta_parameter_positive_gap_verification.py` | `formulas_thickness_mismatch.py`, `thickness_mismatch_mac_tracking.py` | `results/eta_parameter_positive_gap_verification_summary.csv`, `results/eta_parameter_positive_gap_verification_details.csv`, `results/eta_parameter_positive_gap_verification_report.md`, `results/eta_parameter_positive_gap_min_gap_vs_eta.png`, `results/eta_parameter_positive_gap_min_gap_vs_beta_eta_scan.png`, `results/eta_parameter_positive_gap_min_gap_vs_mu_eta_scan.png`, `results/eta_parameter_positive_gap_classification_map.png`, plus selected `Lambda` plots | Strict analytic-only EB audit for eta-, beta-, and mu-parameter positive-gap checks at `epsilon=0.0025`. It solves the first 14 sorted roots, audits adjacent sorted pairs 1-2 through 6-7, tracks first-eight descendants only as label diagnostics, applies local strict root repairs where default sign scans miss close roots, and refines the special `mu` windows near `0.15864`, `0.379`, and `0.716`. Current result: all tested rows are resolved positive gaps with no true sorted-root crossings or unresolved cases; the smallest tested eta-nonzero gap is about `4.0747e-4` for `beta=5 deg`, `eta=0.5`, sorted pair 1-2. |
| Eta-mu-beta sorted-frequency maps | `python scripts/analysis/plot_diagnostic_eta_mu_beta_frequency_maps.py` | `formulas_thickness_mismatch.py` | `results/diagnostic_lambda_eta_beta0_mu0_eps0p0025.*`, `results/diagnostic_lambda_beta_eps0p0025_mu_eta_slices.csv`, `results/diagnostic_lambda_beta_eps0p0025_mu_eta_grid_overview.png`, `results/diagnostic_eta_mu_beta_heatmap_metrics_eps0p0025.*`, and six `results/diagnostic_heatmap_*_eta_mu_eps0p0025.png` heatmaps | Diagnostic-only sorted-frequency mapping for the eta/mu/beta landscape. It plots `Lambda(eta)` at `beta=0, mu=0`, sorted `Lambda(beta)` slices for `mu=0.1,0.2,0.3` and `eta=+-0.1,+-0.2,+-0.3`, and eta-mu heatmaps of minimum adjacent sorted gaps plus beta-sensitivity metrics. Current fallback heatmap grid is eta step `0.025`, mu step `0.05`, beta step `1 deg`, first 10 sorted roots; current run finds min `g_min=0.00838632600513` at `eta=0.5`, `mu=0.35`, `beta=23 deg`, pair 6-7, and the max `S_mean`/`S_max` location at `eta=-0.3`, `mu=0.6`, with no solver warnings or missing roots. These are diagnostic close-approach candidates, not crossing/no-crossing claims. |
| Eta=0.1, mu=0.3 branch-2 beta transition audit | `python scripts/analysis/audit_mode_shape_branch2_beta_transition_eta0p1_mu0p3.py` | `plot_mode_shapes_eta_beta_scan.py`, `formulas_thickness_mismatch.py`, `thickness_mismatch_mac_tracking.py` | `results/mode_shape_branch2_transition_audit_eta_0p1_mu_0p3_eps_0p0025.csv`, `results/mode_shape_branch2_transition_audit_eta_0p1_mu_0p3_eps_0p0025.md`, `results/mode_shape_branch2_transition_*_eta_0p1_mu_0p3_eps_0p0025.png` | Local diagnostic-only audit for the apparent descendant-2 sorted-index change in the eta=0.1, mu=0.3 mode-shape beta scan. It checks the sorted 2-3 gap over beta=4..7 deg at step 0.05, tracks descendants 2 and 3 only as shape-label diagnostics from beta=0, and reports a positive sorted gap with local tracking ambiguity rather than a frequency crossing. |
| Eta=0.1, mu=0.3 fixed sorted 2/3 beta-scan mode shapes | `python scripts/analysis/plot_mode_shapes_eta_beta_scan_sorted_modes.py` | `plot_mode_shapes_eta_beta_scan.py`, `formulas_thickness_mismatch.py`, `thickness_mismatch_mac_tracking.py` | `results/mode_shapes_sorted2_sorted3_eta_0p1_mu_0p3_eps_0p0025/` | Diagnostic-only analytic EB mode-shape plots for fixed sorted indices 2 and 3 over beta=0..15 deg. Selection is the sorted root position independently at each beta, not descendant tracking; the default `short-rod-up` sign convention orients plotted eigenvectors, while MAC and sorted 2-3 gap outputs are metadata for interpretation. |
| Eta=+-0.5 sorted `Lambda(mu)` with single-rod references | `python scripts/analysis/plot_lambda_mu_eta_m0p5_with_single_beam_refs.py --eta 0.5 --betas 0 5 10 15` | `formulas_thickness_mismatch.py`, `FreqMuNet.py`, `solvers.py` | `results/diagnostic_lambda_mu_eta_*eps0p0025_beta_*deg_single_beam_refs.*` | Diagnostic-only sorted-root plots for selected beta values at `epsilon=0.0025`. The coupled curves are sorted frequencies, not descendants. CP/FP and CC/FF single-rod reference families are variable `Lambda(mu)` curves computed separately as `alpha*sqrt(tau_i)/L_i`; they are not horizontal reference lines. For `eta=-0.5`, rod 1 is short/thick and rod 2 is long/thin; for `eta=0.5`, rod 1 is short/thin and rod 2 is long/thick, so reference styles follow the actual thin/thick role rather than rod number. Default PNGs zoom to the coupled spectrum and clip reference segments only for plotting; `--full-scale-refs` writes separate full-scale PNGs. |
| Eta-mu-beta global trend post-processing | `python scripts/analysis/summarize_diagnostic_eta_mu_beta_global_trends.py` | consumes the generated diagnostic map CSVs only | `results/diagnostic_eta_mu_beta_global_trends_by_mu_eps0p0025.csv`, `results/diagnostic_eta_mu_beta_global_trends_by_eta_eps0p0025.csv`, `results/diagnostic_eta_mu_beta_global_trends_eta_sign_abs_eps0p0025.csv`, `results/diagnostic_eta_mu_beta_pair_transition_by_mu_eps0p0025.csv`, `results/diagnostic_eta_mu_beta_beta_at_gmin_histogram_eps0p0025.csv`, `results/diagnostic_eta_mu_beta_sensitivity_branch_dominance_eps0p0025.csv`, `results/diagnostic_eta_mu_beta_global_surrogate_fits_eps0p0025.csv`, `results/diagnostic_eta_mu_beta_global_trends_eps0p0025_report.md`, and eight `results/diagnostic_global_*_eps0p0025.png` plots | CSV-only global post-processing of the sorted-frequency map. It computes by-mu/by-eta quantiles, eta-sign asymmetry, pair transitions, beta-at-minimum histograms, branch-sensitivity dominance, and low-order surrogate fits without recomputing roots or doing local refinement. Current result: mu is the stronger global driver of `g_min`, `S_mean`, and `S_mean_rel`; median `g_min` increases with mu, while median beta sensitivity decreases with mu. |
| Selected-eta `Lambda(beta)` descendant plots | `python scripts/analysis/plot_eta_scan_selected_lambda_beta_descendants.py` | `audit_lambda_beta_sorted_descendant_thickness_mismatch.py`, `formulas_thickness_mismatch.py`, `thickness_mismatch_mac_tracking.py` | `results/eta_scan_selected_lambda_beta_descendants_eps0p0025_mu0.png`, `results/eta_scan_selected_lambda_beta_descendants_eps0p0025_mu0.csv`, `results/eta_scan_selected_lambda_beta_descendants_eps0p0025_mu0_report.md` | Analytic-only 2x3 plot of smooth tracked descendant `Lambda(beta)` curves for eta values `-0.016,-0.004,-0.002,0.002,0.004,0.016` at `epsilon=0.0025`, `mu=0`. Uses shape-MAC continuation from `beta=0`, accepts near-degenerate best-MAC steps as diagnostic metadata instead of falling back to sorted roots, and runs a spike-artifact check on the plotted curves. |
| Refined branch-identity audit | `python scripts/analysis/check_thickness_mismatch_branch_identity_eta_p0p5.py` | none | `results/thickness_mismatch_branch_identity_eta_p0p5_beta15_refined.*` | Local descendant audit for branches 5--7 at eta=0.5. |
| Eta=0.5 crossing audit | `python scripts/analysis/audit_eta0p5_eps0p0025_lambda_mu_crossings.py` | `formulas_thickness_mismatch.py`, `thickness_mismatch_mac_tracking.py` | `results/eta0p5_eps0p0025_beta15_lambda_mu_crossing_audit.*` | Pairwise diagnostic for real descendant-frequency crossings versus finite gaps, sorted-position metadata, and tracking artifacts on the dense `mu=0..0.9` grid. |
| Article beta=5 desc1/desc2 crossing audit | `python scripts/analysis/audit_article_beta5_eta0p5_desc1_desc2_crossing.py` | `formulas_thickness_mismatch.py`, `thickness_mismatch_mac_tracking.py`, article local CSV metadata | `results/article_check_beta5_eta0p5_eps0p0025_mu_desc1_desc2_crossing.*` | Focused analytic-only check for the article `beta=5 deg`, `epsilon=0.0025`, `eta=0.5` `Lambda(mu)` figure. It resolves the apparent desc1/desc2 crossing with a local `mu` refinement and fine root scan, reporting that the refined descendants keep a finite positive gap and no sorted-position swap. |
| Descendant-5 corrected veering shape points | `python scripts/analysis/plot_thickness_mismatch_desc5_veering_mode_shapes.py` | `plot_thickness_mismatch_branch_shapes_vs_eta.py`, existing descendant-6 veering points table | `results/thickness_mismatch_desc5_veering_mode_shapes_beta15_eps0p0025_eta0p5/desc5_near45_rank*.png`, `results/thickness_mismatch_desc5_veering_mode_shapes_beta15_eps0p0025_eta0p5/desc5_near56_rank*.png`, `results/thickness_mismatch_desc5_veering_mode_shapes_beta15_eps0p0025_eta0p5/desc5_veering_points.*` | Analytic-only corrected shape set for the eta=0.5 sorted/descendant identity audit: reuses the six mu/local-gap labels from the descendant-6 table but reconstructs descendant 5, because sorted 5 is descendant 5 for eta=0.5. Descendant-6 outputs are retained separately. |
| Descendant-6 veering shape points | `python scripts/analysis/plot_thickness_mismatch_desc6_veering_mode_shapes.py` | `plot_thickness_mismatch_branch_shapes_vs_eta.py`, existing crossing audit outputs | `results/thickness_mismatch_desc6_veering_mode_shapes_beta15_eps0p0025_eta0p5/desc6_near45_rank*.png`, `results/thickness_mismatch_desc6_veering_mode_shapes_beta15_eps0p0025_eta0p5/desc6_near56_rank*.png`, `results/thickness_mismatch_desc6_veering_mode_shapes_beta15_eps0p0025_eta0p5/desc6_veering_points.*` | Analytic-only TASK 2 wrapper: recomputes dense sorted adjacent gaps for eta=0.5, selects three distinct local minima each for sorted pairs 4-5 and 5-6, then plots descendant 6 shapes with ranked CSV/report metadata. |
| Descendant-6 veering localization audit | `python scripts/analysis/audit_desc6_veering_localization.py` | `plot_thickness_mismatch_branch_shapes_vs_eta.py`, `desc6_veering_points.csv` | `results/thickness_mismatch_desc6_veering_mode_shapes_beta15_eps0p0025_eta0p5/desc6_localization_audit.csv`, `results/thickness_mismatch_desc6_veering_mode_shapes_beta15_eps0p0025_eta0p5/desc6_localization_audit_report.md` | Analytic-only rod-localization audit for the six selected descendant-6 veering shapes; reports displacement, slope, curvature, and tau-weighted EB bending-energy fractions without selecting new mu points. |
| Eta=0.5 analytic/FEM overlay | `python scripts/analysis/plot_thickness_mismatch_fem_comparison_eta_p0p5_beta15.py` | `fem_check_thickness_mismatch_eta_p0p5_beta15.py` | `results/thickness_mismatch_fem_comparison_beta15_eps0p0025_eta_p0p5.*` | Lightweight plotting comparison. The older FEM check keeps the detailed CSV/MAC diagnostics. |
| Eta=0.5 isolated-rod references, beta=15 | `python scripts/analysis/plot_lambda_mu_thickness_mismatch_eta_p0p5_with_isolated_rods.py` | none | `results/thickness_mismatch_lambda_mu_beta15_eps0p0025_eta_p0p5_with_isolated_rods.*` | Uses the clamped-supported / clamped-pinned reference convention and eta-normalized `Lambda_rod`. |
| Eta=0.5 isolated-rod references, beta=45 | `python scripts/analysis/plot_lambda_mu_thickness_mismatch_beta45_eta_p0p5_with_isolated_rods.py` | none | `results/thickness_mismatch_lambda_mu_beta45_eps0p0025_eta_p0p5_with_isolated_rods.*` | Same reference convention as beta=15, different beta only. |
| Eta=0.5 fixed-fixed isolated-rod references, beta=45 | `python scripts/analysis/plot_lambda_mu_thickness_mismatch_beta45_eta_p0p5_with_isolated_rods_fixed_fixed.py` | none | `results/thickness_mismatch_lambda_mu_beta45_eps0p0025_eta_p0p5_with_isolated_rods_fixed_fixed.*` | Uses clamped-clamped / fixed-fixed reference roots. |
| Branch-shape eta overlay | `python scripts/analysis/plot_thickness_mismatch_branch_shapes_vs_eta.py --branch-index 5 --mus 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 --etas -0.5 0 0.5` | `plot_thickness_mismatch_branch5_shapes.py` | `results/thickness_mismatch_desc5_mode_shapes_beta15_eps0p0025/`, `results/thickness_mismatch_desc6_veering_mode_shapes_beta15_eps0p0025_eta0p5/` | Parameterized descendant-shape overlay for several eta values at target mu values, using shape-MAC tracking from `mu=0`; optional compact filenames and aggregate CSV/report outputs support descendant-5 mu sweeps and descendant-6 near-interaction shape diagnostics. |
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

For `Lambda(mu)` diagnostics, keep applicability criteria separated:
Euler--Bernoulli diameter validity uses `2*r_i/l_i <= 0.1`, while
Timoshenko validity is audited with `Omega/Omega_c`. Do not dash or exclude
Timoshenko curves solely because the Euler--Bernoulli diameter criterion fails.

The variable-length Timoshenko implementation is considered
diagnostic-verified only when the limit audit passes all three gates:
`mu=0` consistency with the older equal-rods Timoshenko diagnostic, the
`epsilon -> 0` Euler--Bernoulli limit, and `beta=0` straight-rod
mu-invariance. The current audit output passes these gates for `eta=0` sorted
roots.

The tau-aware `eta != 0` extension in `scripts/lib/variable_length_timoshenko.py`
is considered diagnostic-verified only for the sorted-root checks in
`results/variable_length_thickness_timoshenko_limits_audit.md`: eta=0
regression, eta-to-zero continuity, rod-swap symmetry, the `epsilon -> 0`
Euler--Bernoulli thickness-mismatch limit, and beta=0 straight composite-rod
behavior. Descendant branch identity, energy plots, and FEM validation remain
separate.

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
- `scripts/lib/variable_length_timoshenko.py` contains the diagnostic
  variable-length Timoshenko matrix helpers used by the eta=0 and tau-aware
  eta audits; it does not replace the old Euler--Bernoulli determinant or old
  solvers.

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
