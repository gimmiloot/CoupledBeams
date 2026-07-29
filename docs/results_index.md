# Generated Results Index

`results/` is generated and ignored by Git, apart from its tracked placeholder.
A fresh clone therefore does not contain the output directories listed below.
An absent result path is not automatically a broken documentation reference:
it may be an expected output of a tracked workflow.

Canonical conclusions must be recorded in tracked documentation. Generated
reports, CSV files, plots, solver caches, and external-program artifacts remain
local evidence. Reproduction commands and their assumptions are documented in
the [scripts guide](../scripts/README.md), [workflow status
map](../scripts/STATUS.md), and [thickness-mismatch script
map](../scripts/analysis/thickness_mismatch/README.md). No command in this
index should be run without first reviewing its cost and output contract.

## Result directories

All 23 immediate ignored result directories present in the pre-refactor
snapshot are included here.

| Result directory | Scientific workflow | Status | Canonical report or documentation | Reproduction entry point | Local archive class |
| --- | --- | --- | --- | --- | --- |
| `results/_smoke/` | Small wiring/environment checks from several workflows | `temporary` | Generated reports below `_smoke/`; [smoke convention](../scripts/analysis/thickness_mismatch/README.md#smoke-mode-convention) | documented `--smoke` modes; 3D environment check: `scripts/analysis/thickness_mismatch/audits/check_3d_fem_environment.py` | `temporary-regenerable`; do not delete outside a dedicated cleanup task |
| `results/eb_epsilon_apriori_pilot/` | Original 21-case geometry-only epsilon pilot | `completed-diagnostic` | `results/eb_epsilon_apriori_pilot/analysis/epsilon_apriori_pilot_report.md`; [research plan](thickness_mismatch/eb_safe_spectrum_prefix_research_plan.md) | `scripts/analysis/thickness_mismatch/audits/run_eb_epsilon_apriori_pilot.py` and the CSV postprocessor | `soft-archive-candidate`; preserve provenance |
| `results/eb_epsilon_apriori_pilot_branch_continuation_v1/` | Corrected branch-informed pilot | `completed-diagnostic` | generated `analysis/epsilon_apriori_pilot_report.md`; [research plan](thickness_mismatch/eb_safe_spectrum_prefix_research_plan.md) | `scripts/analysis/thickness_mismatch/audits/audit_eb_timo_branch_continuation_gateway.py --run-pilot` | `preserve-prerequisite` |
| `results/eb_epsilon_apriori_pilot_complete_spectrum_v1/` | Auto-complete-spectrum pilot used by the general audit | `completed-diagnostic` | generated `analysis/epsilon_apriori_pilot_report.md`; [research plan](thickness_mismatch/eb_safe_spectrum_prefix_research_plan.md) | `scripts/analysis/thickness_mismatch/audits/audit_eb_timo_general_spectrum_completeness.py` | `soft-archive-candidate`; superseded for the targeted gateway |
| `results/eb_epsilon_baseline_thresholds/` | Corrected factorized straight-system epsilon thresholds plus preserved legacy cache | `completed-diagnostic` | `results/eb_epsilon_baseline_thresholds/eb_epsilon_baseline_thresholds_report.md`; [research plan](thickness_mismatch/eb_safe_spectrum_prefix_research_plan.md) | `scripts/analysis/thickness_mismatch/audits/audit_eb_epsilon_baseline_thresholds.py` | `preserve-prerequisite`; separate corrected and legacy cache identities |
| `results/eb_epsilon_lower_envelope_step3a/` | Step-3A lower-envelope screen and counterexamples | `closed-research` | `results/eb_epsilon_lower_envelope_step3a/eb_epsilon_lower_envelope_step3a_report.md`; [stage closure](thickness_mismatch/eb_safe_prefix_stage_closure.md) | `scripts/analysis/thickness_mismatch/audits/audit_eb_epsilon_lower_envelope_step3a.py` | `preserve-canonical-evidence` |
| `results/eb_rule_ab_exact_pareto/` | Exact Rules A/B/S search, partitions, predictions, and audits | `closed-research` | `results/eb_rule_ab_exact_pareto/eb_rule_ab_exact_pareto_report.md`; [stage closure](thickness_mismatch/eb_safe_prefix_stage_closure.md) | `scripts/analysis/thickness_mismatch/postprocess/analyze_eb_rule_ab_exact_pareto.py` | `preserve-canonical-evidence` |
| `results/eb_rule_s_cost_break_even/` | Frozen five-case Rule-S engineering cost benchmark | `closed-research` | `results/eb_rule_s_cost_break_even/rule_S_cost_break_even_report.md`; [stage closure](thickness_mismatch/eb_safe_prefix_stage_closure.md) | `scripts/analysis/thickness_mismatch/benchmarks/benchmark_rule_s_cost_break_even.py` | `preserve-canonical-negative-result` |
| `results/eb_timo_branch_continuation_gateway/` | Branch-informed K10/root-11 readiness gateway | `completed-diagnostic` | `results/eb_timo_branch_continuation_gateway/eb_timo_branch_continuation_gateway_report.md`; [research plan](thickness_mismatch/eb_safe_spectrum_prefix_research_plan.md) | `scripts/analysis/thickness_mismatch/audits/audit_eb_timo_branch_continuation_gateway.py` | `preserve-prerequisite` |
| `results/eb_timo_clean_mode_shapes/` | Clean corrected EB/Timoshenko full-shape grids | `active-diagnostic` | `results/eb_timo_clean_mode_shapes/clean_mode_shape_report.md`; [script map](../scripts/analysis/thickness_mismatch/README.md) | `scripts/analysis/thickness_mismatch/shapes/plot_eb_timo_full_mode_shapes_eps0p03_beta45_eta0_modes4_6.py` | `active-local` |
| `results/eb_timo_counterexample_dimensional_frequency_beta/` | Certified dimensional-frequency beta plots for `S3_12` and `S3_14` | `active-diagnostic` | generated `counterexample_dimensional_frequency_beta_report.md`; [frequency-map policy](thickness_mismatch/frequency_plot_computation_policy.md#current-counterexample-figures) | `scripts/analysis/thickness_mismatch/maps/plot_counterexample_dimensional_frequency_beta.py` | `preserve-certified-figure-data` |
| `results/eb_timo_general_spectrum_completeness/` | General spectrum-completeness audit with negative readiness result | `completed-diagnostic` | `results/eb_timo_general_spectrum_completeness/eb_timo_general_spectrum_completeness_report.md`; [research plan](thickness_mismatch/eb_safe_spectrum_prefix_research_plan.md) | `scripts/analysis/thickness_mismatch/audits/audit_eb_timo_general_spectrum_completeness.py` | `soft-archive-candidate`; preserve negative provenance |
| `results/eb_timo_mode_shapes_eps0p03_beta45_eta0_modes4_6/` | Earlier EB/Timoshenko full-displacement shape set | `active-diagnostic` | `results/eb_timo_mode_shapes_eps0p03_beta45_eta0_modes4_6/mode_shape_full_displacement_report.md`; [script map](../scripts/analysis/thickness_mismatch/README.md) | `scripts/analysis/thickness_mismatch/shapes/plot_eb_timo_full_mode_shapes_eps0p03_beta45_eta0_modes4_6.py` | `manual-review`; compare with clean corrected outputs before classification |
| `results/eb_validity_fixed_epsilon_geometry_scan/` | Fixed-epsilon geometry applicability source study | `active-diagnostic` | generated `eb_validity_fixed_epsilon_geometry_scan_report.md`; [stage closure](thickness_mismatch/eb_safe_prefix_stage_closure.md) | `scripts/analysis/thickness_mismatch/audits/audit_eb_validity_fixed_epsilon_geometry_scan.py` | `soft-archive-candidate`; helpers remain reused |
| `results/eb_validity_vs_timoshenko_stage1/` | Stage-1 EB/Timoshenko applicability source study | `active-diagnostic` | generated `eb_validity_vs_timoshenko_stage1_report.md`; [stage closure](thickness_mismatch/eb_safe_prefix_stage_closure.md) | `scripts/analysis/thickness_mismatch/audits/audit_eb_validity_vs_timoshenko_stage1.py` | `soft-archive-candidate`; preserve source evidence |
| `results/eb_vs_timoshenko_3d_validation/` | Straight uniform/stepped and related independent 3D FEM comparisons | `active-diagnostic` | generated case reports; [FEM validation status](thickness_mismatch/fem_validation_status.md) | `validate_eb_timo_3d_beta0_stepped.py` and `validate_eb_timo_3d_beta0_uniform_eps0p05.py` | `active-local`; expensive/manual external-tool provenance |
| `results/eb_vs_timoshenko_lambda_beta_cases/` | Sorted in-plane EB/Timoshenko `Lambda(beta)` maps | `active-diagnostic` | generated `eb_vs_timo_lambda_beta_cases_report.md`; [script map](../scripts/analysis/thickness_mismatch/README.md) | `scripts/analysis/thickness_mismatch/maps/plot_eb_vs_timoshenko_lambda_beta_cases.py` | `active-local-cache` |
| `results/eb_vs_timoshenko_lambda_mu_beta45_eta0_eps_scan/` | Fixed-beta epsilon-family `Lambda(mu)` comparison | `active-diagnostic` | generated `eb_vs_timo_lambda_mu_beta45_eta0_eps_scan_report.md`; [script map](../scripts/analysis/thickness_mismatch/README.md) | `scripts/analysis/thickness_mismatch/maps/plot_eb_vs_timoshenko_lambda_mu_cases.py` with the documented beta/eta/epsilon arguments | `active-local-cache` |
| `results/eb_vs_timoshenko_lambda_mu_cases/` | General sorted EB/Timoshenko `Lambda(mu)` cases | `active-diagnostic` | generated `eb_vs_timo_lambda_mu_cases_report.md`; [script map](../scripts/analysis/thickness_mismatch/README.md) | `scripts/analysis/thickness_mismatch/maps/plot_eb_vs_timoshenko_lambda_mu_cases.py` | `active-local-cache` |
| `results/eb_vs_timoshenko_longitudinal_suspect_modes/` | Longitudinal-character and joint-continuity audit | `active-diagnostic` | no canonical Markdown report detected in the snapshot; [script map](../scripts/analysis/thickness_mismatch/README.md) | `scripts/analysis/thickness_mismatch/audits/audit_longitudinal_suspect_modes_eb_timo.py` | `manual-review` |
| `results/timoshenko_mode_shape_diagnostics/` | Corrected vector/component diagnostics for modes 4--6 | `active-diagnostic` | generated `timoshenko_modes_4_6_shape_diagnostics_report.md`; [script map](../scripts/analysis/thickness_mismatch/README.md) | `scripts/analysis/thickness_mismatch/audits/audit_timoshenko_modes_4_6_shape_diagnostics.py` | `active-local` |
| `results/timoshenko_shape_bug_audit/` | Thin-limit shape/display-transform audit | `active-diagnostic` | generated `timoshenko_shape_bug_audit_report.md`; [script map](../scripts/analysis/thickness_mismatch/README.md) | `scripts/analysis/thickness_mismatch/audits/audit_timoshenko_shape_bug_thin_limit.py` | `preserve-diagnostic-correction-evidence` |
| `results/timoshenko_shape_construction_audit/` | Shape-construction residual, visualization, and provenance audit | `active-diagnostic` | generated `timoshenko_modes456_visualization_report.md`; [script status](../scripts/STATUS.md#manual-review-candidates) | preferred: `scripts/analysis/thickness_mismatch/audits/audit_timoshenko_shape_construction.py` | `manual-review`; adjacent possible-orphan producer must be audited |

The snapshot also contains two ignored files directly under `results/`:
`variable_length_timoshenko_limits_audit.csv` and
`variable_length_timoshenko_limits_audit.md`. They belong to the tau-aware
Timoshenko limit-verification line documented in the
[thickness-mismatch model note](thickness_mismatch/README.md).

## Interpretation and preservation

- `active-diagnostic` directories are local working evidence and are not
  automatic archive candidates.
- `completed-diagnostic` directories preserve prerequisite and reproducibility
  data even when a newer decision stage is canonical.
- `closed-research` directories contain counterexamples, exact decisions, or
  negative engineering results and must be preserved as scientific history.
- `temporary` means regenerable intent, not permission to delete in this task.

The full per-file ignored-data manifest is intentionally kept only in the
external local backup. Public documentation records aggregate result paths and
statuses, not private paths. The [archive policy](archive_policy.md) defines
the gates required before any future local-results move or cleanup.
