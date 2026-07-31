# Script and Workflow Status

This file answers which workflow to use and how older commands should be
interpreted. Detailed parameters and outputs remain in the [scripts
guide](README.md) and the [thickness-mismatch script
map](analysis/thickness_mismatch/README.md).

## Status vocabulary

- `stable-baseline`: preferred maintained baseline.
- `active-diagnostic`: current diagnostic workflow; not automatically an
  article-level or universal conclusion.
- `completed`: the declared finite workflow has reached its result.
- `closed-research`: the question/path is closed under its recorded scope;
  retain the command for provenance.
- `historical`: preserved to reproduce an earlier evidence stage.
- `superseded`: a newer path/report is canonical, but the old command remains
  available.
- `compatibility-wrapper`: stable old/public command forwarding to an existing
  implementation.
- `manual-review`: usage or provenance is not clear enough for archive action.

## Preferred stable entry points

The table groups related commands; it is not a duplicate list of every Python
file.

| Workflow | Preferred entry point | Status | Successor or canonical report | Recommendation |
| --- | --- | --- | --- | --- |
| Baseline analytic sweeps | `scripts/run/run_beta_sweep_mu0_four_radii.py`; `run_mu_sweep_beta0_four_radii.py`; `run_mu_sweep_fixed_beta_four_radii.py`; `run_mu_sweep_four_betas_analytic.py`; `run_lambda_mu_fixed_beta_analytic.py` | `stable-baseline` | [scripts guide](README.md) | Use the documented wrapper matching the requested beta/mu sweep. Preserve the imported root-level implementations. |
| Baseline FEM | `src/my_project/fem/python_fem.py`; `scripts/run/run_branchwise_fem_audit.py` | `stable-baseline` | [project FEM rules](../docs/project_rules.md#fem-comparison-rules) | Keep the production local-to-global transform and assembly convention unchanged. |
| Branch tracking and analytic shapes | `scripts/lib/analytic_branch_tracking.py`; `scripts/run/run_tracked_bending_descendant_shape_ru.py` | `stable-baseline` | [branch rules](../docs/project_rules.md#branch-identity-and-mode-descendants) | Use `branch_id` for identity and record `current_sorted_index` separately. |
| Veering, modal exchange, and localization | `scripts/analyze_target_descendants_beta15_r5.py`; `scripts/analyze_flat_mu_bending_energy.py`; focused audits under `scripts/analysis/` | `active-diagnostic` | [strict assessment](../docs/veering/strict_veering_assessment.md), [slow-evolution assessment](../docs/veering/mu_slow_evolution_assessment.md) | Preserve completed evidence; make strict claims only with tracked branches, paired gaps, and shape/MAC evidence. |
| Thickness-mismatch limits and maps | limit audits plus parameterized map/shape commands listed in the [script map](analysis/thickness_mismatch/README.md#preferred-entry-points) | `active-diagnostic` | [model note](../docs/thickness_mismatch/README.md) | Prefer parameterized entry points and distinguish sorted spectra from descendants. |
| Yartsev Chapter-2 free-free monoclinic rod | `analysis/anisotropic_rods/reproduce_yartsev_fig_2_2.py` | `completed` source-reproduction diagnostic | [single-rod reproduction](../docs/anisotropic_rods/yartsev_ch2_single_rod_reproduction.md): `PASS_WITHIN_GRAPH_RESOLUTION` | Preserve the corrected/printed distinction. This completed gate is not a coupled-rod or production anisotropic workflow. |
| Yartsev Chapter-2 cantilever monoclinic rod | `analysis/anisotropic_rods/reproduce_yartsev_ch2_cantilever.py` | `completed` source-reproduction and boundary-condition diagnostic | [cantilever reproduction](../docs/anisotropic_rods/yartsev_ch2_cantilever_reproduction.md): `BOOK_SLOPE_CLAMP_CONFIRMED` | Supports ordinary full run, smoke, quick boundary gate, and saved-data-only boundary source check. The quick gate is preliminary; the postprocess-only mode does not invoke the root solver. This is not a coupled-rod solver. |
| Yartsev Chapter-2 rigid angular joint | `analysis/anisotropic_rods/pilot_yartsev_ch2_coupled_rods.py` | `completed` small elastic diagnostic — `PASS` | [rigid-joint theory and gate](../docs/anisotropic_rods/yartsev_ch2_rigid_angular_joint.md) | Use for the declared sign/rank/virtual-work/limit gates and `beta=0,30,90 deg` HMS/DX-209 pilot only. It is not a stable baseline, final parameter study, complex-root workflow, EB comparison, or FEM validation. |
| Yartsev Chapter-2 rectangular orthotropic EB validation | `analysis/anisotropic_rods/validate_yartsev_ch2_rectangular_eb.py` | `completed` finite diagnostic — overall `PARTIAL_PASS`; targeted `FAIL_CONVERGENCE_ORDER` | [rectangular EB validation](../docs/anisotropic_rods/yartsev_ch2_rectangular_eb_validation.md) | Original fixed-64 `PARTIAL_PASS` is preserved. Proportional `(64,192)` closes the raw first-three accuracy threshold, but mode 1 violates the unchanged monotonicity allowance at the dense-eigensolver conditioning floor. No model coefficient or threshold changed; do not promote this to a general/off-axis EB model, stable baseline, unequal-thickness or 3D validation. |
| Stage-1 and fixed-epsilon EB applicability | `analysis/thickness_mismatch/audits/audit_eb_validity_vs_timoshenko_stage1.py`; `audit_eb_validity_fixed_epsilon_geometry_scan.py` | `completed` | [stage closure](../docs/thickness_mismatch/eb_safe_prefix_stage_closure.md) | Retain as source-generation provenance; these are not pending selector-development steps. |
| General spectrum-completeness audit | `analysis/thickness_mismatch/audits/audit_eb_timo_general_spectrum_completeness.py` | `superseded` | branch-informed gateway below; [research plan](../docs/thickness_mismatch/eb_safe_spectrum_prefix_research_plan.md#branch-informed-continuation-gateway-to-step-3) | Preserve the important `not_ready_for_step3` result; do not use it as the targeted gateway's current status. |
| Branch-informed completeness gateway | `analysis/thickness_mismatch/audits/audit_eb_timo_branch_continuation_gateway.py` | `completed` | [research plan](../docs/thickness_mismatch/eb_safe_spectrum_prefix_research_plan.md#branch-informed-continuation-gateway-to-step-3) | Canonical targeted K10/root-11 readiness provenance; not a general root-count proof. |
| Step 3A lower-envelope screen | `analysis/thickness_mismatch/audits/audit_eb_epsilon_lower_envelope_step3a.py` | `closed-research` | [stage closure](../docs/thickness_mismatch/eb_safe_prefix_stage_closure.md#43-step-3a-and-rejection-of-the-epsilon-lower-envelope) | Retain `S3_12`/`S3_14` evidence. The generated Step-3B manifest is an unexecuted proposal, not a next step. |
| Exact Rules A/B/S | `analysis/thickness_mismatch/postprocess/analyze_eb_rule_ab_exact_pareto.py` | `completed` | [stage closure](../docs/thickness_mismatch/eb_safe_prefix_stage_closure.md#5-exact-rule-abs-experiment-and-data-lock) | Canonical finite-set rule result. Rule B and Rule S are empirically equivalent on the checked data; no continuous guarantee. |
| Rule-S cost benchmark | `analysis/thickness_mismatch/benchmarks/benchmark_rule_s_cost_break_even.py` | `closed-research` | [stage closure](../docs/thickness_mismatch/eb_safe_prefix_stage_closure.md#12-canonical-cost-result) | Preserve `rule_S_cost_not_beneficial`; do not continue selector refinement without a material change in the engineering problem. |
| Counterexample frequency plots | `analysis/thickness_mismatch/maps/plot_counterexample_dimensional_frequency_beta.py` | `completed` | [frequency-map policy](../docs/thickness_mismatch/frequency_plot_computation_policy.md#current-counterexample-figures) | Use saved certified CSV data and `--plot-only` for presentation changes. |
| Timoshenko shape and suspect-mode diagnostics | `audit_timoshenko_shape_construction.py`; `audit_timoshenko_shape_bug_thin_limit.py`; `audit_timoshenko_modes_4_6_shape_diagnostics.py`; `audit_longitudinal_suspect_modes_eb_timo.py` | `active-diagnostic` | [thickness script map](analysis/thickness_mismatch/README.md#preferred-entry-points) | Use corrected display mappings; older invalidated global Timoshenko vectors are not physical evidence. |
| Out-of-plane EB plus torsion | out-of-plane map scripts and `analysis/thickness_mismatch/audits/compare_out_of_plane_analytic_vs_1d_fem.py` | `active-diagnostic` | [out-of-plane theory](../docs/theory/out_of_plane_eb_torsion.md) | Diagnostic subsystem only; absence of generated output in a clone is not proof of failure. |
| 3D FEM comparison | beta=0 validation wrappers and `scripts/analysis/solid_fem_coupled_equal_rods*.py` | `active-diagnostic` | [FEM validation status](../docs/thickness_mismatch/fem_validation_status.md) | Treat 3D models as independent engineering benchmarks, not tuned replicas of the ideal 1D joint. |
| Article-facing diagnostics | `scripts/analysis/plot_article_fig3_with_fp_and_ff_refs.py`; `validate_article_shape_cases_beta15.py` | `manual-review` | [article promotion workflow](../docs/writing/article_workflow.md) | Results-only diagnostics remain usable only with their documented prerequisites. Local `paper_*` workspaces are not part of this tracked checkout. |

## Active diagnostic workflows

Active work consists mainly of broad EB/Timoshenko maps, thickness-mismatch
branch/FEM checks, corrected Timoshenko shape interpretation, any separately
scoped anisotropic work beyond the completed small rigid-joint pilot,
out-of-plane maps/1D validation, and independent 3D FEM comparisons. Diagnostic
output belongs under `results/` until a canonical tracked document promotes a
conclusion.

## Completed and closed workflows

- The geometry-only epsilon certificate and straight lower-envelope hypothesis
  are `closed-research`; `epsilon_0` remains a screening quantity, not a
  certificate.
- Step 3A is `closed-research` after the confirmed `S3_12` and `S3_14`
  counterexamples.
- Exact Rules A/B/S are `completed`; the finite checked Rule-B safety result is
  retained separately from engineering cost.
- The Rule-S engineering selector is `closed-research` with
  `rule_S_cost_not_beneficial`. Rule S itself was not mathematically refuted.
- Focused veering/localization assessments that already support the tracked
  closure notes are completed provenance, even while the broader direction
  remains diagnostic.

## Historical and superseded workflows

| Workflow | Status | Canonical successor/report | Preservation rule |
| --- | --- | --- | --- |
| Historical fold Rules A--D, Rule A-gap, Rule C, and Rule D | `superseded` | exact A/B/S experiment and [stage closure](../docs/thickness_mismatch/eb_safe_prefix_stage_closure.md) | Keep the CSV-only postprocessor and old reports; do not present them as pending calibration. |
| Geometry-only epsilon pilot | `historical` | corrected branch-informed pilot plus stage closure | Preserve its selected-data provenance and limitations. |
| Pre-correction straight-baseline sign-scan workflow/cache | `superseded` | factorized corrected baseline recorded in the research plan | Never mix legacy and corrected cache identities. |
| General completeness `not_ready_for_step3` workflow | `superseded` for the targeted gate | branch-informed gateway | Retain the negative readiness evidence; the two audits answer different chronological questions. |
| Step-3B proposal | `historical` / unexecuted | Step-3A lower-envelope decision | Do not run or describe it as a pending next step under the closed research scope. |
| Rule-S selector refinement, extra guards, and larger unchanged cost grid | `closed-research` | final stage closure | Reopen only under the documented material-change conditions. |
| Original beta-15 solid-FEM implementation | `historical` | multi-beta `solid_fem_coupled_equal_rods.py` | Preserve until geometry, constraints, and output equivalence are audited. |

## Compatibility wrappers

These seven public paths are intentional wrappers and must remain available
until a separately approved migration proves their commands and outputs:

| Wrapper | Existing implementation |
| --- | --- |
| `scripts/run/run_beta_sweep_mu0_four_radii.py` | `scripts/plot_beta_sweep_mu0_four_radii_compare.py` |
| `scripts/run/run_mu_sweep_beta0_four_radii.py` | `scripts/plot_mu_sweep_beta0_four_radii_compare.py` |
| `scripts/run/run_mu_sweep_fixed_beta_four_radii.py` | `scripts/plot_mu_sweep_beta_fixed_four_radii_compare.py` |
| `scripts/run/run_mu_sweep_four_betas_analytic.py` | `scripts/plot_mu_sweep_radius_fixed_four_betas_analytic.py` |
| `scripts/run/run_branchwise_fem_audit.py` | `scripts/analyze_branchwise_fem_spectrum.py` |
| `scripts/run/run_tracked_bending_descendant_shapes_ru.py` | `scripts/plot_tracked_bending_descendant_shapes_ru.py` |
| `scripts/run/run_analytic_coupled_rods_vibration_shapes_beta15_mu06_eps0025_001_ru.py` | parameterized analytic shape runner/preset |

## Manual-review candidates

`scripts/analysis/thickness_mismatch/audits/audit_timoshenko_modes456_visualization.py`
is `manual-review`. The inventory found neither a Python importer nor a
documentation path reference, but it writes into an adjacent scientific audit
result family. It is not classified as unnecessary and must not be moved or
deleted until its provenance and output relationship are reviewed manually.
