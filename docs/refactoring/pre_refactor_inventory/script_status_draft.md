# Runnable-script status draft

This is a status proposal, not a move plan. `safe_to_run` is deliberately
conservative: “no” means the normal command can solve roots/FEM or write
`results/`; it does not mean the script is defective. Plot-only/CSV-only modes
are listed separately where their contract is documented.

## Main workflows

| Path | Status | Replacement or successor | Canonical document | Safe to run in inventory task | Expected outputs | Root calculation level | Archive candidate | Reason |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `src/my_project/analytic/FreqFromAngle.py` | `stable` | none | `README.md`, theory equations | no | angle-sweep plot/data | EB roots | retain | baseline module and imported entry point |
| `src/my_project/analytic/FreqFromMu.py` | `stable` | wrapper/presentation scripts reuse it | `scripts/README.md` | no | frequency-unit `mu` sweep | EB roots + tracking | retain | baseline and imported by other code |
| `src/my_project/analytic/FreqMuNet.py` | `stable` | `scripts/run/run_lambda_mu_fixed_beta_analytic.py` is convenient CLI, not replacement | `scripts/README.md` | no | `Lambda(mu)` plot/CSV | EB roots + tracking | retain | imported compatibility/public module |
| `src/my_project/fem/python_fem.py` | `stable` | none | root README and AGENTS FEM rule | no | `results/fem_spectrum.csv` | baseline FEM eigensolve | retain | production baseline and 19 importers |
| `scripts/run/run_beta_sweep_mu0_four_radii.py` | `compatibility-wrapper` | `scripts/plot_beta_sweep_mu0_four_radii_compare.py` | `scripts/README.md` | no | beta-sweep results | analytic/FEM | retain wrapper | documented command |
| `scripts/run/run_mu_sweep_beta0_four_radii.py` | `compatibility-wrapper` | `scripts/plot_mu_sweep_beta0_four_radii_compare.py` | `scripts/README.md` | no | mu-sweep results | analytic/FEM | retain wrapper | documented command |
| `scripts/run/run_mu_sweep_fixed_beta_four_radii.py` | `compatibility-wrapper` | `scripts/plot_mu_sweep_beta_fixed_four_radii_compare.py` | `scripts/README.md` | no | fixed-beta mu results | analytic/FEM | retain wrapper | documented command |
| `scripts/run/run_mu_sweep_four_betas_analytic.py` | `compatibility-wrapper` | `scripts/plot_mu_sweep_radius_fixed_four_betas_analytic.py` | `scripts/README.md` | no | analytic mu results | EB roots | retain wrapper | documented command |
| `scripts/run/run_branchwise_fem_audit.py` | `compatibility-wrapper` | `scripts/analyze_branchwise_fem_spectrum.py` | `scripts/README.md` | no | branchwise FEM CSVs | FEM eigensolves | retain wrapper | compatibility and provenance |
| `scripts/run/run_tracked_bending_descendant_shapes_ru.py` | `compatibility-wrapper` | `scripts/plot_tracked_bending_descendant_shapes_ru.py` | `scripts/README.md` | no | tracked-shape figures | EB roots/tracking | retain wrapper | public command |
| `scripts/run/run_analytic_coupled_rods_vibration_shapes_beta15_mu06_eps0025_001_ru.py` | `compatibility-wrapper` | mode-shape runner | `scripts/README.md` | no | one-case shape figures | EB roots | soft archive candidate | one-case preset; wrapper must remain until usage review |
| `scripts/analyze_target_descendants_beta15_r5.py` | `closed-research` | none | `docs/veering/` assessments | no | descendant assessment tables | may consume/compute tracked spectra | soft archive candidate | completed veering evidence line |
| `scripts/analyze_flat_mu_bending_energy.py` | `closed-research` | none | `docs/veering/mu_slow_evolution_assessment.md` | no | arm-energy CSV | FEM/shape postprocess | soft archive candidate | closure evidence; also imported |
| `scripts/analysis/audit_desc6_veering_localization.py` | `closed-research` | none | veering assessments | no | localization CSV/plots | diagnostic roots/shapes | soft archive candidate | one-case research audit |
| `scripts/analysis/check_thickness_mismatch_eta_zero_limit.py` | `stable` | none | `docs/thickness_mismatch/README.md` | no | eta-limit report/CSV | determinant/root checks | retain | mathematical regression workflow |
| `scripts/analysis/audit_variable_length_timoshenko_limits.py` | `stable` | tau-aware companion for eta!=0 is separate | thickness README | no | root limit audit | EB/Timoshenko roots | retain | eta=0 baseline diagnostic |
| `scripts/analysis/audit_variable_length_thickness_timoshenko_limits.py` | `active-diagnostic` | none | thickness README | no | tau-aware limit audit | EB/Timoshenko roots | retain | eta!=0 sorted-root checks |
| `scripts/analysis/thickness_mismatch/audits/audit_eb_validity_vs_timoshenko_stage1.py` | `closed-research` | corrected pilot/gateway/Step-3A evidence chain | safe-prefix closure | no | Stage-1 CSVs/cache/report | heavy EB/Timoshenko | soft archive candidate | source study completed; provenance still canonical |
| `scripts/analysis/thickness_mismatch/audits/audit_eb_validity_fixed_epsilon_geometry_scan.py` | `closed-research` | same evidence chain | safe-prefix closure | no | fixed-epsilon CSVs/cache/report | heavy EB/Timoshenko | soft archive candidate | source study completed; predictors still reused |
| `scripts/analysis/thickness_mismatch/audits/run_eb_epsilon_apriori_pilot.py` | `closed-research` | branch-informed corrected pilot | research plan | no | 21-case pilot | heavy EB/Timoshenko | soft archive candidate | original source generation retained for provenance |
| `scripts/analysis/thickness_mismatch/audits/audit_eb_epsilon_baseline_thresholds.py` | `closed-research` | none | baseline report/research plan | no; plot-only only with saved CSVs | threshold CSVs/cache/plots | heavy factorized/general roots | retain | corrected reference thresholds and legacy audit |
| `scripts/analysis/thickness_mismatch/audits/audit_eb_timo_general_spectrum_completeness.py` | `superseded` | branch-informed gateway for targeted Step 3 | research plan | no; plot-only only with saved data | completeness audit | heavy roots/SVD | soft archive candidate | important negative readiness result; not canonical production solver |
| `scripts/analysis/thickness_mismatch/audits/audit_eb_timo_branch_continuation_gateway.py` | `closed-research` | none | gateway report/research plan | no; plot-only only with saved data | gateway CSV/cache/report | heavy continuation | retain | canonical K10 readiness provenance |
| `scripts/analysis/thickness_mismatch/audits/audit_eb_epsilon_lower_envelope_step3a.py` | `closed-research` | none | Step-3A report and closure | no; plot-only is zero-root | Step-3A audit/results | heavy continuation in normal mode | retain | canonical counterexample evidence |
| `scripts/analysis/thickness_mismatch/postprocess/analyze_eb_safe_prefix_certification.py` | `superseded` | exact A/B/S experiment | historical research-plan section | CSV-only but writes results | fold certification outputs | zero roots | soft archive candidate | historical Rules A–D workflow |
| `scripts/analysis/thickness_mismatch/postprocess/analyze_eb_epsilon_apriori_pilot.py` | `closed-research` | exact A/B/S for final selector decision | research plan | CSV-only but writes results | pilot analysis | zero roots | soft archive candidate | geometry-only screening completed |
| `scripts/analysis/thickness_mismatch/postprocess/analyze_eb_rule_ab_exact_pareto.py` | `closed-research` | none | exact report and closure | no; normal mode writes results | exact A/B/S CSVs/report | zero roots; EB reconstruction/SVD | retain | canonical Phase-I result |
| `scripts/analysis/thickness_mismatch/benchmarks/benchmark_rule_s_cost_break_even.py` | `closed-research` | none | cost report and closure | no; plot-only only with saved data | five-case cost results | Timoshenko roots in normal mode | retain | canonical negative engineering result |
| `scripts/analysis/thickness_mismatch/maps/plot_counterexample_dimensional_frequency_beta.py` | `closed-research` | none | computation policy/Step-3A | no; plot-only is zero-root | two PDFs and audit CSVs | heavy continuation in normal mode | retain | certified S3 visualization |
| `scripts/analysis/thickness_mismatch/maps/plot_eb_vs_timoshenko_lambda_beta_cases.py` | `active-diagnostic` | none | thickness script map | no; plot-only with cache/data | beta maps/cache/report | EB/Timoshenko roots | retain | current sorted-frequency diagnostic and imported helper |
| `scripts/analysis/thickness_mismatch/maps/plot_eb_vs_timoshenko_lambda_mu_cases.py` | `active-diagnostic` | none | thickness script map | no; plot-only with cache/data | mu maps/cache/report | EB/Timoshenko roots | retain | current diagnostic |
| `scripts/analysis/plot_mode_shapes_eta_beta_scan.py` | `active-diagnostic` | proposed future parameterized sorted/descendant CLI | thickness script-map TODO | no | descendant shape figures/CSV | EB roots/tracking | hard archive only after successor audit | parallel implementation |
| `scripts/analysis/plot_mode_shapes_eta_beta_scan_sorted_modes.py` | `active-diagnostic` | same proposed CLI | thickness script-map TODO | no | sorted shape figures/CSV | EB roots | hard archive only after successor audit | parallel implementation |
| `scripts/analysis/thickness_mismatch/maps/plot_out_of_plane_lambda_beta_mu_eta.py` | `active-diagnostic` | none | out-of-plane theory note | no | out-of-plane maps/report | out-of-plane roots | retain | generated canonical outputs absent locally |
| `scripts/analysis/thickness_mismatch/maps/plot_out_of_plane_mode_character_beta.py` | `active-diagnostic` | none | out-of-plane theory note | no | mode-character maps | out-of-plane roots + 1D FEM matrices | retain | diagnostic only |
| `scripts/analysis/thickness_mismatch/audits/compare_out_of_plane_analytic_vs_1d_fem.py` | `active-diagnostic` | none | out-of-plane theory section 11 | no | validation CSV/report | analytic roots + 1D FEM | retain | independent 1D validation workflow |
| `scripts/analysis/thickness_mismatch/audits/compare_full_spectrum_analytic_vs_3d_fem.py` | `active-diagnostic` | none | thickness script map | no | analytic union/3D comparison | analytic roots; optional 3D FEM | retain | scaffold, not article validation |
| `scripts/analysis/solid_fem_coupled_equal_rods.py` | `active-diagnostic` | none | FEM validation status | no | multi-beta solid-FEM outputs | external Gmsh/CalculiX FEM | retain | active baseline diagnostic |
| `scripts/analysis/solid_fem_coupled_equal_rods_beta15.py` | `historical` | multi-beta script | FEM validation status | no | original beta-15 outputs | external FEM | hard archive after dependency/output audit | reproduction path for original output |
| `scripts/analysis/solid_fem_coupled_equal_rods_point_joint.py` | `active-diagnostic` | none | FEM validation status | no | point-joint/matching outputs | external FEM | retain | current joint diagnostic |
| `scripts/analysis/plot_article_fig3_with_fp_and_ff_refs.py` | `one-off` | no in-repository successor | `scripts/README.md` | no | results-only article diagnostic | analytic roots/data | manual review required | imports absent local article workspace |

## Import audit summary

The complete per-file role/status data are in
`manifests/python_file_inventory.csv`; exact edges are in
`manifests/import_edges.csv`.

Forty runnable-looking files are also imported. They must not be moved as
standalone CLIs without an import migration:

```text
scripts/analysis/audit_lambda_beta_sorted_descendant_thickness_mismatch.py
scripts/analysis/audit_mu_scan_eta0_first6_rearrangement.py
scripts/analysis/check_analytic_shape_in_fem_residual.py
scripts/analysis/compare_analytic_fem_tracked_descendant_shape.py
scripts/analysis/compare_coupled_equal_rods_eb_timoshenko.py
scripts/analysis/compare_single_rod_eb_timoshenko.py
scripts/analysis/plot_desc05_full_shapes_beta15_eps_sweep.py
scripts/analysis/plot_mode_shapes_eta_beta_scan.py
scripts/analysis/plot_thickness_mismatch_branch_shapes_vs_eta.py
scripts/analysis/thickness_mismatch/audits/analyze_universal_eb_validity_parameters_stage1.py
scripts/analysis/thickness_mismatch/audits/audit_eb_epsilon_baseline_thresholds.py
scripts/analysis/thickness_mismatch/audits/audit_eb_epsilon_lower_envelope_step3a.py
scripts/analysis/thickness_mismatch/audits/audit_eb_timo_branch_continuation_gateway.py
scripts/analysis/thickness_mismatch/audits/audit_eb_timo_general_spectrum_completeness.py
scripts/analysis/thickness_mismatch/audits/audit_eb_validity_fixed_epsilon_geometry_scan.py
scripts/analysis/thickness_mismatch/audits/audit_eb_validity_vs_timoshenko_stage1.py
scripts/analysis/thickness_mismatch/audits/audit_timoshenko_shape_construction.py
scripts/analysis/thickness_mismatch/audits/run_eb_epsilon_apriori_pilot.py
scripts/analysis/thickness_mismatch/benchmarks/benchmark_rule_s_cost_break_even.py
scripts/analysis/thickness_mismatch/maps/plot_counterexample_dimensional_frequency_beta.py
scripts/analysis/thickness_mismatch/maps/plot_eb_vs_timoshenko_lambda_beta_cases.py
scripts/analysis/thickness_mismatch/maps/plot_lambda_mu_beta90_eta0p1_with_single_rod_refs.py
scripts/analysis/thickness_mismatch/postprocess/analyze_eb_epsilon_apriori_pilot.py
scripts/analysis/thickness_mismatch/postprocess/analyze_eb_rule_ab_exact_pareto.py
scripts/analysis/thickness_mismatch/postprocess/analyze_eb_safe_prefix_certification.py
scripts/analysis/track_lambda_eta_thickness_mismatch.py
scripts/analyze_branchwise_fem_spectrum.py
scripts/analyze_flat_mu_bending_energy.py
scripts/compare_beta0_analytic_vs_fem.py
scripts/plot_beta_sweep_mu0_compare.py
scripts/plot_beta_sweep_mu0_four_radii_compare.py
scripts/plot_mu_sweep_beta_fixed_four_radii_compare.py
scripts/plot_mu_sweep_beta0_four_radii_compare.py
scripts/plot_mu_sweep_radius_fixed_four_betas_analytic.py
scripts/plot_tracked_bending_descendant_shapes_ru.py
scripts/sweep_grid_policy.py
src/my_project/analytic/FreqFromAngle.py
src/my_project/analytic/FreqFromMu.py
src/my_project/analytic/FreqMuNet.py
src/my_project/fem/python_fem.py
```

The root-level `Freq*.py` modules and root-level `scripts/*.py` presentation
programs are therefore compatibility modules as well as runnable files. The
seven thin `scripts/run/` wrappers are explicit compatibility entry points.

Ninety-six runnable files have no Python importer but do have documentation
path/CLI references. The CSV inventory is the authoritative full list; these
include the 3D FEM audits, eta/crossing audits, map scripts, and documented
limit checks. `imported_by_count=0` is expected for such commands.

Only one runnable file has neither an importer nor a detected documentation
path reference:

```text
scripts/analysis/thickness_mismatch/audits/audit_timoshenko_modes456_visualization.py
```

It is a possible orphan, not an archive decision.

Potential duplicate entry-point families are the seven wrapper/implementation
pairs, the two beta-scan mode-shape programs, the beta-15 and multi-beta solid
FEM programs, and the historical fold certification versus exact A/B/S
postprocessors. Their output contracts differ and must be compared before any
consolidation.

