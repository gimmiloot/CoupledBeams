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

## FEM-3A blocked first execution evidence

`results/nlsp_nonlinear_dynamic_3d_fem_pilot/a69310e3bb30bab7/` retains the
[blocked dynamic pilot](numerics/nlsp_nonlinear_dynamic_3d_fem_pilot.md): pinned
source/solver/code/environment hashes, installed-version protocol evidence,
frozen-p64 preflight and acceleration/safety checks, both inspected input decks,
one actual failed linear input/stdout/stderr/resource record, explicit attempt
ledger and read-only crash audit. No accepted preload/transient trajectory,
1D dynamic history, energy curve or comparison is fabricated. The historical
FEM-2R and every source result remain unchanged; no automatic rerun occurs.

## FEM-2R controlled continuation evidence

`results/nlsp_nonlinear_static_3d_fem_resume/210b74b8b166997c/` contains the
[completed qualified static comparison](numerics/nlsp_nonlinear_static_3d_fem_validation.md#fem-2r-controlled-continuation-after-input-serialization-failure):
pinned failed-parent/source manifests, reused1D profile hashes, serialization
checks/authorization ledger, six new actual INP/DAT/FRD/STA/logs, nodal U/RF/S/E,
material-section recovery41/81, independent reactions/current-moment balance,
linear/NL/correction tables/profiles, observed mesh/signal uncertainty, resource
counters and three PDF/PNG figures. All six cases and bounded diagnostics pass;
overall FEM2_DIAGNOSTIC_COMPLETE_WITH_QUALIFICATIONS. Historical failed6714 bundle
is unchanged; no source regeneration, fresh1D equilibrium or dynamic job.

## FEM-2 static preflight and stopped first job

`results/nlsp_nonlinear_static_3d_fem/6714bd9f2778e6d7/` preserves the
[partial static result](numerics/nlsp_nonlinear_static_3d_fem_validation.md):
frozen load/source manifests, p48/p64 linear/NL fields, residual/Hessian/reaction
and work diagnostics, attempted medium linear INP/DAT/STA/stdout/stderr and
original execution-code provenance. The input parser rejected an overlength
*STATIC number before equilibrium; no usable 3D static field or correction exists.
No historical mesh/result was regenerated and no hidden solver retry occurred.

## FEM-1R single refined-level evidence

`results/nlsp_linear_rectangular_3d_fem_refinement/63d44daae533389c/` contains
[the continuation](numerics/nlsp_linear_rectangular_3d_fem_validation.md#fem-1r-one-additional-mesh-refinement):
parent manifest/hash links,one.020 C3D10 mesh/24 complete eigenvectors,input/output,
four-mesh/signed-model-difference CSV,MAC/section/c/warping diagnostics,resources
and2 PDF/PNG figures. All8 preset mesh checks pass; this does not replace the
historical three-grid PARTIAL report or validate nonlinear V0.

## FEM-1 full linear family evidence

`results/nlsp_linear_rectangular_3d_fem/4262efa427b03dad/` contains the
[thick rectangular3D comparison](numerics/nlsp_linear_rectangular_3d_fem_validation.md):
pre-FEM frozen geometry/coefficients/window and1D completeness/profile data,
three audited C3D10 meshes and actual GEO/MSH/INP/DAT/FRD/logs,24 nodal vectors
per mesh, shape-only assignment/section diagnostics, raw/comparison CSV,
mesh convergence, axial effective contraction,3 PDF/PNG figures and49 tests.
No prior result is overwritten. Five modes remain MESH_UNRESOLVED;
identification PASS is not nonlinear V0 validation.

## Result directories

`results/nlsp_planar_physical_sanity_checks/888042e17cfc315a/` contains the
[bounded physical sanity diagnostic](theory/weakly_nonlinear_planar_physical_sanity_checks.md):
immutable source/action hashes, one new p64 half-amplitude tight history toT1,
all8 leading/asymptotic and amplitude comparisons, snapshots, retained strain/
resultant maxima, end forces/moments and finite-p/truncation balance accounting,
energy/mass/safety, CSV/JSON and3 PDF/PNG figures. Original execution code/manifest
and unchanged half-history hash are preserved after provenance/report-only
refinements. Overall DIAGNOSTIC_COMPLETE_WITH_QUALIFICATIONS; no independent
half-case p/time certification or physical-validation claim.

`results/planar_prepared_one_T1/795dcb14d3cd3a55/` contains the
[one-linear-period continuation](theory/planar_prepared_initial_state.md#prepared-one-period-feasibility):
immutable-source hashes and exact reused q0, shared95057 timestamps including all
old short times, three complete memory-mapped q/v histories, snapshots/observations,
energy/Loewner mass/safety and actual steps/counters; all8 spatial/temporal CSV/JSON,
instantaneous/cumulative difference curves, spatial argmax, full-horizon scales,
window tables, prefix regression and3 PDF/PNG figures. Exactly3 ODE; no new MP,
projection/BVP/eigen/symbolic audit. Primary numerical work232.845s, integration
150.293s. Execution/feasibility COMPLETED_EXPLORATORY_NOT_CERTIFIED, prefix/temporal/
energy-mass PASS, spatial PARTIAL. Original executionc505209d58fe74c0 is preserved
in execution manifest/summary/code snapshots: the subsequent coverage/caption
metadata revision changed no numerical arrays and performed zero integrations.
The old284a bundle is unchanged. Cached compute/report/plot perform zero new
numerical preparation or integrations; no automatic full5T1 extension.

`results/planar_prepared_feasibility/284a4039177391d1/` contains the separate
[precision/feasibility continuation](theory/planar_prepared_initial_state.md#prepared-precision-feasibility):
verified immutable sources, independent MP quadrature/projection probes,
original/new projection metrics, initial-only common endpoint-L2 coordinates,
explicit strict table and EXPLORATORY_NOT_CERTIFIED decision, exactly3 complete
0.1T1 histories, actual internal dt/counters, all8 spatial/temporal CSV/JSON,
energy/Loewner mass/safety, code/config snapshots and3 PDF/PNG figures.

Initial representation PASS; strict verification PARTIAL; exploratory runs
COMPLETED; spatial PARTIAL (theta_t max only); temporal PASS. Numerical cost
53.49s includes17.34s of integration, zero new BVP/model eigensolves. Matching
compute/report/plot does zero new preparation/integration; no full5T1 or p96
nonlinear run. Generated precision source inputs are preserved under
`results/planar_prepared_feasibility/precision_evidence/77d6a6db28bf677c/` and
copied into this bundle with their own manifest/hash checks. The older result
entries below retain their historical statuses.

`results/planar_prepared_initial_state/5ea8d41faf8ede54/` contains the
[prepared initial-state gate](theory/planar_prepared_initial_state.md): source
hashes, all saved-mode stat/harm coordinates and physical profiles, derivative
convergence, exact A/B endpoint audit, common U_star/C_star and quintic Theta3,
finite-amplitude residuals, projections32/48/64, initial energy/mass bounds,
old eight-component short norms on actual timestamps and two PDF/PNG figures.

Preparation/through-cubic-order compatibility PASS; common projection PARTIAL;
new temporal/spatial short checks NOT_RUN; overall pilot PARTIAL. Additional
relative strong/weak numerical checks remain qualified at p48/p64, with fixed
2e-12 gates and strict XFAIL tests. Primary audit3.857s,6 spectrum restorations,
4 direct validation solves,0 M-H/Timoshenko eigensolves and0 ODE integrations.
Matching compute/report/plot do zero BVP/eigen/history/ODE evaluations; report
never extends missing old prefixes or promotes old PARTIAL. That historical bundle contains no new trajectory
or full5T1 result. Early development preparation bundles f6c06a90f979185c
and b392d104d3391942 are retained separately, without overwriting any historical
source bundle.

`results/planar_second_order_axial_response/b3ea4eb6ac95d6e1/` contains the
[exact-time leading axial diagnostic](theory/planar_second_order_axial_response.md):
audited forcing expressions, common continuous bending background, complete
M-H spectral states, independent assembly/quadrature/exponential controls,
physical convergence and L2 projection data, historical nonlinear comparisons,
actual timestamps, sampling qualifications and three PDF/PNG figures.
Primary p16/24/32/48/64 coverage and one conditional p96 refinement are retained.
Charged numerical work was 373.599s within1200s: six primary M-H spectral
decompositions, two checkpoint restorations and zero ODE integrations.
The interrupted extension's spectral checkpoint is reused without a second p96
eigendecomposition; its metadata failure and resume provenance remain visible.

The diagnostic is COMPLETE, spatial convergence PARTIAL. Full-interval norms
cover 0...5T1; resolved leading-response traces use 0...0.1T1. Historical
small-amplitude p24 and nonlinear p48 comparisons retain their actual shorter
prefixes. Neither historical nonlinear bundle below is recomputed or promoted
to PASS. Matching compute and report/plot-only perform zero eigendecompositions,
exact-time response evaluations, symbolic derivations and ODE integrations.

`results/weakly_nonlinear_planar_recovery/054874a4a4c9c9ff/` contains the
[targeted continuation](theory/weakly_nonlinear_planar_time_pilot.md): validated
historical manifests, physical error localization and L2 projection/tail
data, continuous initial-compatibility evidence, preserved baseline helper,
pointwise equivalence, profiling, three short integration controls, fixed
budget decision and two PDF/PNG diagnostic figures. The original execution
bundle `db47d4efb6941bed` is retained. The transparent post-execution cache
revision records plot lookup corrections and safety/timestamp fixes in the
unexecuted full-run path. Original numerical data are unchanged; no controls
or trajectories were reintegrated.

The p48 full-run forecast was 863.33s against 859.59s remaining; charged
profiling/integration time was 40.41s. No full p48 trajectory was started:
`REFINEMENT_DEFERRED_BY_BUDGET`; solver recovery remains PARTIAL. Matching
`diagnose_weakly_nonlinear_planar_rod.py --compute` and report/plot-only read
saved evidence with zero integrations. The first-pilot bundle below, its
manifest and actual incomplete-case timestamps are unchanged.

`results/weakly_nonlinear_planar_time_pilot/<fingerprint>/` contains the
[first four-field time pilot](theory/weakly_nonlinear_planar_time_pilot.md):
atomic per-case coefficient/velocity NPZ histories, initial shape/projections,
linear controls, convergence, quartic energy/mass/domain diagnostics and up
to three PDF/PNG figures. Fixed numerical budgets preserve completed cases
if a later gate is unresolved. Identity includes protected model/reference
hashes, geometry, basis, quadrature, initial shape, integrator settings and
versions. Matching compute and plot/report-only perform zero integrations,
root solves and symbolic derivations. Old bundles are not overwritten.

Current bounded bundle: `c97287772bc461ef`, six full5T1 histories and one
budget-stopped prefix. The [tracked numerical report](theory/weakly_nonlinear_planar_time_pilot.md)
records overall PARTIAL despite passing temporal/energy checks; p24→32
spatial convergence fails for u,c and several velocities. Three figures use
only the two completed final-degree amplitude histories.

`results/weakly_nonlinear_spatial_rod/<fingerprint>/` contains the
[completed seven-field action audit](theory/weakly_nonlinear_spatial_rod.md):
exact rational quartic action, two cubic residual representations, all21
coefficient differences, supplied-draft comparison, boundary/energy/symmetry
identities, manufactured jets/amplitude orders, limited linear/split profiles
and source/input/code/version/Git provenance. `current.json` selects a matching
bundle. `--compute` reuses identity/artifact-checked results; `--report-only`
performs zero derivations and roots. A tracked generated appendix preserves
the full expansion. Source page images remain in sibling `source_audit/`.
No nonlinear time evolution, threshold, modal classification or parameter map.

`results/mindlin_herrmann_timoshenko_lambda_beta_large_checks/356c3be4953268f2/`
contains the completed [production geometry checks](theory/mindlin_herrmann_timoshenko_lambda_beta_large_checks.md):
fixed Lambda reference,37-angle spectra, shared baseline, paired-section
counts/brackets/quality/Gram coefficients, direct/stepped/swap controls,
68 independent QR reviews, performance and long CSV. Two3-panel figures,
PDF+PNG, use independently sorted positions1–12; guard13 is saved only.
Run `python scripts/analysis/verify_mindlin_herrmann_timoshenko_lambda_beta_large_checks.py --compute`;
redraw via `--plot-only results/mindlin_herrmann_timoshenko_lambda_beta_large_checks/356c3be4953268f2`.
Matching cache/plot-only performs zero roots. All previous bundles remain
unchanged; prototype132358083da36937 is retained, final bundle adds Gram
and independent anomaly reviews. No theory comparison/tracking/third figure.

`results/coupled_longitudinal_theory_thickness_screening/7be0fce968fd2b35/`
contains the completed [thickness screening](theory/coupled_longitudinal_theory_thickness_screening.md):
exact five-thickness/three-angle input and scaling audit, all45 inventories,
15 direct beta0 frequency/profile checks,180 comparison rows, full geometric
overlaps, contraction norms, gaps, profiles/counts/brackets/diagnostics and
three figures. Sub-cutoff certificates cover12+guard13; configured3.75
ceiling does not assert coverage of unsupported optical tails. Run
`python scripts/analysis/screen_coupled_longitudinal_theory_thickness.py --compute`;
redraw with `--plot-only results/coupled_longitudinal_theory_thickness_screening/7be0fce968fd2b35`.
Fingerprint/artifact checks give zero-root matching reuse; the previous
hierarchy/MH bundles remain immutable. No energy classes or across-case tracking.

`results/coupled_longitudinal_theory_hierarchy_screening/9a35ff23c3f43c66/`
contains the completed [bounded hierarchy screening](theory/coupled_longitudinal_theory_hierarchy_screening.md):
24 certified model/angle cases, comparator gates, 12+guard13 profiles,
96-row frequency/difference table, all fixed-beta 12x12 geometric overlaps,
theta small-norm records, contraction norms/D_c, adjacent gaps, diagnostics
and three saved-data plots. Run
`python scripts/analysis/screen_coupled_longitudinal_theory_hierarchy.py --compute`;
matching identity/artifact hashes reuse data with zero roots. Redraw with
`--plot-only results/coupled_longitudinal_theory_hierarchy_screening/9a35ff23c3f43c66`.
Old MH references are immutable. Earlier local implementation attempts remain
under the same result root for provenance; the named bundle is canonical.
No energy classification, tracking or applicability certification.

`results/mindlin_herrmann_timoshenko_general_beta_joint/<fingerprint>/`
contains the completed [general-frame gate](theory/mindlin_herrmann_timoshenko_rigid_joint.md#9-general-angle-geometry-действующий-project-contract):
transforms/ranks/duality, frozen beta0/zero-limit/swap/reflection metrics,
certified arm-pole catalogs and Schur counts, full coupled brackets/SVD,
frequency and full local-state profiles. The fixed5/45/90 equal-arm G20
pilot reports12+guard13; bounded complete inventories23/24/24. Run
`python scripts/analysis/verify_mindlin_herrmann_timoshenko_general_beta_joint.py --compute`.
Manifest checks input/code/source/frozen/artifact hashes; matching reuse
performs zero root evaluations. No parameter map, hierarchy or tracking;
old beta0 bundle3059d70b1b50ea2e is retained unchanged.

`results/mindlin_herrmann_timoshenko_beta0_joint/<fingerprint>/` is the
completed [reduced joint transparency gate](theory/mindlin_herrmann_timoshenko_rigid_joint.md):
beta0 only, three fixed splits of G20,total L=1, immutable direct reference
reuse, inventories/brackets/counts, matrices, virtual work, frequencies and
all eight state profiles, MAC/L2/c/R/interface/arm-swap diagnostics. Run
`python scripts/analysis/verify_mindlin_herrmann_timoshenko_beta0_joint.py --compute`.
Manifest validates inputs/source/code/reference/artifact hashes; matching
reuse performs zero root evaluations. Source page snapshots are in sibling
`source_audit/`. All separate gates PASS; no nonzero-angle evidence implied.

All 23 immediate ignored result directories present in the pre-refactor
snapshot are included here, together with later registered diagnostic result
families.

| Result directory | Scientific workflow | Status | Canonical report or documentation | Reproduction entry point | Local archive class |
| --- | --- | --- | --- | --- | --- |
| `results/mindlin_herrmann_timoshenko_single_rod/` | One finite G20 straight CC rod, selected Jang project closure, elementary/planar RL/MH hierarchy | finite spectrum `PASS`; hierarchy `PARTIAL_PASS` for resolved-contraction clamp interpretation | [canonical §14--20](theory/mindlin_herrmann_timoshenko_single_rod.md#14-production-formulation-decision) | `scripts/analysis/verify_mindlin_herrmann_timoshenko_single_rod.py --compute`; matching hash-validated reuse | Exact parameters, source/input/code hashes, roots/brackets/conditioning/count certificates, independent QR-expm checks, mass-normalized profiles and hierarchy CSV; source_audit/ holds page images. No map or coupled model. |
| `results/mindlin_herrmann_timoshenko_literature/` | One-rectangle Rucka/Jang source dispersion and one controlled comparison; separate `rectangular_preset_audit/` page/provenance evidence | `completed` diagnostic — `MHTIM_VARIANT_DEPENDENT`; new preset gate `RECTANGULAR_MH_PRESET_MAPPING_UNRESOLVED`; production coefficients unresolved | [canonical audit](theory/mindlin_herrmann_timoshenko_single_rod.md) and §13 | `scripts/analysis/reproduce_mindlin_herrmann_timoshenko_literature.py --compute --case all --jang-kappa 5/6`; plot-only reads saved data. New gate requires only source-check/tests, no scientific compute or plot | Old bundles retained. New evidence records Ng normal-block difference and missing Fernandes PDF; no production dispersion run. Jang kappa remains explicit conditional input, not source/default. |
| `results/joint_review/circular_eb_spring_general_spectrum/` | Circular EB low spectrum, forms and bounded kappa continuation | `completed`; 6/6 local descendants, 6+guard prefix only | [durable scientific note](laminated_beams/circular_eb_rotational_spring_rigid_limit.md) | `scripts/analysis/joint_review/check_circular_eb_spring_spectrum.py`; spectrum / shapes / kappa-continuation modes in [script inventory](../scripts/README.md#circular-eb-springrigid-reviewer-diagnostic) | Preserve original endpoints/direct failures, continuation evidence and higher-spectrum qualification. Generated local evidence may be absent in a fresh clone. |
| `results/joint_review/rotational_spring_rigid_trend_pilot/` | Earlier close-candidate audit and predeclared geometry retries | `historical`; stopped unresolved, not current preferred workflow | [history and limitations](laminated_beams/circular_eb_rotational_spring_rigid_limit.md#история-остановок-и-спектральный-результат), [K01/K02](memory/knowledge.md#eb-joint-k01) | Historical one-shot provenance in local diagnostics; no permanent runner | Preserve candidates and audit/retry reports; no physical-multiplicity conclusion. Ignored local evidence may be absent in a fresh clone. |
| `results/timoshenko_bishop_single_rod/` | Exact displacement/energy audit of one rectangular rod | `completed`; `COMBINED_KINEMATICS_NOT_UNIQUELY_DEFINED` | [candidate fields, mixed terms and hard gate](theory/timoshenko_bishop_single_rod.md) | `scripts/analysis/audit_timoshenko_bishop_single_rod.py --compute` | Preserve source images/text, content-addressed audit/manifest and targeted-test XML. Algebra only; no spectrum, cache reads or figures. |
| `results/bishop_literature/` | Fixed circular longitudinal literature reproduction | `completed`; numerical checks pass, source-print/figure qualifications retained | [canonical report](theory/bishop_literature_reproduction.md) | `scripts/analysis/reproduce_bishop_literature.py --compute --case all` | Preserve source audit and fingerprint bundles; `current.json` selects current data. Manifests, source precision, frequencies, profiles, residual/energy/Gram checks, independent 50/70 dps and two publication-comparison figures. |
| `results/laminated_beams/inplane_kelvin_voigt_eb_rlb_complex_confirmation/` | Fixed R0–R4 weak-damping cross-theory confirmation | `COMPLEX_CONFIRMATION_COMPLETED` | [D22/K23 report](laminated_beams/inplane_kelvin_voigt_eb_rlb_complex_confirmation.md) | `scripts/analysis/laminated_beams/confirm_inplane_kelvin_voigt_eb_rlb.py --compute` (missing-only) | Preserve six new roots, five comparison rows, ten-state diagnostics, complex forms and K12/K15/K19/K22 provenance; no auxiliary roots or solver changes |
| `results/laminated_beams/inplane_kelvin_voigt_rlb_elastic_screening/` | Sparse real RLB screening and same-angle comparison with read-only EB K15 | `completed` — 24 states / 24 confirmed pairs | [D21/K22 report](laminated_beams/inplane_kelvin_voigt_rlb_elastic_screening.md) | `scripts/analysis/laminated_beams/screen_inplane_kelvin_voigt_rlb_elastic.py --compute` | Preserve modal/comparison CSV, forms, diagnostics and source hashes; 14 reused target states, 10 new, no positive-d roots |
| `results/laminated_beams/inplane_kelvin_voigt_rlb_solver_architecture/` | Bounded RLB production routing regression and exact EB coefficient limit | `completed` technical validation — `RLB_KV_PRODUCTION_PASS` | [tracked D20/K21 report](laminated_beams/inplane_kelvin_voigt_rlb_solver_architecture.md) | `scripts/analysis/laminated_beams/check_inplane_kelvin_voigt_rlb_solver_architecture.py --compute` | Preserve 11 control rows, matrix/limit diagnostics and source hashes; no new physical parameter study |
| `results/laminated_beams/inplane_kelvin_voigt_targeted_weak_damping_completion/` | Complete the original three-state/two-d comparison | `SIX_STATE_COMPARISON_COMPLETED`; separate full C qualification | [six states and commands](laminated_beams/inplane_kelvin_voigt_targeted_weak_damping_completion.md) | `scripts/analysis/laminated_beams/complete_inplane_kelvin_voigt_weak_damping.py --compute` (missing-only) | `preserve-diagnostic-evidence`; two new roots, four reused rows, complex forms and protected K15–K18 provenance |
| `results/laminated_beams/inplane_kelvin_voigt_solver_architecture/` | Fixed reduced/full EB KV regression | `completed-technical`; reduced PASS, full with local raw rank qualification C | [controls, scope and commands](laminated_beams/inplane_kelvin_voigt_solver_architecture.md) | `scripts/analysis/laminated_beams/check_inplane_kelvin_voigt_solver_architecture.py --compute` (missing-only) | `preserve-diagnostic-evidence`; 12 control rows, source hashes, diagnostics and normalized two-arm forms |
| `results/laminated_beams/inplane_kelvin_voigt_ac_diagnostics/` | Exact complex symmetry blocks for the same K16 A/C .001 | `completed-diagnostic`; `FULL_TRANSFER_RECOVERY_CONDITIONING` for both, original K16 qualifications retained | [report, residuals and commands](laminated_beams/inplane_kelvin_voigt_ac_diagnostics.md) | `scripts/analysis/laminated_beams/diagnose_inplane_kelvin_voigt_ac.py --compute` (missing-only) | `preserve-diagnostic-evidence`; local CSV/JSON/NPZ, protected K15/K16 sources |
| `results/laminated_beams/inplane_kelvin_voigt_targeted_weak_damping/` | Three prescribed ACTIVE EB seeds, two small d values | `PARTIAL_NUMERICAL_QUALIFICATIONS`: two accepted, two rejected candidates, two targets not run | [report, diagnostics and commands](laminated_beams/inplane_kelvin_voigt_targeted_weak_damping.md) | `scripts/analysis/laminated_beams/check_inplane_kelvin_voigt_targeted_weak_damping.py --compute` | `preserve-diagnostic-evidence`; all attempts retained, prior K12/K15 sources immutable |
| `results/laminated_beams/inplane_kelvin_voigt_elastic_screening/` | Four-angle EB elastic participation and symmetry screening | `completed` — 24 accepted states | [screening report and commands](laminated_beams/inplane_kelvin_voigt_elastic_screening.md) | `scripts/analysis/laminated_beams/screen_inplane_kelvin_voigt_elastic.py --compute` (missing-only) | `preserve-diagnostic-evidence`; local CSV/JSON/NPZ; K11/K12 reused read-only; no new positive-d roots |
| `results/laminated_beams/inplane_kelvin_voigt_literature_benchmarks/` | Failla Table1 / Hong Tables2–3 source-specific external KV checks | First pass `PARTIAL` unchanged; D13 second pass `PASS_WITH_SOURCE_PRINT_QUALIFICATIONS` | [tracked comparison, rounding and commands](laminated_beams/inplane_kelvin_voigt_literature_benchmarks.md) | `scripts/analysis/laminated_beams/benchmark_inplane_kelvin_voigt_literature.py --second-pass` | `preserve-diagnostic-evidence`; separate `*_second_pass*` and precision-audit files; local CSV/JSON/NPZ, not included in a fresh clone |
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
| `results/eb_timo_counterexample_dimensional_frequency_beta/` | Certified dimensional-frequency beta plots for `S3_12` and `S3_14` | `active-diagnostic` | generated `counterexample_dimensional_frequency_beta_report.md`; [frequency-map policy](numerics/frequency_map_computation_policy.md) | `scripts/analysis/thickness_mismatch/maps/plot_counterexample_dimensional_frequency_beta.py` | `preserve-certified-figure-data` |
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
| `results/anisotropic_rods/yartsev_ch2_free_free/` | Yartsev Chapter-2 one-rod corrected/printed comparison and Figure-2.2 graph-resolution gate | `completed` source reproduction | generated `single_rod_reproduction_report.md`; [tracked free-free note](anisotropic_rods/yartsev_ch2_single_rod_reproduction.md) | `scripts/analysis/anisotropic_rods/reproduce_yartsev_fig_2_2.py` | Local generated evidence; not guaranteed in Git or a fresh clone |
| `results/anisotropic_rods/yartsev_ch2_cantilever/` | Full orientation/length cantilever reproduction and diagnostic evidence for both clamp variants | `completed` source reproduction / diagnostic evidence | generated `cantilever_reproduction_report.md`; [tracked cantilever note](anisotropic_rods/yartsev_ch2_cantilever_reproduction.md) | `scripts/analysis/anisotropic_rods/reproduce_yartsev_ch2_cantilever.py` | Local generated evidence; not guaranteed in Git or a fresh clone |
| `results/anisotropic_rods/yartsev_ch2_cantilever_quick_gate/` | Preliminary elastic sensitivity of the two cantilever clamp variants | `completed` preliminary boundary-sensitivity diagnostic | generated `quick_boundary_gate_report.md`; [tracked cantilever note](anisotropic_rods/yartsev_ch2_cantilever_reproduction.md) | cantilever CLI with `--quick-boundary-gate` | Local generated evidence; the quick gate does not replace source reproduction |
| `results/anisotropic_rods/yartsev_ch2_cantilever_boundary_source_check/` | Saved-data-only Figure-2.8 frequency/loss comparison used to identify the book clamp | `completed` source-boundary decision evidence | generated `boundary_source_check_report.md`; [tracked cantilever note](anisotropic_rods/yartsev_ch2_cantilever_reproduction.md) | cantilever CLI with `--postprocess-boundary-source-check` | `BOOK_SLOPE_CLAMP_CONFIRMED`; local generated evidence, not guaranteed in Git |
| `results/anisotropic_rods/yartsev_ch2_coupled_joint_pilot/` | Ideal rigid angular-joint sign, virtual-work, limit, straight-rod-equivalence, and small elastic spectrum gate | `completed` diagnostic pilot — `PASS` | generated `coupled_joint_pilot_report.md`; [tracked rigid-joint note](anisotropic_rods/yartsev_ch2_rigid_angular_joint.md) | `scripts/analysis/anisotropic_rods/pilot_yartsev_ch2_coupled_rods.py` | Local generated evidence; first six roots plus seventh guard at `beta=0,30,90 deg`; not a parameter study or stable baseline |
| `results/anisotropic_rods/yartsev_ch2_rectangular_eb_validation/` | `theta=0` rectangular EB/Saint-Venant exact, unequal-length, Timoshenko-limit, and independent 1D-FEM gate | `completed` finite diagnostic — overall `PARTIAL_PASS`; targeted `FAIL_CONVERGENCE_ORDER` | original generated `rectangular_eb_validation_report.md` plus `targeted_refinement_report.md`; [tracked validation note](anisotropic_rods/yartsev_ch2_rectangular_eb_validation.md) | `scripts/analysis/anisotropic_rods/validate_yartsev_ch2_rectangular_eb.py` | Local ignored evidence; original fixed-64 `PARTIAL_PASS` is preserved. Raw proportional `(64,192)` closes the first-three target (`6.18e-6`), but mode 1 violates the unchanged monotonic-refinement allowance at the conditioning floor; no coefficient or threshold changed. |

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
