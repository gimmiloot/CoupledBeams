# Archive-candidate draft

Nothing was moved, renamed, or deleted. The word “archive” below means a later
status/navigation decision unless a separate dependency-audited `git mv` is
explicitly approved. Every row requires manual review.

## A. Do not archive

| Path | Current role | Reason | Dependencies | Canonical replacement | Risk of archiving | Required checks before move |
| --- | --- | --- | --- | --- | --- | --- |
| `docs/theory/equations.tex` | source-of-truth theory | verified formulas and unknown ordering | all analytic work | none | critical scientific corruption | no move proposed |
| `docs/theory/assumptions.md` | source-of-truth assumptions | notation, limits, physical interpretation | theory/code/tests | none | critical | no move proposed |
| `docs/literature/source_index.md` | source provenance | includes literature warnings | theory audit | none | high | no move proposed |
| `src/my_project/analytic/formulas.py` | baseline determinant | 30 importers; frozen signs/order | analytic scripts/tests | none | critical | no move proposed |
| `src/my_project/analytic/solvers.py` | baseline solver layer | 24 importers and numerical contracts | analytic scripts/tests | none | critical | no move proposed |
| `src/my_project/analytic/formulas_thickness_mismatch.py` | eta determinant/geometry | 50 importers; eta=0 and mass invariants | thickness workflows | none | critical | no move proposed |
| `scripts/lib/analytic_branch_tracking.py` | branch identity | canonical branch-selection helper | plots/shapes/tests | none | critical | no move proposed |
| `scripts/lib/branch_informed_spectrum_continuation.py` | completeness continuation | canonical targeted K10 provenance | gateway/Step-3A | none | critical | no move proposed |
| `scripts/lib/general_spectrum_completeness.py` | general audit machinery | negative readiness and completeness provenance | tests/audits | none | high | retain even if status becomes historical |
| `scripts/lib/variable_length_timoshenko.py` | shared Timoshenko implementation | 16 importers; eta-aware model | maps/audits/tests | none | critical | no move proposed |
| `src/my_project/fem/python_fem.py` | production baseline FEM | transform convention and 19 importers | comparisons/tests | none | critical | no move proposed |
| `tests/` tracked files | regression gates | canonical mathematical/branch/FEM checks | all refactors | none | critical | no move proposed |
| `docs/thickness_mismatch/eb_safe_prefix_stage_closure.md` | scientific closure | final safe-prefix status | navigation/reports | none | high | no move proposed |
| `results/eb_epsilon_lower_envelope_step3a/` | S3 provenance | counterexamples and verification audits | closure/exact phase | none | critical evidence loss | keep local backup and canonical report |
| `results/eb_rule_ab_exact_pareto/` | exact A/B/S provenance | thresholds/equivalence/locked validation | closure/cost proposal | none | critical evidence loss | keep local backup and canonical report |
| `results/eb_rule_s_cost_break_even/` | negative cost provenance | closes production-selector path | closure | none | critical evidence loss | keep local backup and canonical report |

## B. Soft archive candidates

These files stay in place but may later receive `historical`,
`closed-research`, or `superseded` status.

| Path | Current role | Reason | Dependencies | Canonical replacement | Risk of archiving | Required checks before move |
| --- | --- | --- | --- | --- | --- | --- |
| `scripts/analysis/thickness_mismatch/audits/audit_eb_timo_general_spectrum_completeness.py` | general completeness audit | completed `not_ready_for_step3` stage | imported by tests/gateway-related code | branch-informed gateway for targeted Step 3, but not equivalent generally | high | retain code; status only; verify importers and negative result links |
| `scripts/analysis/thickness_mismatch/postprocess/analyze_eb_safe_prefix_certification.py` | historical Rules A–D folds | superseded for final decision | reused helpers and tests | exact A/B/S report for final phase | high | freeze historical outputs and ensure exact phase does not depend on hidden semantics |
| `scripts/analysis/thickness_mismatch/postprocess/analyze_eb_epsilon_apriori_pilot.py` | geometry-only pilot postprocessor | scientific line closed | exact phase consumes corrected pilot, not this conclusion alone | stage closure | medium | preserve pilot provenance and documentation links |
| `scripts/analysis/thickness_mismatch/audits/audit_eb_validity_vs_timoshenko_stage1.py` | Stage-1 source study | completed source generation | imported by later audits | corrected gateway/closure are later evidence, not code replacement | high | status only until all helper importers are removed |
| `scripts/analysis/thickness_mismatch/audits/audit_eb_validity_fixed_epsilon_geometry_scan.py` | fixed-epsilon source study | completed source generation | predictor helpers reused by exact phase/cost | no full replacement | high | status only; extract helpers only with regression tests |
| `scripts/analyze_target_descendants_beta15_r5.py` | veering assessment generator | assessment concluded | docs/veering data | tracked closure tables | medium | reproduce tracked CSVs and retain command provenance |
| `scripts/analyze_flat_mu_bending_energy.py` | localization postprocessor | assessment concluded | imported by other scripts | tracked veering evidence | high | importer audit and output-hash comparison |
| `scripts/analysis/audit_desc6_veering_localization.py` | one-case localization audit | closed research question | result-only workflow | veering assessment | medium | confirm no article dependency |
| `scripts/run/run_analytic_coupled_rods_vibration_shapes_beta15_mu06_eps0025_001_ru.py` | one-case compatibility preset | configuration can eventually become a preset | shape runner | parameterized shape CLI | medium | preserve command/output names through wrapper test |
| `scripts/analysis/solid_fem_coupled_equal_rods_beta15.py` | original beta-15 3D FEM | multi-beta implementation exists | docs/results reproduction | `solid_fem_coupled_equal_rods.py` | high | compare mesh, constraints, output schema, and beta-15 results |

## C. Hard archive candidates after dependency audit

Physical moves are not approved here. A future move must use `git mv` and keep
a wrapper or migration note where documented/imported paths remain public.

| Path | Current role | Reason | Dependencies | Canonical replacement | Risk of archiving | Required checks before move |
| --- | --- | --- | --- | --- | --- | --- |
| root-level presentation family `scripts/plot_*_four_radii_compare.py` and `scripts/plot_mu_sweep_radius_fixed_four_betas_analytic.py` | implementation behind `scripts/run/` wrappers | could live under a presentation/legacy namespace | several are imported; documented outputs | existing `scripts/run/` wrappers only after wrapper inversion | high | full import graph, CLI snapshot, output hashes, retain thin compatibility modules |
| `scripts/analysis/solid_fem_coupled_equal_rods_beta15.py` | historical one-beta implementation | likely subsumed by multi-beta implementation | direct documentation and result provenance | `solid_fem_coupled_equal_rods.py` | high | exact geometry/deck/result comparison and reproduction note |
| `scripts/analysis/plot_mode_shapes_eta_beta_scan.py` plus sorted companion | two parallel mode-selection CLIs | shared plotting/reconstruction code is duplicated | one file is imported; docs list both | future parameterized `selection=sorted|descendant` CLI (not yet implemented) | critical branch-identity risk | build successor first; regression-test branch IDs, signs, normalization, filenames, CSV schema |
| possible orphan `audit_timoshenko_modes456_visualization.py` | visualization audit | no detected importer or doc reference | adjacent shape-construction data may be implicit | none | high provenance uncertainty | manual author review, result-file search, output reproduction, then wrapper/migration note |

## D. Local results archive candidates

| Path | Current role | Reason | Dependencies | Canonical replacement | Risk of archiving | Required checks before move |
| --- | --- | --- | --- | --- | --- | --- |
| `results/_smoke/` | temporary smoke outputs | 1,094 files / 205,691,561 bytes; reproducible intent | tests/manual diagnostics may reference examples | none | medium | separate cache from unique reports; restore test; no active process |
| legacy/cache subtrees of `results/eb_epsilon_baseline_thresholds/` | preserved pre-correction cache | large historical cache alongside corrected source | completeness provenance | corrected versioned cache/report | high | identify immutable legacy directories and hash them separately |
| `results/eb_epsilon_apriori_pilot/` | original pilot | superseded for final evidence by corrected branch-informed pilot | exact phase metadata/provenance | corrected pilot plus comparison audits | high | keep manifest/comparison; verify no exact-phase input points to original only |
| closed Step-3A/exact/cost directories | canonical closure evidence | stable and no longer active computation | closure documentation | none | critical | archive only as a complete immutable bundle; leave tracked/local index and restore instructions |
| active map/mode-shape caches | generated accelerators | potentially regenerable | plot-only workflows | source CSVs/report | high | classify cache vs source rows per workflow; never bulk-move entire directory |

The active 814 MB 3D-validation subtree is not a local archive recommendation;
it remains manual review because independent regeneration may require external
Gmsh/CalculiX inputs and substantial cost.

## E. Unknown/manual review

| Path | Current role | Reason | Dependencies | Canonical replacement | Risk of archiving | Required checks before move |
| --- | --- | --- | --- | --- | --- | --- |
| `private/` | private local workspace | contents intentionally ignored; ownership not inferred | possible article/source inputs | none | critical confidentiality/data loss | content-owner review; local-only encrypted/controlled backup policy |
| absent `paper_dorofeev_style/` and thickness article workspaces | external/ignored article workspaces referenced by docs | not present, so status cannot be verified | article-facing scripts | unknown | high | locate authoritative copies and compare paths before documentation changes |
| `.agents/` | empty local directory | purpose unclear | none detected | none | low | owner/tooling review |
| `results/eb_vs_timoshenko_longitudinal_suspect_modes/` | active-looking results without Markdown report | 25 files but no canonical report detected | suspect-mode audit code | none | high | identify generating command and decision record |
| possible orphan visualization audit | undocumented runnable | no detected path reference | adjacent audit modules | none | high | author/provenance review |

