# Residual Refactoring Audit

## Scope

This is the closing, read-only technical-debt audit after the first guarded code
refactor. It asks whether any additional abstraction is necessary before future
book work; similarity alone is not treated as a reason to unify code. No Python
source, test, script, formula, determinant, solver, FEM implementation, result,
private file, or cache was changed or executed for scientific computation.

The prompt named a second
`pre_refactor_inventory/duplicate_function_candidates.csv`, but that top-level
file does not exist. The actual and complete source used here is the tracked
[manifest copy](pre_refactor_inventory/manifests/duplicate_function_candidates.csv).
The review covers its 15 groups and eight supplemental families required by the
audit: compact number formatting, scalar grids, sign helpers, bisection/Brent
wrappers, SVD minima, local/global geometry, plotting utilities, and operation
counters.

Detailed evidence is recorded in:

- [residual refactor candidates](manifests/residual_refactor_candidates.csv);
- [possible-orphan audit](manifests/orphan_script_audit.csv);
- [manual path-reference review](manifests/manual_path_reference_review.csv).

## Repository state

- Repository root: `<repository-root>`.
- Branch: `refactor/common-thickness-scaling`.
- Initial HEAD: `6892bf9529e374250f69b77fbed5f7dc26b01164`.
- Initial worktree: exactly the four expected anisotropic-reservation document
  changes.
- Phase-0 reservation commit: `9682083` (`docs: reserve anisotropic rods
  research direction`).
- Baseline `e550a39cd94d9fb7f1711172693595c9f2eb1eea` is an ancestor of the
  current HEAD.
- Tag `pre-common-code-refactor-2026-07-29` still resolves to that baseline.
- The audit itself remains uncommitted; no push or pull request was performed.

The pre-audit integrity baseline was re-counted as 9,317 files under
`results/`, 52 under `private/`, 91 root-cache files, and 9,459 ignored files.
The closing check must reproduce those counts and show no tracked changes in
those locations.

## Previous completed refactors

- Repository backup, safety tag, dedicated refactor branch, and public
  [inventory](pre_refactor_inventory/README.md) are complete.
- Documentation and navigation were reorganized without moving scientific
  workflows or generated evidence.
- Mass-preserving thickness value algebra was centralized in
  [`thickness_scaling.py`](../../src/my_project/analytic/common/thickness_scaling.py).
  EB and Timoshenko validation, square-root backends, dataclasses, return types,
  exceptions, imports, and downstream scientific contracts remain separate.
- The anisotropic-rods direction has a reserved
  [documentation boundary](../anisotropic_rods/README.md), but no theory,
  scientific scope, API, scaffold, or implementation has been started.

## Duplicate-candidate review

Twenty-three groups were inspected from their actual implementations. The
review compared numerical formula, input domain, validation and exception
behavior, return type, multiplicity semantics, quality contract,
cache/provenance role, production or diagnostic status, downstream callers,
and expected evolution for every row. The
decision distribution is:

| Decision category | Groups | Count |
| --- | --- | ---: |
| `resolved-by-common-thickness-scaling` | mass-preserving tau value algebra | 1 |
| `safe-pure-utility-candidate` | `.10g` number token; `:g` filename token | 2 |
| `similar-but-different-contract` | tau dataclasses; local thickness; diameter/length; two CSV families; compact formatting; bisection/Brent wrappers | 7 |
| `scientifically-coupled-do-not-unify` | shape normalization; three root/spectral groups; mode CLI; circular section; SVD minima; geometry transforms; counters | 9 |
| `duplication-cheaper-than-abstraction` | general/gateway CLI definitions; determinant-local sign helpers | 2 |
| `historical-only` | old isolated-rod plotting pair | 1 |
| `manual-review` | scalar grid endpoint behavior | 1 |
| `already-shared` | none beyond the completed thickness extraction | 0 |

The completed thickness refactor removed only the common denominator
`1 + 2*mu*eta + eta**2` and the two tau divisions. The legacy EB
`ThicknessMismatchFactors` and Timoshenko `TauFactors` must not be collapsed:
their input validation, exception text, square-root implementation, dataclass
identity, public import path, and downstream geometry/matrix roles differ.

Two exact pure-utility clones are real:

1. `scripts/lib/diagnostic_common.py::number_token` and
   `scripts/lib/thickness_mismatch_diagnostic_helpers.py::number_token`;
2. `scripts/lib/analytic_branch_tracking.py::filename_number_token` and
   `scripts/lib/tracked_bending_descendant_shapes.py::filename_number_token`.

They are low-risk future candidates, but not worthwhile standalone tasks now.
Each has a tiny implementation, established importer paths, and filename
provenance implications; extracting it would create more churn than benefit.
If a future change already touches both copies, preserve the old imports as
wrappers and snapshot finite, exponent, signed-zero, plus-sign, and nonfinite
tokens together with representative output filenames.

The apparent pure scalar-grid duplicate is not confirmed. For
`start=0, stop=1, step=0.6`, `inclusive_grid` returns
`[0.0, 0.6, 1.0, 1.2]`, whereas the beta-grid helper clamps to
`[0.0, 0.6, 1.0]`; they agree only when the step lands on the endpoint. This is
a contract question, not a refactor opportunity, and remains manual review.

Other similarities are deliberately retained:

- local epsilon and local thickness-to-length helpers return quantities with
  different dimensions and structures;
- equal-rod diameter/length is a special-case scalar, while the mismatch
  helper is eta-aware and per-arm;
- CSV writers differ in field discovery, empty-row behavior, scientific
  formatting, newline bytes, JSON/bool/nonfinite encoding, and provenance;
- compact formatting uses `.12g` and configurable nonfinite text in one path,
  versus `.16e` and typed safe-prefix serialization in another;
- the general audit and branch gateway CLIs share option words but own
  different cache, isolation, validation, and algorithm contracts;
- a circular-section calculation in two models is shared physics, not a safe
  common properties API.

## Root and spectral helpers that must remain separate

The following are high- or critical-risk scientific contracts, even where
their control flow contains scanning, sign tests, bisection, deduplication, or
SVD:

| Family | Why it remains separate | Risk |
| --- | --- | --- |
| EB versus Timoshenko root scans | fixed manual bisection/list contract versus growing search range, nonfinite handling, Brent warnings, ndarray and NaN-fill contract | critical |
| Out-of-plane scan | different determinant family, validation, deduplication tolerance, and warning/completeness semantics | critical |
| Root deduplication | simple positive-finite values versus source-preserving factorized candidates, minimum singular values, and cross-family multiplicity | critical |
| Shape normalization | model-selected scale and returned metadata versus caller-provided scale/floor after sign normalization | high |
| SVD minima | raw determinant mode recovery, row-normalized completeness acceptance, block-scaled factorized multiplicity, and canonical branch/block SVD answer different questions | critical |
| Mode-shape CLI | descendant branch identity versus fixed current sorted index; unification could silently change branch semantics | critical |
| Local/global transforms | display-field basis and transverse sign versus production FEM DOF, stiffness, mass, and energy mapping | critical |
| Operation counters | different primitives, cache scopes, force-verification stages, and benchmark CSV meanings | high |

The scalar `sgn` expression is identical in several determinant-local bodies,
but extracting it would edit frozen scientific code for no behavioral or
maintenance gain. No common eigenvalue framework, shared root scan, shared SVD
fallback, or generic scientific counter is justified.

## Possible orphan script

The audited path is
[`audit_timoshenko_modes456_visualization.py`](../../scripts/analysis/thickness_mismatch/audits/audit_timoshenko_modes456_visualization.py).
It exists, is tracked, has a `__main__` guard, and has no conventional importer;
the diagnostic import-smoke test does dynamically import it. `git log --follow`
shows that its 843 lines were introduced together in commit
`1daf01d6f90a9d615351f5b3485f3dc713c38357`; blame shows no later line-level
history.

It is not a harmless visualization wrapper. It solves sorted Timoshenko roots
and reconstructs modes, so it was not run. Its unique contribution is a
pairwise comparison of full displacement fields against plotted centerlines
for modes 4--6, plus component and amplified-transverse figures. All five
declared normal outputs exist locally, as does the smoke report. The generated
report records a substantive diagnostic result: mean full-field similarity
about `0.1078` despite plotted-centerline similarity about `0.9943`.

Adjacent shape-construction, modes-4--6 diagnostic, and clean-plot workflows
are successors for presentation, but none replaces that quantitative
comparison. The script is therefore not an orphan and must not be deleted.
Recommended status: `soft-archive`, retained at its current path as unique
scientific provenance and not extended with new functionality. It is not
article-facing in tracked navigation. A future hard-archive review would first
need to check external citations to its exact output filenames.

## Soft archive consistency

The status sources are broadly consistent with the preservation-first
[archive policy](../archive_policy.md). Remaining ambiguity is navigational;
physical moves would add path and provenance risk without benefit.

| Workflow | Current visible status | Recommended status | Successor/report | Documentation action later |
| --- | --- | --- | --- | --- |
| Stage-1 applicability source study | `completed` | completed / soft archive | Stage-1 closure reports | Keep completion explicit; no new features. |
| Fixed-epsilon source study | `completed` | completed / soft archive | fixed-epsilon reports | Keep as reproduction path. |
| Original epsilon pilot and postprocessor | historical / closed | historical soft archive | later epsilon workflows | Preserve commands and outputs. |
| General completeness audit | `not_ready`, superseded for the targeted gate | superseded / soft archive | branch gateway and exact later gates | Retain negative-result provenance. |
| Historical Rules A--D postprocessor | superseded | superseded / soft archive | exact A/B/S workflow | Do not extend. |
| Step-3A baseline, exact A/B/S, Rule-S cost benchmark | completed / closed | completed provenance / soft archive | safe-prefix closure reports | Preserve frozen samples and negative cost result. |
| Counterexample dimensional plot | `completed` in script status; `active-diagnostic` in results index | completed provenance | counterexample report | Later harmonize `docs/results_index.md`; no move. |
| Focused veering generators, energy and localization audits | direction active, finite studies completed | completed provenance within active direction | veering status and figure index | Distinguish closed workflow from open research direction. |
| Original beta-15 solid FEM workflow | historical | historical soft archive | multi-beta solid FEM workflow | Keep old path for reproducibility. |
| Old `plot_freq_mu_vs_fem.py` one-off | legacy/unclear | superseded / soft archive after status review | maintained analytic/FEM comparison | Add a successor note only. |
| `run_fem_case.py`, `run_program1.py`, `run_program2.py` | compatibility/placeholder, manual review | retain historical / manual review | maintained FEM entry points | Do not create a hard archive without ownership review. |
| One-case vibration wrapper | compatibility wrapper | retain compatibility path | multi-case maintained runner | Do not add features; preserve CLI. |
| Modes-4--6 visualization audit | possible orphan / manual review | soft archive | adjacent shape diagnostics; unique report retained | Mark unique provenance and computational nature. |
| Older isolated-rod plot pair | historical-looking flat scripts | historical soft archive | newer parameterized reference workflows | Stop extending; no plotting framework extraction. |

## Manual path-reference review

All 26 non-`results/` references flagged for contextual review were classified;
none points into `private/` and no missing path was created.

| Classification | Count | Disposition |
| --- | ---: | --- |
| `historical path` | 15 | Preserve changelog/proposal/source provenance. |
| `optional external prerequisite` | 6 | Clarify ownership only when external article work resumes. |
| `manual unresolved` | 2 | Locate the dissertation source; clarify conceptual `scripts/reports`. |
| `superseded path` | 1 | Use the canonical tracked-descendant runner in new navigation. |
| `typo` | 1 | Later correct the veering dataset relative link to `../data/dataset_index.md`. |
| `tracked existing path` | 1 | The Step-3A manifest exists under the audits data subtree. |

The detailed manifest distinguishes the literal scanner resolution from known
tracked replacements. Automatic link rewriting is rejected because changelog
history, external workspaces, proposals, and contextual shorthand have
different meanings. The forbidden source documents remain unchanged; future
documentation corrections are recommendations only.

## src versus scripts/lib boundary

Location under `scripts/lib/` is not itself architectural debt. Production
equations, solvers, FEM, and the stable shared value algebra belong in `src/`;
diagnostic quality machinery, research-specific branch workflows, and figure
geometry can remain script-owned until a stable public contract actually
emerges.

| Module | Current location | Scientific or non-scientific | Importers | Production or diagnostic | Stable contract | Move risk | Recommended action |
| --- | --- | --- | ---: | --- | --- | --- | --- |
| thickness mismatch formulas | `src/my_project/analytic/formulas_thickness_mismatch.py` | scientific | 50 | production | yes | high | location acceptable |
| analytic solvers | `src/my_project/analytic/solvers.py` | scientific | 24 | production | yes | critical | must remain separate |
| Python FEM | `src/my_project/fem/python_fem.py` | scientific | 19 | production | yes | critical | location acceptable |
| common thickness algebra | `src/my_project/analytic/common/thickness_scaling.py` | scientific value algebra | 2 direct | production shared core | narrow and tested | high if expanded | location acceptable |
| variable-length Timoshenko | `scripts/lib/variable_length_timoshenko.py` | scientific | 16 | diagnostic model | established, not public-core | critical | move only if future work promotes a stable model |
| analytic branch tracking | `scripts/lib/analytic_branch_tracking.py` | scientific | 10 | canonical diagnostic selection | stable identity contract | critical | move only if touched by future public-package work |
| general completeness | `scripts/lib/general_spectrum_completeness.py` | scientific quality contract | 9 | diagnostic | workflow-stable | critical | must remain separate |
| branch continuation | `scripts/lib/branch_informed_spectrum_continuation.py` | scientific quality contract | 7 | diagnostic | workflow-stable | critical | must remain separate |
| straight factorized roots | `scripts/lib/straight_factorized_roots.py` | scientific | 3 | diagnostic/completeness | model-specific | critical | must remain separate |
| analytic shapes and MAC tracking | `scripts/lib/analytic_shapes.py`; `scripts/lib/mac_branch_tracking.py` | scientific | 16 each | diagnostic infrastructure | established | high | location acceptable |
| in-plane display geometry | `scripts/lib/in_plane_shape_geometry.py` | scientific presentation | 15 | diagnostic/display | established sign contract | critical | must remain separate from FEM |
| diagnostic common utilities | `scripts/lib/diagnostic_common.py` | non-scientific | 3 | diagnostic | small/local | low | location acceptable |

There is no pre-book architecture blocker and no case for a mass migration.

## Candidate next refactors

No next code refactor is recommended. The two confirmed filename-token clones
are recorded as `future-refactor-candidate` with timing `only-if-touched`, not
as mandatory work. Their measurable benefit is only one implementation point
if both call sites are already changing; otherwise a new shared API and wrappers
would cost more review than the duplication.

If either candidate is reopened, stop immediately if byte-for-byte filename
snapshots, public import wrappers, signed-zero/nonfinite behavior, or existing
output provenance cannot be preserved. It must be its own small commit after
pre/post snapshots and focused utility tests. It must not share a commit with
scientific code or anisotropic preparation.

## Refactors explicitly rejected

- A common eigenvalue-solver framework or unification of all root scans.
- A common SVD-minimum, multiplicity, completeness, or branch-tracking layer.
- Moving all `scripts/lib/` modules into `src/`.
- Mass-moving flat scripts or physically archiving closed workflows.
- Hard-archiving closed research when soft status and successor navigation are
  sufficient.
- Deleting negative results, historical reproduction paths, or the possible
  orphan script.
- Creating `RodEffectiveProperties`, an anisotropic API, or an anisotropic code
  scaffold before the book and a separately approved scientific model.
- Renaming `my_project`.
- Extracting sign, plotting, CLI, CSV, or counter abstractions merely to reduce
  line count.

## Final decision

`refactoring_complete_for_current_scope`

The remaining similarities either have different scientific, numerical,
quality, cache, provenance, or serialization contracts, or are too small to
justify churn. Soft archive classification is sufficient; no hard archive or
path move is needed. The repository is ready to stop technical refactoring
before future book reading. New common APIs at this point would be speculative.

No EB or Timoshenko root calculation, anisotropic calculation, FEM calculation,
continuation, mode reconstruction, or cost benchmark was run. The electronic
book was not inspected, searched, inferred, or analyzed.

## Conditions for reopening technical refactoring

Reopen only when at least one of these evidence-based conditions holds:

- a future change already touches every copy of one of the two exact pure
  filename utilities and focused byte-level regression snapshots are available;
- a diagnostic module gains a documented stable public contract and multiple
  production importers, making a narrow `src/` migration measurable;
- a concrete defect demonstrates inconsistent behavior between intended-equivalent
  contracts, with pre/post scientific gates defined before editing;
- book-derived scientific requirements are explicitly approved after reading,
  including physical meaning, notation, isotropic limits, and regression gates.

Code resemblance, directory preference, line-count reduction, or anticipated
anisotropy is not sufficient to reopen this phase.
