# Research Directions and Status

This index is the public map of the repository's scientific directions. It
describes the evidence visible in the tracked checkout; generated outputs and
local article workspaces may be absent from a fresh clone.

## Status vocabulary

- `stable-baseline`: verified foundation used by current workflows.
- `active-research`: an open scientific question with ongoing planned work.
- `active-diagnostic`: implemented diagnostic work that is not a final model
  or article-level conclusion.
- `completed`: the stated finite study or verification stage is complete.
- `closed-negative-result`: the investigated path answered its engineering or
  scientific question negatively and should not be extended without a
  material change in scope.
- `historical`: retained to reproduce the evidence chain, not as a preferred
  current workflow.
- `superseded`: a newer workflow or report is canonical, while the older path
  remains for provenance and compatibility.
- `planned`: scoped only as future work; no implementation is implied.
- `manual-review`: local state, provenance, or ownership is not sufficiently
  established for a stronger public status.

## Research directions

| Research direction | Status | Main question | Canonical documentation | Main implementation | Current conclusion |
| --- | --- | --- | --- | --- | --- |
| Baseline isotropic coupled beams | `stable-baseline` | In-plane Euler--Bernoulli frequencies of two rigidly coupled circular rods; influence of `beta`, `mu`, and `epsilon`; comparison with the baseline FEM and single-rod references | [equations](theory/equations.tex), [assumptions](theory/assumptions.md), [project rules](project_rules.md) | [`formulas.py`](../src/my_project/analytic/formulas.py), [`solvers.py`](../src/my_project/analytic/solvers.py), [`python_fem.py`](../src/my_project/fem/python_fem.py) | The determinant, signs, unknown order, normalization, and FEM transform convention are frozen baselines. |
| Branch identity and spectral tracking | `stable-baseline` | Preserve descendant identity separately from the branch's current sorted position, including close roots and MAC-based continuation | [project rules](project_rules.md), [script status](../scripts/STATUS.md) | [`analytic_branch_tracking.py`](../scripts/lib/analytic_branch_tracking.py), [`branch_informed_spectrum_continuation.py`](../scripts/lib/branch_informed_spectrum_continuation.py) | `branch_id` is continuation identity; `current_sorted_index` is metadata. Low-MAC assignments are not canonical without an accepted diagnostic. |
| Veering, quasi-degeneracy, modal exchange, and localization | `active-diagnostic` | Determine whether close spectral interactions support strict veering or only slower modal-character reorganization/localization | [terminology](veering/terminology.md), [strict assessment](veering/strict_veering_assessment.md), [slow-evolution assessment](veering/mu_slow_evolution_assessment.md) | tracked-branch, shape-MAC, and arm-energy workflows listed in [script status](../scripts/STATUS.md) | Strict claims require tracked branches, a local paired gap, non-crossing evidence, and mode-shape/MAC evidence. Present conclusions are deliberately cautious. |
| Thickness-mismatch model and mass-preserving radii | `active-diagnostic` | Extend the equal-radius geometry with `eta` while preserving total mass and the isotropic baseline limit | [model note](thickness_mismatch/README.md) | [`formulas_thickness_mismatch.py`](../src/my_project/analytic/formulas_thickness_mismatch.py), [`variable_length_timoshenko.py`](../scripts/lib/variable_length_timoshenko.py) | The `eta=0`, mass, swap, and selected sorted-root limits are checked. Branch/FEM validation for nonzero `eta` remains diagnostic. |
| Broad EB/Timoshenko and FEM comparison | `active-diagnostic` | Quantify where Euler--Bernoulli and Timoshenko spectra diverge and compare them with independent 1D/3D FEM evidence where applicable | [thickness-mismatch navigation](thickness_mismatch/README.md), [FEM status](thickness_mismatch/fem_validation_status.md), [frequency-map policy](thickness_mismatch/frequency_plot_computation_policy.md) | comparison maps and FEM audits in [script status](../scripts/STATUS.md) | Finite maps and validation cases exist, but they are not a universal applicability certificate or a full 3D validation of the ideal point joint. |
| `K=10` spectrum completeness and branch-informed gateway | `completed` | Establish reliable first-ten sorted roots with multiplicity and a right-hand completeness guard | [research plan](thickness_mismatch/eb_safe_spectrum_prefix_research_plan.md), [stage closure](thickness_mismatch/eb_safe_prefix_stage_closure.md) | general completeness and branch-informed gateway helpers/audits | Root 11 is the mandatory guard for roots 1--10. The targeted gateway resolved its declared dataset; this is finite numerical evidence, not a root-count theorem. |
| Geometry-only epsilon and Step 3A | `closed-negative-result` | Can the straight baseline or `epsilon_0` certify a safe EB spectrum prefix over the checked nonbaseline geometries? | [stage closure](thickness_mismatch/eb_safe_prefix_stage_closure.md) | Step-3A audit listed in [script status](../scripts/STATUS.md) | `epsilon_0` is not a certificate. `S3_12` and `S3_14` are confirmed finite-screen counterexamples; Step 3B was not needed to reject the checked lower-envelope hypothesis. |
| Exact Rules A/B/S | `completed` | Select finite-set EB-only prefix rules without observed false-safe cases on the declared partitions | [research plan](thickness_mismatch/eb_safe_spectrum_prefix_research_plan.md), [stage closure](thickness_mismatch/eb_safe_prefix_stage_closure.md) | exact A/B/S postprocessor listed in [script status](../scripts/STATUS.md) | Rule B degenerates to shear-only Rule S on all 49 checked geometries. The finite-sample safety result is not a continuous-domain guarantee. |
| Rule-S engineering selector | `closed-negative-result` | Does EB selection plus a Timoshenko suffix reduce cost relative to direct reliable Timoshenko `K=10`? | [stage closure](thickness_mismatch/eb_safe_prefix_stage_closure.md) | frozen Rule-S cost benchmark listed in [script status](../scripts/STATUS.md) | `rule_S_cost_not_beneficial`. This closes the current engineering-selector path; it does not mathematically refute Rule S. |
| Historical Rules A--D and pre-correction safe-prefix workflows | `historical` / `superseded` | Preserve the calibration and completeness evidence that led to the exact/branch-informed workflow | [historical plan](thickness_mismatch/eb_safe_spectrum_prefix_research_plan.md), [archive policy](archive_policy.md) | historical postprocessors and source-generation audits listed in [script status](../scripts/STATUS.md) | Retained for provenance and reproduction, not as current next steps. |
| Out-of-plane EB plus Saint-Venant torsion | `active-diagnostic` | Characterize the separate out-of-plane bending/torsion spectrum and compare it with an independent 1D continuum FEM | [theory note](theory/out_of_plane_eb_torsion.md) | out-of-plane solvers, maps, and FEM audit listed in [script status](../scripts/STATUS.md) | Determinant sanity and a 1D validation workflow exist. Generated outputs are not present in every checkout, and no full 3D validation claim is made. |
| Tracked article-promotion workflow | `active-research` | Promote reviewed diagnostics into external-facing figures and claims without contaminating canonical theory or diagnostic outputs | [article workflow](writing/article_workflow.md) | article-facing diagnostic scripts listed in [script status](../scripts/STATUS.md) | Promotion is an explicit review step; diagnostic output is not article evidence by default. |
| Local or historical article workspaces | `manual-review` | Locate the authoritative manuscript workspaces referenced by historical documentation | [refactoring status](refactoring/README.md) | no tracked workspace in this checkout | The referenced `paper_*` directories are absent from the tracked checkout. Their local ownership and current status must be resolved outside this index. |
| Anisotropic rods | `planned — deferred until refactoring completion` | Reserve a future direction without defining its scientific scope before refactoring and source study | [reserved direction](anisotropic_rods/README.md) | none | A documentation location is reserved. The physical model, source interpretation and implementation will be defined only after the current refactoring and after the user's electronic book has been studied and reproduced. |

When status records conflict, the newer stage-closure note or canonical report
takes priority. Mathematical source-of-truth priority remains the policy in
`AGENTS.md` and [project rules](project_rules.md).
