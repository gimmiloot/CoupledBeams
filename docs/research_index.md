# Research Directions and Status

This index is the public map of the repository's scientific directions. It
describes the evidence visible in the tracked checkout; generated outputs and
local article workspaces may be absent from a fresh clone.

## Production M-H/Timoshenko Lambda(beta) implementation checks (2026-10-07)

The [large geometry check](theory/mindlin_herrmann_timoshenko_lambda_beta_large_checks.md)
is `completed`, MHTIM_LAMBDA_BETA_LARGE_CHECKS=PASS. One fixed canonical
Lambda scale,3 length cases and3 thickness-contrast cases on37 angles,
with shared baseline:185 unique map cases +20 sparse swap controls.
All sorted12+guard13, beta0 direct/stepped references, baseline regression
and signed force/rotation/c/R checks pass.200 numerical bracket predictors,
5 seed scans, no fallback;68 curvature flags pass independent same-angle QR.
Two figures are implementation evidence, not tracked branches, physical
sensitivity/applicability or3D joint validation. Source/physics unchanged.

## Coupled longitudinal-theory screening (2026-10-07)

Subsequent [bounded thickness screening](theory/coupled_longitudinal_theory_thickness_screening.md)
is `completed`: scaling audit and all45 inventories PASS, fixed-case overlaps
and screening COMPLETE. h/h0=1/1.25/1.5/1.75/2, beta0/45/90 only, same
material/width/length/kappa. Max E/MH0.27081%, RL/MH0.47929%; all best
correspondences diagonal. Growth is not a universal h^2 spectral law;
end layers and changing sorted prefix matter. The close beta90,s_h1.25
positions9–10 pair is diagnostic context, not a veering claim. MH remains
a reference; no range extension/applicability threshold/tracking/FEM.

The [bounded screening](theory/coupled_longitudinal_theory_hierarchy_screening.md)
is `completed`: elementary and planar Love comparator gates PASS;
unchanged production MH references reused with hashes. Eight fixed G20
angles cover independently sorted12+guard13. Max model differences are
0.14103% E/MH and0.17514% RL/MH, all best geometric overlaps diagonal.
MH is a reference, not exact truth; cross-sectional contraction clamps
are not identical across theories. No energy classification, across-beta
tracking, applicability threshold, angle refinement or broader study.

## Single-rod M-H/Timoshenko source audit (2026-10-06)

Current [general-frame structural gate and fixed pilot](theory/mindlin_herrmann_timoshenko_rigid_joint.md#9-general-angle-geometry-действующий-project-contract):
**MHTIM_GENERAL_BETA_JOINT_GATE=PASS**. One common assembly recovers the
frozen beta0 gate exactly and passes project geometry, duality, rank8,
zero-limit/right-angle/swap/reflection checks. Only one G20 equal-arm
5/45/90 pilot is computed: exact energy-count inventories23/24/24 cover
first12+guard13 with full mode/resultant residuals. Reduced c-continuity
remains a model closure, not 3D elasticity. No angle map, hierarchy,
applicability claim or across-beta tracking follows from this result.

Subsequent [reduced rigid-joint beta0 gate](theory/mindlin_herrmann_timoshenko_rigid_joint.md):
**MHTIM_BETA0_JOINT_GATE=PASS**. Published common nodal contraction/rotation
closure is qualified as a variational reduced 1D model, not a direct finite
3D joint derivation. Three artificial splits reproduce the unchanged direct
G20 fixed--fixed rod: 7 MH/11 bending roots, combined12+guard13, c/R
transmission, mass-normalized shapes and reflection. Nonzero-angle spectra
and nonzero-angle validation were outside that beta0 stage; the subsequent
limited general-angle stage is separately recorded above.

Current separate decision: **PRODUCTION_MHTIM_FORMULATION_SELECTED**,
JANG_BARE_ISOTROPIC_REDUCED_MH_TIMOSHENKO, project rectangular K=5/6
(PRODUCTION_MHTIM_KAPPA_RESOLVED). The [finite single-rod gate](theory/mindlin_herrmann_timoshenko_single_rod.md#14-production-formulation-decision)
passes primary/independent roots, BC/energy and saturated min-max count
checks on the existing normalized G20 input. Numerical hierarchy passes;
HIERARCHY_SINGLE_ROD=PARTIAL_PASS records that reduced axial theories omit
the separately resolved contraction clamp. Fernandes now confirms the
alternative Ng factors; that normal block/preset is not production.
The single-rod stage selected no joint; the subsequent beta0-only assembly
above preserves its arm theory and qualifies the adopted reduced closure.
Earlier source-audit statuses below are historical answers to different questions.

The subsequent [rectangular prescription gate](theory/mindlin_herrmann_timoshenko_single_rod.md#13-production-rectangular-m-h-correction-prescription)
is `completed` as a finite source audit, with
`RECTANGULAR_MH_PRESET_MAPPING_UNRESOLVED`. Ng directly confirms the S1/S2
stiffness/inertia roles, but its normal Lamé coefficients differ from the
current reduced block. The second new PDF is Elishakoff--Tharu's circular
SSRN review, not the expected Fernandes 2022 full text, which was not found.
No production preset was adopted; this does not alter the earlier results.

The [canonical audit](theory/mindlin_herrmann_timoshenko_single_rod.md) is
`completed` as a diagnostic: `MHTIM_VARIANT_DEPENDENT`. Source energies,
local block structure, exact family mapping and low-frequency limits are
verified. Rucka's 100--120 kHz mode-count statements pass; figure agreement
is qualitative. Jang Fig.9(a) is conditional on explicit kappa (5/6 control),
whose numeric source value remains unstated. Published prescriptions are
`MH_SOURCE_VARIANTS_NOT_EQUIVALENT`; production coefficients remain
`PRODUCTION_MH_COEFFICIENTS_UNRESOLVED`. No angular joint or coupled-beam
implementation is implied. Bishop remains a standalone reference and its
closed combined-kinematics result is retained.

## Literature preparation (2026-10-05)

- [Mindlin–Herrmann + Timoshenko source map](literature/mindlin_herrmann_timoshenko_sources.md):
  literature preparation `completed`: three new local full texts (Rucka,
  Jang–Park–Lee, Liu 2021), extended reading of Banerjee 2019, source-specific
  factors and distinct reduced constitutive blocks documented. M-H axial +
  Timoshenko bending is the **current candidate** for the combined in-plane
  model; implementation/validation were pending at registration. The subsequent
  isolated diagnostic source audit is recorded above; production remains pending.
  Bishop is retained as a standalone reference/diagnostic theory, rather
  than the preferred production candidate for that combination. The prior
  [Timoshenko–Bishop audit](theory/timoshenko_bishop_single_rod.md) retains
  `COMBINED_KINEMATICS_NOT_UNIQUELY_DEFINED`; no hybrid closure or angular-joint
  conditions have been selected. Source validations apply to their own problems.
- [Longitudinal rod models and published control problems](literature/longitudinal_rod_models_sources.md):
  six sources in the initial registration. The subsequent narrow Rayleigh–Bishop
  [literature reproduction](theory/bishop_literature_reproduction.md) is
  `completed` with separate numerical/source-print statuses: internal checks
  pass, only one of five Marais frequencies matches print precision, and the
  rounded Popov table gives nu=.336842781730. Conditions for the project's
  angular joint remain unchosen.
- [Nonlinear in-plane and out-of-plane motions](literature/nonlinear_inplane_outofplane_sources.md):
  five local sources registered. Derivation of the project's nonlinear
  equations is deferred pending discussion with the supervisor.

The new M-H/Timoshenko registration is documentation only; the earlier
longitudinal reproduction and kinematic audit remain separate stages. Existing EB/RLB,
damping, anisotropic-rod and other research statuses below are unchanged.

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
| General-angle M-H/Timoshenko reduced joint | `completed` finite diagnostic | Does the project angle realization preserve duality/symmetries and a bounded finite inventory? | [general joint gate](theory/mindlin_herrmann_timoshenko_rigid_joint.md#9-general-angle-geometry-действующий-project-contract) | [general CLI](../scripts/analysis/verify_mindlin_herrmann_timoshenko_general_beta_joint.py), [common helper](../scripts/lib/mindlin_herrmann_timoshenko_joint.py) | All10 gates PASS; exact beta0, rank8, small angles/right angle and swap/reflection; fixed5/45/90 pilot12+guard13. No 3D joint proof, applicability/hierarchy study, sweep or tracking. |
| M-H/Timoshenko reduced rigid joint: beta0 transparency | `completed` finite diagnostic | Does published common-DOF assembly preserve one homogeneous direct rod? | [joint theory and gate](theory/mindlin_herrmann_timoshenko_rigid_joint.md) | [beta0 CLI](../scripts/analysis/verify_mindlin_herrmann_timoshenko_beta0_joint.py), [joint helper](../scripts/lib/mindlin_herrmann_timoshenko_joint.py) | Seven separate gates PASS; splits .5/.35/.65, 12+guard13 and full c/R form checks. Reduced closure adopted, not derived from 3D elasticity; beta!=0 not computed or validated. |
| Single rectangular rod: M-H axial + Timoshenko bending | `completed` finite diagnostic; production Jang closure/K selected | Verify a finite CC rod and elementary/planar RL/MH hierarchy | [source audit and finite gate](theory/mindlin_herrmann_timoshenko_single_rod.md) | [finite CLI](../scripts/analysis/verify_mindlin_herrmann_timoshenko_single_rod.py), [unchanged source CLI](../scripts/analysis/reproduce_mindlin_herrmann_timoshenko_literature.py), [M-H helper](../scripts/lib/mindlin_herrmann_longitudinal.py) | Finite MH spectrum PASS; hierarchy PARTIAL_PASS for contraction-clamp interpretation. Source variants remain different, source Jang kappa unstated. Ng/Fernandes is alternative only; no M-H coupled rods/joint selection. |
| Rotational-joint continuation: circular EB reviewer diagnostic | `completed` finite diagnostic | Do six low modes approach exact RIGID as rotational stiffness grows? | [scientific note and qualifications](laminated_beams/circular_eb_rotational_spring_rigid_limit.md) | [existing-kernel/generic-solver orchestration](../scripts/analysis/joint_review/check_circular_eb_spring_spectrum.py) | Six local descendants confirmed; rigid equivalence only for 6+guard 7. Positions 11–12 and seed-06 nonmonotonic rotation remain qualified. Transmission/equilibrium and thin-joint asymptotic justification remain an open theoretical question; no real-joint calibration. |
| Single rectangular rod: Timoshenko + Rayleigh–Bishop kinematics | `completed` finite audit; diagnostic-only | Does a common displacement/energy field yield the two unchanged linear subsystems? | [kinematics and energy audit](theory/timoshenko_bishop_single_rod.md) | [exact algebra CLI](../scripts/analysis/audit_timoshenko_bishop_single_rod.py), [targeted tests](../tests/test_timoshenko_bishop_single_rod.py) | `COMBINED_KINEMATICS_NOT_UNIQUELY_DEFINED`: centered cross coefficients vanish in two candidates, but raw self terms fail the unchanged Timoshenko limit. A relaxed hybrid needs explicit extra closure; no combined spectrum or joint conditions. |
| Local longitudinal Rayleigh–Bishop literature controls | `completed` finite study; diagnostic-only | Reproduce Marais section 4 and Popov–Sadovsky (5),(6),(9)–(15) | [canonical report](theory/bishop_literature_reproduction.md) | [bounded CLI](../scripts/analysis/reproduce_bishop_literature.py), [module](../scripts/lib/bishop_longitudinal.py) | Source audit, energy/orthogonality, independent 50/70 dps pass. Marais print match 1/5; Popov rounded-table ranking Rayleigh–Love/Bishop/wave, with figure/reference qualifications. No angular joint or rectangular-system conclusions. |
| Representative complex EB/RLB damping confirmation | `completed`; `PAUSED_FOR_SUPERVISOR_DIRECTION` | Does elastic G-ratio predict actual weak-damping zeta-ratio? | [D22/K23 technical report](laminated_beams/inplane_kelvin_voigt_eb_rlb_complex_confirmation.md); [supervisor synthesis](laminated_beams/inplane_kelvin_voigt_research_status_for_supervisor.md) | [fixed six-target orchestration](../scripts/analysis/laminated_beams/confirm_inplane_kelvin_voigt_eb_rlb.py) | Six new roots at d=.001, four reused states; five ratios agree within .001904%. R3 contrast and positive R4 correction confirmed. Further direction awaits discussion. |
| Sparse elastic RLB damping participation vs EB | `completed` — `SCREENING_AND_COMPARISON_COMPLETED` | Compare changes of frequency and rotational damping predictor at four fixed angles | [D21/K22 screening](laminated_beams/inplane_kelvin_voigt_rlb_elastic_screening.md) | [bounded reduced screening/comparison](../scripts/analysis/laminated_beams/screen_inplane_kelvin_voigt_rlb_elastic.py) | 24 RLB states, 24 same-angle form matches to K15. At 75 degrees, pair 05 changes Omega by −2.075% and G by −74.49%; 12 exact inactive states. No new positive-d roots or cross-beta tracking. |
| RLB rotational KV production routing | `completed` technical validation — `RLB_KV_PRODUCTION_PASS` | Extend the existing reduced/full dispatcher without starting a new physical study | [D20/K21 architecture](laminated_beams/inplane_kelvin_voigt_rlb_solver_architecture.md) | [shared dispatcher](../scripts/lib/inplane_kelvin_voigt_solver.py) | K12 RLB active/inactive controls, full/reduced equivalence, exact invS=J=0 EB limit and one unequal-arm elastic K11 control pass. EB regressions preserved; no sparse RLB screening or asymmetric viscous roots. |
| Local weak-damping parity of a KV eigenbranch | `completed` theoretical note | Explain odd decay and even frequency without equal-arm symmetry | [D19/K20 parity note](laminated_beams/inplane_kelvin_voigt_weak_damping_parity.md) | Analytic proof from the existing full boundary matrix; no new implementation | A simple isolated branch obeys z(-d)=-conj(z(d)); a/d and zeta/d corrections and frequency shift are O(d²). Covers real conservative EB/project RLB with unequal arms; K19 is an illustration, no new calculations or strong-damping claim. |
| Completion of three-state weak-damping KV comparison | `completed` finite study | Compare elastic damping predictors with six physical states | [D18/K19 completion](laminated_beams/inplane_kelvin_voigt_targeted_weak_damping_completion.md) | [two-target reduced continuation](../scripts/analysis/laminated_beams/complete_inplane_kelvin_voigt_weak_damping.py) | Only A/C .005 newly computed; four rows reused. Ranking A>B>C confirmed, largest relative predictor departure at .005 belongs to C (.3811%). Full C raw rank qualification retained; no solver changes or further parameters. |
| EB KV solver routing | `completed` technical regression | Route exactly identical arms through eta± blocks and preserve the unequal-arm full solver | [D17/K18 architecture](laminated_beams/inplane_kelvin_voigt_solver_architecture.md) | [production dispatcher](../scripts/lib/inplane_kelvin_voigt_solver.py) | Seven saved controls, 12 path evaluations; reduced passes, full C retains a raw rank qualification. A/C full physical residuals now pass unchanged gates. No new physical points. |
| KV A/C symmetry diagnostic | `completed` bounded diagnostic | Explain the saved A/C physical failures and C's second-singular flag | [D16/K17 report](laminated_beams/inplane_kelvin_voigt_ac_diagnostics.md) | [same-target half-system diagnostic](../scripts/analysis/laminated_beams/diagnose_inplane_kelvin_voigt_ac.py) | Consistent analytic half forms pass unchanged physical gates at the saved roots; transfer/reaction conditioning identified. K16 and full singular flag retained; no beta trigger or .005 continuation. |
| Three targeted weak-damping KV states | `completed` bounded attempt — `PARTIAL_NUMERICAL_QUALIFICATIONS` | Compare full complex roots with three K15 ACTIVE elastic predictors at d=.001/.005 | [six-target report](laminated_beams/inplane_kelvin_voigt_targeted_weak_damping.md) | [bounded K12 solver orchestration](../scripts/analysis/laminated_beams/check_inplane_kelvin_voigt_targeted_weak_damping.py) | D15/K16: two INTERMEDIATE roots accepted; STRONG/WEAK first targets retain physical gate failures after two attempts, second targets not run. No extra parameters or whole-model validation. |
| Rotational KV elastic screening | `completed` bounded diagnostic | Identify structural inactivity and first-order damping predictors from elastic forms | [24-state EB screening](laminated_beams/inplane_kelvin_voigt_elastic_screening.md) | [real EB screening entry point](../scripts/analysis/laminated_beams/screen_inplane_kelvin_voigt_elastic.py) | D14/K15: H/L/L/H, mu=0, kappa=1, beta=0/5/45/75; 12 inactive and 12 active states. K12 read-only slope control agrees; no new positive-d roots or cross-beta identity. |
| External KV literature benchmarks | `completed` bounded diagnostic — `PASS_WITH_SOURCE_PRINT_QUALIFICATIONS` | Reproduce Failla Table1 and Hong Tables2–3 before expanding the rotational KV study | [source mappings and chronological comparison](laminated_beams/inplane_kelvin_voigt_literature_benchmarks.md) | [source-specific benchmark](../scripts/analysis/laminated_beams/benchmark_inplane_kelvin_voigt_literature.py) | D13 separately accepts Table2 equation-level consistency and obtains five Table3 PRINT_MATCH + SOLVER_PASS roots. K13 PARTIAL and Failla/Table2 printed mismatches remain. No whole-model validation or new parameter study. |
| Baseline isotropic coupled beams | `stable-baseline` | In-plane Euler--Bernoulli frequencies of two rigidly coupled circular rods; influence of `beta`, `mu`, and `epsilon`; comparison with the baseline FEM and single-rod references | [equations](theory/equations.tex), [assumptions](theory/assumptions.md), [project rules](project_rules.md) | [`formulas.py`](../src/my_project/analytic/formulas.py), [`solvers.py`](../src/my_project/analytic/solvers.py), [`python_fem.py`](../src/my_project/fem/python_fem.py) | The determinant, signs, unknown order, normalization, and FEM transform convention are frozen baselines. |
| Branch identity and spectral tracking | `stable-baseline` | Preserve descendant identity separately from the branch's current sorted position, including close roots and MAC-based continuation | [project rules](project_rules.md), [script status](../scripts/STATUS.md) | [`analytic_branch_tracking.py`](../scripts/lib/analytic_branch_tracking.py), [`branch_informed_spectrum_continuation.py`](../scripts/lib/branch_informed_spectrum_continuation.py) | `branch_id` is continuation identity; `current_sorted_index` is metadata. Low-MAC assignments are not canonical without an accepted diagnostic. |
| Frequency-map computation policy | `stable-baseline` | How ordinary frequency maps are computed without conflating plot production with research-grade spectral certification | [project-wide policy](numerics/frequency_map_computation_policy.md) | Model-specific existing solvers and local policy instances | Ordinary maps use `fast_plot`; `certified_audit` is triggered only by an explicit scientific purpose or an unresolved numerical event. |
| Veering, quasi-degeneracy, modal exchange, and localization | `active-diagnostic` | Determine whether close spectral interactions support strict veering or only slower modal-character reorganization/localization | [terminology](veering/terminology.md), [strict assessment](veering/strict_veering_assessment.md), [slow-evolution assessment](veering/mu_slow_evolution_assessment.md) | tracked-branch, shape-MAC, and arm-energy workflows listed in [script status](../scripts/STATUS.md) | Strict claims require tracked branches, a local paired gap, non-crossing evidence, and mode-shape/MAC evidence. Present conclusions are deliberately cautious. |
| Thickness-mismatch model and mass-preserving radii | `active-diagnostic` | Extend the equal-radius geometry with `eta` while preserving total mass and the isotropic baseline limit | [model note](thickness_mismatch/README.md) | [`formulas_thickness_mismatch.py`](../src/my_project/analytic/formulas_thickness_mismatch.py), [`variable_length_timoshenko.py`](../scripts/lib/variable_length_timoshenko.py) | The `eta=0`, mass, swap, and selected sorted-root limits are checked. Branch/FEM validation for nonzero `eta` remains diagnostic. |
| Broad EB/Timoshenko and FEM comparison | `active-diagnostic` | Quantify where Euler--Bernoulli and Timoshenko spectra diverge and compare them with independent 1D/3D FEM evidence where applicable | [thickness-mismatch navigation](thickness_mismatch/README.md), [FEM status](thickness_mismatch/fem_validation_status.md), [frequency-map policy](numerics/frequency_map_computation_policy.md) | comparison maps and FEM audits in [script status](../scripts/STATUS.md) | Finite maps and validation cases exist, but they are not a universal applicability certificate or a full 3D validation of the ideal point joint. |
| `K=10` spectrum completeness and branch-informed gateway | `completed` | Establish reliable first-ten sorted roots with multiplicity and a right-hand completeness guard | [research plan](thickness_mismatch/eb_safe_spectrum_prefix_research_plan.md), [stage closure](thickness_mismatch/eb_safe_prefix_stage_closure.md) | general completeness and branch-informed gateway helpers/audits | Root 11 is the mandatory guard for roots 1--10. The targeted gateway resolved its declared dataset; this is finite numerical evidence, not a root-count theorem. |
| Geometry-only epsilon and Step 3A | `closed-negative-result` | Can the straight baseline or `epsilon_0` certify a safe EB spectrum prefix over the checked nonbaseline geometries? | [stage closure](thickness_mismatch/eb_safe_prefix_stage_closure.md) | Step-3A audit listed in [script status](../scripts/STATUS.md) | `epsilon_0` is not a certificate. `S3_12` and `S3_14` are confirmed finite-screen counterexamples; Step 3B was not needed to reject the checked lower-envelope hypothesis. |
| Exact Rules A/B/S | `completed` | Select finite-set EB-only prefix rules without observed false-safe cases on the declared partitions | [research plan](thickness_mismatch/eb_safe_spectrum_prefix_research_plan.md), [stage closure](thickness_mismatch/eb_safe_prefix_stage_closure.md) | exact A/B/S postprocessor listed in [script status](../scripts/STATUS.md) | Rule B degenerates to shear-only Rule S on all 49 checked geometries. The finite-sample safety result is not a continuous-domain guarantee. |
| Rule-S engineering selector | `closed-negative-result` | Does EB selection plus a Timoshenko suffix reduce cost relative to direct reliable Timoshenko `K=10`? | [stage closure](thickness_mismatch/eb_safe_prefix_stage_closure.md) | frozen Rule-S cost benchmark listed in [script status](../scripts/STATUS.md) | `rule_S_cost_not_beneficial`. This closes the current engineering-selector path; it does not mathematically refute Rule S. |
| Historical Rules A--D and pre-correction safe-prefix workflows | `historical` / `superseded` | Preserve the calibration and completeness evidence that led to the exact/branch-informed workflow | [historical plan](thickness_mismatch/eb_safe_spectrum_prefix_research_plan.md), [archive policy](archive_policy.md) | historical postprocessors and source-generation audits listed in [script status](../scripts/STATUS.md) | Retained for provenance and reproduction, not as current next steps. |
| Out-of-plane EB plus Saint-Venant torsion | `active-diagnostic` | Characterize the separate out-of-plane bending/torsion spectrum and compare it with an independent 1D continuum FEM | [theory note](theory/out_of_plane_eb_torsion.md) | out-of-plane solvers, maps, and FEM audit listed in [script status](../scripts/STATUS.md) | Determinant sanity and a 1D validation workflow exist. Generated outputs are not present in every checkout, and no full 3D validation claim is made. |
| Tracked article-promotion workflow | `active-research` | Promote reviewed diagnostics into external-facing figures and claims without contaminating canonical theory or diagnostic outputs | [article workflow](writing/article_workflow.md) | article-facing diagnostic scripts listed in [script status](../scripts/STATUS.md) | Promotion is an explicit review step; diagnostic output is not article evidence by default. |
| Local or historical article workspaces | `manual-review` | Locate the authoritative manuscript workspaces referenced by historical documentation | [refactoring status](refactoring/README.md) | no tracked workspace in this checkout | The referenced `paper_*` directories are absent from the tracked checkout. Their local ownership and current status must be resolved outside this index. |
| Anisotropic rods | `active-research` | Validate the first ideal rigid angular-joint model and its rectangular orthotropic EB endpoint after the completed Chapter-2 single-rod source gates | [direction status](anisotropic_rods/README.md), [notation translation](anisotropic_rods/yartsev_ch2_notation_translation.md), [rigid-joint gate](anisotropic_rods/yartsev_ch2_rigid_angular_joint.md), [rectangular EB validation](anisotropic_rods/yartsev_ch2_rectangular_eb_validation.md) | [`yartsev_ch2_monoclinic_rod.py`](../scripts/lib/yartsev_ch2_monoclinic_rod.py), [`yartsev_ch2_coupled_rods.py`](../scripts/lib/yartsev_ch2_coupled_rods.py), [`yartsev_ch2_rectangular_eb.py`](../scripts/lib/yartsev_ch2_rectangular_eb.py) | The rigid-joint pilot passed. The finite `theta=0` rectangular EB/exact/unequal-length/1D-FEM validation remains `PARTIAL_PASS`: proportional refinement closes the original first-three accuracy threshold with raw error `6.18e-6`, but the targeted status is `FAIL_CONVERGENCE_ORDER` because mode 1 reaches a dense-eigensolver conditioning floor and violates the unchanged monotonicity allowance. No model coefficient or threshold changed. This is not a final coupled model, stable baseline, off-axis study, production API, unequal-thickness study, or 3D validation. |

When status records conflict, the newer stage-closure note or canonical report
takes priority. Mathematical source-of-truth priority remains the policy in
`AGENTS.md` and [project rules](project_rules.md).
