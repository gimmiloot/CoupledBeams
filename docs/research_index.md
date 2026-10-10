# Research Directions and Status

This index is the public map of the repository's scientific directions. It
describes the evidence visible in the tracked checkout; generated outputs and
local article workspaces may be absent from a fresh clone.

## FEM-3C: bounded numerical robustness and dissertation verification (2026-10-10)

[Technical continuation](numerics/nlsp_nonlinear_dynamic_3d_fem_pilot.md#fem-3c-numerical-robustness-and-dissertation-verification)
and [scientific summary](numerics/nlsp_straight_rod_3d_fem_verification_summary.md)
keep one straight-rod task and all historical evidence. The four actual C1 jobs
pass preload/release/recovery and the preregistered temporal/spatial guide.
Dt=5.40090e-10 and Dh=2.10750e-8 are .00191202 and .07460934 of the fixed
2.824717e-7 baseline; interpolation comparability passes. The updated fine
1D/3D evolving-w difference is 7.76472% on its own scale, not an artificially
preserved 7.30%. These observed changes are not strict FEM error bounds.
Both full-period 1D trajectories complete, with all8 spatial PARTIAL (2/8).
Both medium full-period 3D jobs reach T1; nonlinear w difference is 9.45%,
evolving-correction difference 10.07% on their full-horizon scales. Absolute
model differences grow; quarter-period robustness is not full-period certification.
The limited straight-rod verification is complete with qualifications. [D16](memory/decisions.md#nlsp-d16) and [K17](memory/knowledge.md#nlsp-k17)
records the separate conditional authorization, preserved physics/load and
no-retry stop. This does not validate every nonlinear seven-field coupling,
angular joints or out-of-plane stability; native energy remains qualified.

## FEM-3B: evolution of NL-minus-L beyond the static initial offset (2026-10-09)

[Scoped continuation](numerics/nlsp_nonlinear_dynamic_3d_fem_pilot.md#fem-3b-nonlinear-correction-evolution-and-longer-horizon)
uses the same frozen preload/release problem. The old .05T1 result is diagnosed
read-only; one p64 and optional saved-state p48 trajectory to .5T1 guide a
predeclared choice of .25T1 or .5T1 before either new medium 3D job. Total
correction, evolving correction and the additional 1D initial-state/same-state
nonlinear decomposition are reported separately. Native energy bookkeeping
remains qualified; no independent 3D time/mesh certification is implied.
Both native jobs reach .25T1 with502 matching frames. Evolving w maxima are
3.58562e-6/3.86809e-6, with7.30261% max model difference; the signal exceeds
observed recovery/output differences. Energy and full1D spatial remain PARTIAL.
[Decision D15](memory/decisions.md#nlsp-d15) and [evidence K16](memory/knowledge.md#nlsp-k16)
preserve all earlier decisions and qualifications.

## FEM-3AR: short independent free-motion pilot completed with qualifications (2026-10-09)

[Controlled continuation](numerics/nlsp_nonlinear_dynamic_3d_fem_pilot.md#fem-3ar-controlled-continuation-after-native-elke-output-path-failure)
separately authorizes the corrected output route after the historical native
ELKE failure. Both medium preload+dynamic jobs and the retained p64 linear/NL
references reach .05T1. Actual preload repetition, instantaneous release and
restoring motion pass. Native energy remains PARTIAL owing to source-localized
bookkeeping/pseudo-time output semantics. The2.24247% sampled w-correction
difference mostly reflects the inherited static offsets, not certified evolution
of a dynamic nonlinear correction. No new physics/mesh/load or accuracy sweep;
[decision D14](memory/decisions.md#nlsp-d14), [evidence K15](memory/knowledge.md#nlsp-k15).

## FEM-3A: source/preflight ready, native execution blocked (2026-10-09)

[Dynamic pilot](numerics/nlsp_nonlinear_dynamic_3d_fem_pilot.md) preserves the
h=.10 statics/geometry/physics and documents installed2.22 sequential preload,
OP=NEW/zeroGRAV/STEP release, zero physical initial velocities and ALPHA=0.
Frozen p64 coordinates and released1D acceleration checks pass. Release/velocity
protocol statuses remain PARTIAL because actual3D evidence is absent. First linear
medium3D job terminated with access violation0xC0000005; no accepted static or
transient state can be verified. No NL job,1D trajectory or extra fixture/retry
ran. Overall BLOCKED_BY_SOLVER; not evidence against V0 or general absence of
CalculiX dynamic capability. Read-only source/binary tracing
localizes first LINEAR STATIC ELKE reading freed veold; a corrected output-only
preview remains NOT_RUN. Actual release/transfer/energy/comparison are NOT_RUN.
FEM-2R static/old strict/scoped statuses remain; stop without further execution.

## FEM-2R: completed static diagnostic with qualifications (2026-10-09)

[Controlled continuation](numerics/nlsp_nonlinear_static_3d_fem_validation.md#fem-2r-controlled-continuation-after-input-serialization-failure)
reuses the saved1D equilibria/load and three meshes after the separate explicit
execution authorization. All six new linear/NLGEOM static jobs reach full load,
complete output and independent reaction-balance gates. Refined Delta w=-9.16268e-6
versus1D-8.87008e-6; negative signs agree, full-profile difference3.19353% of the
common correction scale. The signal exceeds observed mesh/recovery/rounding
measures; these are not exact continuum/complete solver-error bounds. Overall
FEM2_DIAGNOSTIC_COMPLETE_WITH_QUALIFICATIONS, not PHYSICAL_NONLINEAR_VALIDATION_PASS
for all V0 coefficients. Historical failed attempt/strict PARTIAL remain.
No new mesh,1D/modal/ODE solve or FEM-3; stop after the bounded report.

## FEM-2: bounded nonlinear static comparison (2026-10-09)

[Static report](numerics/nlsp_nonlinear_static_3d_fem_validation.md) fixes the same
L=1,b=.20,h=.10 rod and a dead global transverse load selected before FEM.
1D p48/p64 linear and quartic nonlinear statics pass their equilibrium gates;
midspan w=.005 / .00499112992069312 and Delta w=-8.8700793068e-6.
The first medium linear 3D attempt failed while parsing a numeric *STATIC field,
before equilibrium. No medium NL/fine/refined jobs ran; 3D comparison and static
mesh convergence remain NOT_RUN, overall PARTIAL. This is a new input-formatting
issue, not physical evidence against V0. FEM-3 is not justified or authorized
by the incomplete comparison. Historical scope/statuses remain unchanged.

## FEM-1R: bounded extra mesh (2026-10-09)

[Continuation](numerics/nlsp_linear_rectangular_3d_fem_validation.md#fem-1r-one-additional-mesh-refinement)
loads immutable1D/three-grid evidence and adds exactly one .020 C3D10 mesh/24-mode
job. All8 shape identities and.1% preset mesh changes pass; frequencies decrease
regularly across four grids. Updated axial/bending/twist differences.47666%/
1.36648-1.74302%/4.40336-4.97569% retain model/end/warping qualifications. Linear
baseline supports a separately authorized limited FEM-2 test; no nonlinear V0
validation, new1D solve or automatic fifth grid/FEM-2 follows.

## Full-family thick rectangular FEM-1 (2026-10-08)

The [linear3D report](numerics/nlsp_linear_rectangular_3d_fem_validation.md)
uses h=.10,b=.20,L=1 chosen BEFORE FEM: first MH acoustic is position8. Eight1D
modes within omega3.4651 cover both bending planes, twist and axial contraction.
Three real C3D10 meshes each yield24 complete modes; all8 shapes match uniquely,
with no additional3D mode in the window. Geometry/1D/mesh/execution/identification
PASS; mesh convergence and all-family quantitative comparison PARTIAL (5/8 exceed
predeclared.1% medium/fine criterion). Observed differences about.49% axial,
1.40-1.81% bending,4.46-5.04% twist are qualified by mesh/end/warping effects.
No coefficient fitting, backup thickness, fourth mesh or nonlinear step followed.

## Local 3D FEM readiness audit (2026-10-08)

The [readiness report](numerics/nlsp_3d_fem_environment_readiness.md) confirms
Gmsh4.15.2 and already unpacked CalculiX2.22 x64 run locally, even though PATH/
GMSH_EXE/CCX_EXE are unset. Historical linear solid jobs are retained. Static and
direct transient NLGEOM support is confirmed in the local2.22 manual; project
nonlinear workflows and the new rectangular G20 test remain NEEDS_EXECUTION_TEST.
No mesh/job/install/model change occurred. One bounded linear rectangular test
is proposed, not executed or automatically authorized. NLSP physical sanity
remains DIAGNOSTIC_COMPLETE_WITH_QUALIFICATIONS; other scoped stops are unchanged.

## Bounded planar physical sanity (2026-10-08)

The [physical consistency note](theory/weakly_nonlinear_planar_physical_sanity_checks.md)
compares the saved p64 .05 one-T1 history with leading second-order profiles and
one new p64 .025 tight trajectory. Normalized motions show expected2/4 amplitude
patterns; leading-deviation reductions are approximately4. Exact sign symmetry,
formal classical stretching and retained-strain diagnostics pass. Reaction
accounting includes finite-p strong residuals and the explicit degree4 rotational
truncation term. Overall **DIAGNOSTIC_COMPLETE_WITH_QUALIFICATIONS**: no physical
validation, no independent half-case p/time certificate, strict/historical PARTIAL
retained. Physics/RHS/basis/BC unchanged; no further calculation is authorized.

## Prepared movement through one linear period (2026-10-08)

The [one-T1 continuation](theory/planar_prepared_initial_state.md#prepared-one-period-feasibility)
reuses exactly the saved common initial coefficients and time settings at p48/p64.
All3 runs reachT1=19.791162590151373 in EXPLORATORY_NOT_CERTIFIED mode, with no
new projection/MP/BVP/eigen/symbolic audit. Old0.1T1 prefixes match; temporal8/8
and energy/mass/safety pass. Spatial remains PARTIAL:7/8 pass, theta_t max relative
1.590e-4 exceeds1e-4. Absolute spatial maxima grow roughly1.94–2.73x; changes of
full-horizon characteristic denominators are recorded separately. The sampled
results support bounded practical computability, not continuous-PDE convergence,
a nonlinear periodic orbit or out-of-plane/Floquet stability. Physics/IC/BC/basis
and all historical strict/PARTIAL qualifications remain; no automatic extension.

## Historical prepared-state precision / short feasibility continuation (2026-10-08)

The [continuation evidence](theory/planar_prepared_initial_state.md#prepared-precision-feasibility)
separates strict verification from explicitly authorized exploratory execution.
Exact Gram/analytic moments and one initial-only endpoint-constrained L2 policy
preserve the frozen common physical state and full independent four-field space.
Both p48/p64 projections pass the original1e-6 gate. Independent MP45/70 checks
localize float64 relative strong/weak discrepancy to stored Gauss data; its
original2e-12 gate remains FAIL, so strict verification is PARTIAL.

**NLSP_PREPARED_FEASIBILITY_RUN=COMPLETED_EXPLORATORY_NOT_CERTIFIED.** All three
p48/p64 short controls reach0.1T1 with positive mass and accepted safety/energy.
Temporal comparison passes all8 components; spatial comparison remains PARTIAL:
7/8 pass, theta_t relative max1.405e-4 exceeds1e-4. Charged numerical work53.49s,
including17.34s of integration, zero new BVP/model eigensolves. This establishes
bounded finite-dimensional computability, not strict continuous-PDE accuracy.
Historical prepared/zero-u-c PARTIAL reports and LONG closed / KV paused /
angular same-clamp reference unavailable remain. No longer run or new stage
is selected. Use the explicit `--compute --feasibility` mode; default strict
workflow and cached report/plot behavior are preserved.

## Historical prepared planar initial-state gate (2026-10-08)

The [focused preparation note](theory/planar_prepared_initial_state.md)
separates constant/second-harmonic and free parts of the saved leading axial
response. Stat/harm profiles, physical derivatives and endpoint jets converge
on p64/p96; one common numerical reference produces a distinct O2 initial
state and a quintic theta correction with acceleration compatibility through
cubic order. Finite-amplitude higher-order residuals are retained.

**NLSP_PREPARED_INITIAL_STATE_PILOT=PARTIAL.** Neither permitted pair
p32/p48 or p48/p64 admits the common four-field initial state under the preset
1e-6 projection/jet policy. Both new short checks are NOT_RUN, with0 new ODE
and0 M-H/Timoshenko eigensolves. Additional floating strong/weak checks retain
an unresolved relative2e-12 qualification at p48/p64; no gate is weakened.
The old zero-axial case, its histories and nonlinear PARTIAL statuses remain.
No improvement of nonlinear convergence is established. V0/RHS/BC/basis and
coefficients are unchanged; LONG stays closed, EB/RLB-KV paused and the angular
same-clamp out-of-plane reference unavailable. No further strategy is selected.

## Exact-time leading axial response diagnostic (2026-10-08)

The [second-order axial note](theory/planar_second_order_axial_response.md)
extracts the leading forced u2,c2 response from the audited quartic action,
using the same continuous first Timoshenko eigenpair at every Shen resolution.
The amplitude convention is epsilon_a=A/h0, W=h0*w_hat and
Theta=h0*theta_hat; physical axial fields are epsilon_a² times this response.
All 2(p−1) discrete M-H coordinates are retained in the exact-time matrix
function. No new ODE integration, nonlinear feedback or production closure
is introduced, and the initial compatibility mismatch remains present.

Derivation, forcing/quadrature and exact-time controls pass. The bounded
p16/24/32/48/64 program and one conditional p96 clarification complete the
diagnostic, while **NLSP_SECOND_ORDER_SPATIAL_CONVERGENCE=PARTIAL** remains
separate from **NLSP_SECOND_ORDER_AXIAL_DIAGNOSTIC=COMPLETE**. Comparisons
with historical nonlinear trajectories use their actual complete intervals
or saved prefixes. Exact-time refers to the finite-dimensional forced system;
it does not establish a continuum or full nonlinear exact solution.
On the refined p64→96 pair only u2 meets both 1e−3 norm gates; c2 and both
velocities remain unresolved. The velocity difficulty is therefore already
present without time-integration error, and this does not assign its sole cause.
Both previous nonlinear PARTIAL statuses remain unchanged; full nonlinear
spatial convergence is unresolved. LONG remains closed, EB/RLB-KV paused,
and the angular same-clamp reference unavailable. No further stage is selected.

## Targeted planar solver diagnosis and recovery (2026-10-07)

The [diagnostic continuation](theory/weakly_nonlinear_planar_time_pilot.md)
reads the immutable first-pilot histories before any new integration. Physical
L2 projection separates unresolved spatial tails from different evolution in
the common space; the latter dominates the full-interval differences, without
establishing a phase-only explanation. The continuous initial form has an
axial acceleration trace mismatch of order A² at the clamps. This is a
qualification of boundary smoothness, not a code-error or invalid-IVP claim.

The existing helper now evaluates energy, gradient and Hessian only when
requested; action, variable mass, inertial terms, initial fields and BC are
unchanged. Pointwise equivalence and identical-settings old/new p32 controls
on 0...0.1T1 pass. The guarded p48 full-run forecast exceeds the remaining
fixed 900s budget, so **NLSP_PLANAR_P48_SPATIAL_CHECK=
REFINEMENT_DEFERRED_BY_BUDGET** and **NLSP_PLANAR_SOLVER_RECOVERY=PARTIAL**.
Exactly three short controls and no new full trajectory were performed.
The original two-amplitude pilot and its incomplete small-amplitude control
retain their historical PARTIAL status; no further refinement is selected.

## First four-field planar nonlinear time pilot (2026-10-07)

The [bounded numerical note](theory/weakly_nonlinear_planar_time_pilot.md)
extends the accepted cubic model to a fixed-fixed straight G20 initial-value
problem, with independent u,w,theta,c, two prescribed small amplitudes and
spatial/time convergence. It discretizes the quartic action, retains variable
mass and compares against continuous and semidiscrete linear references.
Actual trajectory/convergence statuses are recorded in that note; an algebra
PASS alone is not trajectory validation. No angular or out-of-plane stability
study is authorized; LONG remains closed and EB/RLB-KV remains paused.

**Bounded result: NLSP_PLANAR_TIME_PILOT=PARTIAL.** Both final p32 histories
reach5T1; linear and temporal controls pass, but u/c and several velocities
remain spatially unresolved. The final neighboring-p small-amplitude case
has a saved2.648T1 prefix after the fixed budget stop. Energy/mass remain
safe on computed histories; no complete nonlinear validation is claimed.

## Seven-field spatial nonlinear action audit (2026-10-07)

The [NLSP model note](theory/weakly_nonlinear_spatial_rod.md) records the
explicitly adopted nonlinear reduced V0 separately from frozen linear models.
One isolated helper/CLI compares variation of the quartic action with an
independent cubic expansion of the exact balances. All21 coefficient identities
and all21 supplied-draft parts agree; boundary/energy/reflection/planar/axial
checks pass. Manufactured jets verify fourth-order residual truncation;
bounded old-operator and straight-split controls use section-rotation clamps.
This is mathematical/model verification, not nonlinear physical validation.
No evolution, threshold, maps or spring/KV transfer. The longitudinal-model
selection remains CLOSED in its adopted 1D scope; EB/RLB-KV remains paused.
[Decision](memory/decisions.md#nlsp-d01), [evidence](memory/knowledge.md#nlsp-k01).

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
