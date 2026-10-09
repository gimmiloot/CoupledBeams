# FEM-3A: static preload, instantaneous release and short free motion

2026-10-09. Initial main HEAD `bcf5d96f5ede0b5deca994272f5d9787d129d203`,
clean working tree and index. [Authorized scope D13](../memory/decisions.md#nlsp-d13).
This is a bounded computational pilot for transfer/release/free-motion semantics,
not full nonlinear physical validation or a temporal/spatial convergence study.
Actual executed evidence and qualifications are recorded separately below.

**BLOCKED_BY_SOLVER:** the sole attempted linear preload+dynamic job terminated
with Windows access violation0xC0000005. No accepted static/dynamic output was
obtained; actual preload transfer, free motion and energy are unverified. The
nonlinear3D job and both1D trajectories were not started. The protocol/time
sections below describe checked inputs and installed-source semantics, not
successful dynamic execution.

## 1. Motivation and frozen scope

[FEM-2R](nlsp_nonlinear_static_3d_fem_validation.md#fem-2r-controlled-continuation-after-input-serialization-failure)
completed a qualified static comparison: both models reduce bending deflection,
with3.19353% common-scale max correction-profile difference. That result does not
prove all V0 coefficients or variable inertia. Its historical status remains
FEM2_DIAGNOSTIC_COMPLETE_WITH_QUALIFICATIONS. FEM-3A separately asks whether the
same static states can initialize short freely moving linear/nonlinear systems.

Same monolithic L=1,b=.20,h=.10; E=rho=1,nu=.3,kappa=5/6. Both full3D end faces
are fixed, lateral faces free, no internal joint. Positive1D w is negative global
Y. Static gravity remains g=.0014224751066856333,
q=rho*A0*g=2.844950213371267e-5. After release, gravity and every other external
force are zero, supports unchanged, physical initial velocities zero. No load,
initial amplitude/profile, constitutive coefficient or physical damping is fitted.

V0/quartic potential, cubic coordinate residuals, variable mass and all inertial
terms/Jacobian, four independent1D fields and Shen-Legendre space are unchanged.
1D essential values u=w=theta=c=0; no new slope constraints, inextensibility,
static condensation or prescribed c=-nu*u_s. The3D StVK law and full-face clamp
remain distinct from the reduced model; effective c is only a section diagnostic.

## 2. Immutable source data and initial states

| Source | Bundle | Pinned manifest SHA256 |
|---|---|---|
| FEM-2R | `results/nlsp_nonlinear_static_3d_fem_resume/210b74b8b166997c/` | `2ba6f23d41d49275c4e7e15b2c947258671c2a9b1b9a18a2bb8ef5d8bbfb1151` |
| Historical1D static FEM-2 | `results/nlsp_nonlinear_static_3d_fem/6714bd9f2778e6d7/` | `e0a972bf852c64658202bc5e385fabbede8ea5f5396248603e0ed6fe76d3d62c` |
| FEM-1 / medium mesh | `results/nlsp_linear_rectangular_3d_fem/4262efa427b03dad/` | `9f2d5139b84b2aa133b20d9a7cae806bac085c178fba506e087cae2941da1d90` |
| Frozen action | `results/weakly_nonlinear_spatial_rod/a9cedd4b6de99295/` | `04346246a1bf88a59d3d5a5055849313e3febe09f62f6929dd5bbcec7522f505` |

Source manifests/artifact hashes are checked before execution. No old mesh,
static equilibrium, FEM eigenfrequencies, nonlinear expansion or state is
regenerated to repair a cache hash. Saved p64 coordinates q_linear/q_nonlinear
are loaded directly with their old resting-mass whitening, no new L2 projection.
The full252 coordinates are four independent63-coefficient blocks, nq=129.
Each model/case starts at its own linear or nonlinear static equilibrium;
the slight L/NL initial-displacement difference is part of the defined pair.
Initial velocities are zero. The cases are not amplitude/phase aligned.

## 3. Installed-version protocol evidence

The installed CalculiX2.22 manual
`D:/PHD/CalculiX-Windows-master/src/downloads/ccx_2.22.pdf`, SHA256
`56963f827422ec7663cf218b60fffded19fd6ccebab793d2ccba667227d19d39`,
and local2.22 source/example archive members establish the protocol read-only.
Manual printed pages472/477/479 and599-601 describe bodyloads, direct dynamic
steps and STEP amplitude;535/537 describe initial velocities;490-496 and557-563
cover field/energy output. Source excerpts/hashes are retained as evidence,
without new literature search, installation or execution of historical examples.

A separate job for each L/NL case contains a static preload followed immediately
by a direct dynamic step, keeping nodal displacements, integration-point stresses
and internal state inside one solver execution. Saved FEM-2R U alone is not used
as a fictitious complete restart. Actual preload U/S/E/RF must reproduce the
corresponding saved medium state before the next production case is admitted.

The release explicitly clears previous bodyloads with DLOAD OP=NEW, then sets
the same GRAV magnitude to zero under STEP AMPLITUDE=STEP. Local dloads.f/bodyadd.f/
steps.f/tempload.f establish that the old magnitude is cleared and the zero
replacement is active from dynamic time0, without a full-period ramp-down.
No missing GRAV card is treated as proof that the old load disappeared.
Linear uses STEP NLGEOM=NO; nonlinear uses NLGEOM. Both use direct DYNAMIC,
not modal dynamics, with explicit ALPHA=0 and zero initial velocity cards.

Local initialconditionss.f/static solver/end-of-step routines establish zero
physical velocities and persistence of preload state. Native acceleration
initialization includes a tiny numerical time1.23571113e-20, so the internal
velocity after initialization may differ from exact zero by O(dt_tiny*a).
This is an output-accuracy qualification, not a prescribed physical velocity
or additional velocity boundary constraint.

## 4. Integration and energy policy

Fixed1D frequency from saved FEM-1 is omega1=.6054167303477958,
T1=2*pi/omega1=10.37828159055014. Authorized dynamic horizon is
.05*T1=.518914079527507; displayed tau=t/T1. This is the period scale of the
h=.10 linear1D beam, not the old h=.05 period or an inferred nonlinear period.

|3D policy quantity|Preset value|
|---|---:|
| Initial dynamic increment | T1/4000 = .00259457039763754 |
| Maximum increment | T1/2000 = .00518914079527507 |
| Minimum increment | 1e-4 of initial = 2.59457039763754e-7 |
| Maximum increments | 2000 |
| Direct integration alpha | 0 |
| Requested native output frequency | every accepted increment |

ALPHA=0 gives Newmark beta=.25,gamma=.5 and removes the requested HHT algorithmic
dissipation; it is not proof of zero time-discretization error. Installed native
initialization regularizes acceleration through
`[M+beta*(dt_initial/10)^2*(1+alpha)*K]*a0=Fext(0)-Fint(q0)`.
This built-in computational operation is retained; ALPHA=0 does not remove it.
No arbitrary new time refinement/damping or modification of native source occurs.

The planned, unexecuted1D linear reference uses exact-in-time evolution of the
complete finite-dimensional K,M0 system; any full252-coordinate spectral factorization is not modal truncation
or a new physical root search. Nonlinear1D uses one existing Radau integration
with analytic variable-mass Jacobian and the unchanged tight rtol/relative-atol
policy1e-10. Componentwise atol and cutoff-based max_step are evaluated for the
current h=.10 coefficients/characteristic scale, not copied as dimensional
numbers from h=.05. Execution remains EXPLORATORY_NOT_CERTIFIED/admitted=False;
old strict float64 strong/weak threshold2e-12 PARTIAL is retained.

After release,1D energy is .5*v.T*M(q)*v+V4(q), or the corresponding quadratic
linear energy.3D requests internal+kinetic energy, including energy output in
the preload step so history arrays exist. The potential of removed GRAV is not
included in free mechanical energy. Native energy-balance normalization differs
from normalization by the physical free initial energy; both are distinguished.
No energy projection or energy-based classification is used.

## 5. Transient metadata and physical recovery

The scoped reader distinguishes STATIC step1 from DYNAMIC step2, actual
increment, U/V/RF/S/E/energy datasets and total native time. Dynamic time is
`t_total-actual_static_end`. A static final output is stored as the physical
origin/continuity evidence, never relabelled a native dynamic t=0 frame; actual
positive-time dynamic samples retain their own metadata. Missing prefixes are
not filled, last states repeated or histories extrapolated.

Streaming fixed-width FRD parsing retains one nodal block at a time and checks
real node IDs/completeness/duplicates; modal readers are not used for transient
output. DAT U/RF are primary when available; V is requested in NODE FILE only
because the installed structural NODE PRINT does not support that variable.
No native structural acceleration output is assumed. Early U/V and the released
1D acceleration identity provide qualified restoring-sign evidence.

The same FEM-2R material-section recovery is used:41 weighted original-x slabs,
finite polar section rotation, selected81-section sensitivity, local coordinates
and positive-w sign. Effective c remains a thickness-strain projection diagnostic,
not a Cartesian FEM degree of freedom. Small recovered v/Phi/psi and section
residuals remain numerical diagnostics, without artificial planar constraints.

At actual3D timestamps, compare linear motion, total nonlinear motion and each
model's own NL-minus-L correction separately, retaining initial static offsets.
Use physical max/L2 differences and declared full-interval characteristic/fixed
scales, never instantaneous zeros; no amplitude/phase/time shift is fitted.
One medium mesh and one3D time policy cannot establish strict dynamic correction
resolution, full spatial/temporal convergence or PHYSICAL_DYNAMIC_VALIDATION_PASS.

## 6. Execution limits and stop

Only two medium C3D10 production jobs (5649nodes/3120elements), one nonlinear1D
p64 integration and the full linear semidiscrete reference are authorized.
A tiny load-release fixture is permitted only if the installed-version protocol
cannot otherwise be resolved; its use/counter must be explicit. Per3D job
<=1200s/4GiB, one thread, total numerical stage<=3600s, sequential execution.
A failed source/protocol/input/preload/output/solver gate preserves evidence and
stops the next production case; no hidden retry by changing code hash.

A new scoped [CLI](../../scripts/analysis/pilot_nlsp_nonlinear_dynamic_3d_fem.py)
and [config](../../data/input/nlsp_nonlinear_dynamic_3d_fem_pilot.json) route existing
physics/static/resource functions plus two narrow saved-state/transient helpers.
The distinct static-to-dynamic/release/time-metadata contract warrants a reusable
entry point under Script Proliferation Control; no general FEM framework is added.
The actual failed attempt and available preflight evidence are recorded below;
no missing trajectories, comparison or energy diagnostics are fabricated.

LONG CLOSED; EB/RLB-KV PAUSED_FOR_SUPERVISOR_DIRECTION; angular same-clamp reference
UNAVAILABLE; FEM-1R linear PASS, FEM-2R and planar physical sanity
DIAGNOSTIC_COMPLETE_WITH_QUALIFICATIONS, prepared strict PARTIAL are unchanged.
After the bounded report stop: no fine/refined or extra timestep dynamics, full
T1, new amplitude/geometry, angular joint, periodic orbit or Floquet study.


## 7. Available source/input and1D acceleration preflight

The separate bundle is `results/nlsp_nonlinear_dynamic_3d_fem_pilot/a69310e3bb30bab7/`.
Both intended L/NL input decks pass native finite numeric-field checks: maximum
width18, preset T1/duration/increments, correct linear/NL routing, unchanged
preload physics, explicit zero GRAV/STEP release, ALPHA=0 and zero velocity cards.
These checks establish the intended input contract, not successful execution.

The saved p64 initial-state check is successful without an equilibrium/root/ODE
solve. Raw and physical-profile reproduction errors and essential endpoint
residuals are zero in both cases. After release, `M*a0+grad(V)` relative residual
is1.12470e-16. The loaded equilibrium residual is1.76298e-14 for L and1.42785e-13
for NL, below the unchanged preflight gate. This semidiscrete released identity
does not repair or raise the historical independent strong/weak qualification.

| Saved1D initial quantity | L | NL |
|---|---:|---:|
| Midspan w | .004999999999999883 | .004991129920693128 |
| Midspan restoring acceleration | -.001418192734439063 | -.001418192734439028 |
| min(1+c) | 1 | .999981554828575 |
| Relative mass lower bound | 1 | .999963109997374 |
| Relative mass condition upper bound | 1 | 1.00003689136355 |
| Initial velocity coordinates | all zero | all zero |

The prepared1D tight time settings are rtol1e-10, componentwise504-entry atol
(min2.57172250e-16,max3.23507724e-13), max_step=.00720939298770064; all exact
entries/scales are stored in `one_d_preflight.json`. They were NOT used in an
ODE run. Preflight counters: RHS0,Jacobian0,linear eigendecompositions0. No full
linear exact-time factorization/trajectory was executed. Positive initial mass
and restoring acceleration are available evidence only at the restored states,
not along a nonexistent trajectory.

## 8. Actual first native attempt and hard-gate stop

One actual CalculiX2.22 command, in `cases/linear/`, was the unchanged binary
`results/_smoke/3d_fem_environment_check/calculix_2p22/CalculiX-2.22.0-win-x64/bin/ccx.exe motion`.
The exact input SHA256 is
`a38fb741b7638103bedee3893280266d32053cdae1b50d66804dd99369650998`.
It contains both intended preload and dynamic steps, but execution fails before
any accepted output can be verified. No input/physics change, new mesh, second
production case or hidden retry follows this failure.

| Attempt / evidence | Actual result |
|---|---|
| Production case | linear medium preload+dynamic, ordinal1 |
| Native return code | 3221225477 = 0xC0000005 access violation |
| Native wall time | 2.015955s |
| Peak working set | 25,305,088bytes = 24.13MiB |
| Recorded attempted-job stage | 2.047564s /3600s |
| Native stdout / stderr | 390bytes, version/build banner /0bytes |
| DAT / STA / CVG | all empty |
| FRD | 7bytes, no usable nodal block |
| Accepted static/dynamic increments or frames | none available |
| Actual static final time / load factor | not established |
| Actual dynamic start/end | not established; target.518914079527507 not verified |
| Nonlinear3D job | NOT_RUN |
| Nonlinear1D ODE / linear trajectory | NOT_RUN |
| Technical release fixture | NOT_RUN, zero calls |

This is a native process exception, not an observed nonlinear divergence,
resource-limit stop or parser failure on valid output. The resource wrapper did
not kill the process; the small memory/time figures are below the preset limits.
stdout buffering means the short banner cannot locate the crash at input reading
or before mechanical equilibrium. No displacement/reaction/stress/strain frame
or completed preload state is inferred from it.

Read-only Windows Application-event evidence identifies fault module ccx.exe,
offset0x2e1af6. PE/COFF lookup of31,242 symbols gives nearest preceding function
`resultsmech_` at RVA0x2dac10, fault offset+0x6ee6. This narrows the mechanical
results path; it does not identify an exact source line/invalid pointer by itself.
`crash_binary_audit.json` preserves this check without disassembly, rebuilding
or another solver call. Read-only source/binary tracing localizes the cause to the kinetic-energy
output path in the first LINEAR STATIC step: ccx_2.22.c1030-1037 frees veold;
resultsini.c313-331 enables ikin=1 for requested ELKE despite STATIC;
linstatic.c823-838 passes the dangling velocity pointer; resultsmech.f459-464
reads it. The fault RVA lies before the next resultstherm_ symbol and its MOVSD
read follows argument32 veold through the preserved stack/register trace.
This establishes a native output-path use-after-free, not an unknown physical
failure of the preload/release law. It still does not recover a completed
preload state or confirm that the dynamic step was reached.

A minimal current-generator preview requests only ELSE/ENER in the first linear
STATIC step while retaining ELKE in DYNAMIC, where velocities are allocated.
Geometry, gravity, constraints, release/ALPHA and physical initial states are
unchanged. That corrected preview is NOT_RUN; the actual failed deck and original
execution-code provenance remain authoritative. The user stop rule, "If any
mandatory gate fails, retain actual outputs and stop before the next production
job", prevents a hidden corrected retry within this task.

The input/source protocol still supports OP=NEW/zero GRAV/STEP semantics and
physical zero velocities, but actual transfer/release was not exercised to a
verifiable output. In particular no new preload U/S/E/RF or support result is
available to compare against FEM-2R, whose old static fields remain valid source
evidence. An absent native t=0 frame is not replaced by an invented transient.

## 9. Missing motion, correction and energy evidence

Neither linear nor nonlinear3D reached a verifiable0.05T1 endpoint. The workflow
stopped at the prescribed first production gate, before the next NL case and
before1D time-reference generation. There are no actual section displacement/
velocity histories or dynamic NL-minus-L corrections to report. Existing static
Delta w is not relabelled a measured dynamic correction.

| Requested comparison | Outcome |
|---|---|
| Reproduced native static preload | NOT_RUN / unverified |
| Actual complete load removal / zero initial native V | unverified; source/input definition only |
| Linear and nonlinear short free motion | not obtained |
|1D versus3D trajectories / nonlinear corrections | NOT_RUN |
| Native internal/kinetic energy or damping audit | NOT_RUN; no energy data |
| Temporal/spatial dynamic convergence | not established |

No dynamic plot, energy curve, damping diagnosis or fictional reached prefix is
created. Zero figures are appropriate for the failed execution; the planned
three-figure limit is not a requirement to manufacture missing observations.
Even a future successful one-mesh/one-timestep pilot would need separate accuracy
evidence; this failure adds no physical verdict about V0/inertia, nonlinear
periodic orbits, angular joints, out-of-plane stability or Floquet.

## 10. Reproducibility and qualifications

```powershell
python scripts/analysis/pilot_nlsp_nonlinear_dynamic_3d_fem.py --preflight
python scripts/analysis/pilot_nlsp_nonlinear_dynamic_3d_fem.py --run-pilot
python scripts/analysis/pilot_nlsp_nonlinear_dynamic_3d_fem.py --report-only results/nlsp_nonlinear_dynamic_3d_fem_pilot/a69310e3bb30bab7
python scripts/analysis/pilot_nlsp_nonlinear_dynamic_3d_fem.py --plot-only results/nlsp_nonlinear_dynamic_3d_fem_pilot/a69310e3bb30bab7
```

The commands replay the recorded failed authorization, not a new native attempt.
Cache identity includes source manifests/coordinates/mesh/action, installed
protocol, binary/runtime DLLs, config/model/time policies, code/environment and
explicit authorization; the failure ledger prevents a code-hash-triggered retry.
Full input, stdout/stderr, empty native files, resource/attempt records,
preflight and read-only crash evidence remain in the new bundle. Historical
FEM-1/FEM-1R/FEM-2/FEM-2R results/manifests/scripts are unchanged.

The distinction is explicit: (A) source/input semantics and1D state/acceleration
preflight are checked, but actual native preload transfer/feasibility remains
blocked; (B) observed short dynamic agreement cannot be assessed; (C) strict
3D time/space accuracy and nonlinear inertia validation are not established.
The concrete execution blocker, not a newly failed mathematical threshold,
prevents the pilot goal from being achieved. Further execution requires a separate
explicit technical decision. No next research stage or retry is selected here.


## 11. Final statuses and verification

| Status | Outcome |
|---|---|
| NLSP_FEM3A_SOURCE_PRESERVATION | PASS |
| NLSP_FEM3A_LOAD_RELEASE_PROTOCOL | PARTIAL |
| NLSP_FEM3A_STATIC_STATE_TRANSFER | NOT_RUN |
| NLSP_FEM3A_INITIAL_VELOCITIES | PARTIAL |
| NLSP_FEM3A_1D_DYNAMIC_REFERENCE | NOT_RUN |
| NLSP_FEM3A_3D_LINEAR_DYNAMIC | FAIL |
| NLSP_FEM3A_3D_NONLINEAR_DYNAMIC | NOT_RUN |
| NLSP_FEM3A_TRANSIENT_FIELD_RECOVERY | NOT_RUN |
| NLSP_FEM3A_ENERGY_DIAGNOSTICS | NOT_RUN |
| NLSP_FEM3A_SHORT_RESPONSE_COMPARISON | NOT_RUN |
| Overall | BLOCKED_BY_SOLVER |

Protocol/initial-velocity PARTIAL distinguishes supported installed-source and
input definitions (plus exactly zero1D initial coordinates) from absent actual
3D release/velocity evidence. Static transfer is NOT_RUN because no usable
preload output was obtained. Linear dynamic FAIL denotes the failed combined
production attempt; it does not prove that the dynamic step was entered.
The old strict float64 qualification remains PARTIAL separately.

The output-only current generator fix is retained under
`remediation_preview/linear_corrected_NOT_RUN.inp`, with its separate input gate
PASS/NOT_RUN. `post_attempt_changes.json` preserves original execution-code
SHA256 `68c46dc9728f12e9f176dcc338fc19828fdca2efa81a014885ace73717f61157`
and the later reporting/generator snapshot. Failed actual input/native logs,
one-attempt ledger and execution-code evidence are not replaced by the preview.
No further solver or fixture call occurred after failure.

**102 new focused tests PASS (12.89s)** and **42 selected historical regressions
PASS (2.87s)**, total144. Tests use synthetic fixtures and saved data with zero
real3D jobs/nonlinear integrations. They cover frozen-state loading, discrete
acceleration, complete linear reference using synthetic systems, transient
metadata/precision/missing/duplicate handling, releases/velocities/routing,
linear-static ELKE regression and failure-ledger replay. No old threshold/test
was relaxed; tests of a helper's synthetic arithmetic are not new production
trajectories. JUnit evidence is saved as `targeted_tests.xml` and
`historical_targeted_tests.xml` in this bundle. Final cache/link/hash/HEAD/index
and `git diff --check` evidence is retained separately after final review.


Final cache verification confirms all four matching routes (`--preflight`,
cached `--run-pilot`, `--report-only`, `--plot-only`) with numerical entry points
forbidden: zero CCX/Gmsh/1D ODE/static/root/BVP/eigen/symbolic calls. Failed input
and native log hashes are identical after replay; no histories, energies or
figures are fabricated. The failed-authorization guard survives the current
output-generator/reporting code-hash change and returns the preserved blocker,
not a corrected execution. Evidence is stored in the final cache-check artifact.
Old D/K exact byte prefixes, the entire FEM-2 canonical report and unique D13/K14
anchors were checked; `git diff --check` PASS. Root's final comprehensive source/
link/HEAD/index audit is retained separately with the completed bundle manifest.
