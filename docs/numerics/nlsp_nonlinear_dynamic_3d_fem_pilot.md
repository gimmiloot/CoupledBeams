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


<a id="fem-3ar-controlled-continuation-after-native-elke-output-path-failure"></a>

## FEM-3AR — controlled continuation after native ELKE output-path failure

2026-10-09. This is a new explicit authorization after the historical failure,
[NLSP-D14](../memory/decisions.md#nlsp-d14), with its own continuation namespace.
Initial clean main HEAD `51f9d69f13743c65b78ab763769fe8d8d31425a0`.
The FEM-3A record above remains BLOCKED_BY_SOLVER: its failed native job is not
reclassified, overwritten or replaced by the successful interpretation of a
corrected deck. The new request permits at most two sequential production
preload+dynamic jobs and, after both pass their actual gates, one nonlinear
1D trajectory and the complete exact-time semidiscrete linear reference.

### Continuation contract and immutable evidence

Parent `results/nlsp_nonlinear_dynamic_3d_fem_pilot/a69310e3bb30bab7/` has
manifest SHA256
`bf217d9c77d42170d1021b55c61873f8b2a94ed1fa67950d5ff9ed82c5cd257f`.
Its original failed input, logs, attempt ledger, native cause/source audit,
execution-code snapshots, corrected NOT_RUN preview and saved-state preflight
are checked under the parent's own hashes. FEM-2R, FEM-2, FEM-1 and the audited
action are verified under their existing manifests, not regenerated because
of a new CLI hash. The saved medium C3D10 mesh and p64 initial coordinates are
reused; no Gmsh, modal, static-only or new static-equilibrium calculation is
part of the continuation.

The separate [continuation CLI](../../scripts/analysis/resume_nlsp_nonlinear_dynamic_3d_fem.py)
and [config](../../data/input/nlsp_nonlinear_dynamic_3d_fem_resume.json)
change the authorization and immutable-parent/replay contract. This warrants a
thin orchestration entry point under Script Proliferation Control: merging a
new execution into the historically blocked authorization would make replay
unsafe. Generation, monitored native execution, parsers, material-section
recovery, saved-state loading, linear reference, nonlinear Radau and comparison
are reused from FEM-3A/FEM-2R. No second physics or general FEM solver is added.
An authorization-bound ledger remains authoritative after code changes; a
failure cannot be retried by changing a fingerprint.

### Output-only remediation and source gate

The previously established crash is a native velocity-buffer use-after-free
in LINEAR STATIC when ELKE forces kinetic-energy evaluation. The corrected
preload requests internal ELSE/ENER and excludes ELKE there; DYNAMIC retains
internal and kinetic outputs. The original bad input, saved corrected preview
and current generated linear deck are compared before production execution.
Only include-path context and the documented output correction may differ.
The linear check also recognizes explicit NLGEOM=NO as linear, rather than
mistaking the presence of the string NLGEOM for nonlinear routing.

Installed CalculiX 2.22 source semantics are separately checked for nonlinear
STATIC: that route retains allocated, initially zero velocities rather than
freeing the buffer as the linear route does. Its existing ELKE request therefore
does not reproduce the identified dangling-pointer path. This source finding
is a safety gate for the existing output prescription; it is not a new physical
model or a replacement for checking the real nonlinear native job.

Finite numeric fields remain bounded by native limits, including each separate
20-character field, with the same material, supports, gravity and static
controls. The release still uses a second same-job STEP with OP=NEW, zero GRAV,
AMPLITUDE=STEP and ALPHA=0. No source/binary rebuild, solver installation,
modified mass, artificial damping or adjusted integration parameter is used.

### Frozen state, time and comparison policy

All quantities in Sections 1–5 remain unchanged: L=1, b=.20, h=.10,
E=rho=1, nu=.3, kappa=5/6, g=.0014224751066856333 and
q=2.844950213371267e-5 before release. Physical w is minus global Y;
full solid end faces stay fixed. After release external gravity is zero.
Each linear/nonlinear calculation starts at its own equilibrium with zero
physical initial velocities, including its static stress state in the same
3D job. The 1D p64 coordinates are loaded directly without another projection,
Newton equilibrium or amplitude correction.

The horizon is .05*T1=.518914079527507, with saved omega1=.6054167303477958
and T1=10.37828159055014. Initial, minimum and maximum native increments,
2000-increment ceiling, output frequency and ALPHA are unchanged. The 1D
nonlinear reference keeps the existing tight componentwise tolerance policy,
variable mass, inertial terms and analytic Jacobian. Execution remains
EXPLORATORY_NOT_CERTIFIED, admitted=False; strict float64 2e-12 PARTIAL is
retained separately. No additional time or mesh convergence series is authorized.

Actual STATIC pseudo-time and positive DYNAMIC times are kept distinct.
The final preload may supply the physical initial observation but is never
relabelled a native dynamic t=0 frame. Comparison uses real matched native
timestamps, original material sections and the prescribed 41/81 recovery
sensitivity. Missing frames or unmatched times are not interpolated or filled.
L/NL corrections retain their initial static offset. Effective c is diagnostic
only; it is not identical to the reduced M-H coordinate. Phase, amplitude,
material coefficients and time scale are not fitted.


### Successful native execution and parser-only remediation

The corrected linear job is recorded in the separate continuation bundle
`results/nlsp_nonlinear_dynamic_3d_fem_resume/f6ac2f2f38510893/`.
It returns zero with real STATIC and DYNAMIC output through the prescribed
horizon, unlike the historical access violation. The first extraction encounters
an independent I/O limitation: the old reader assumes fixed absolute agreement
of printed total time minus printed step time. Native STA formatting such as
`0.151891E+01` has a decimal quantum of 1e-5; a valid pair of rounded columns can
therefore differ from an exact offset by more than that old assumption permits.

A minimal repair is confined to
[the existing transient output reader](../../scripts/lib/nlsp_fem3a_transient_output.py).
It retains raw printed tokens and derives their rounding intervals from the
actual lexical decimal precision, then checks interval compatibility and unique
increment association. It performs no interpolation, phase/time fitting,
scientific-gate change or modification of solver input. Incorrect or ambiguous
matches beyond the combined printed intervals are still rejected.
The same successful native files are reparsed without another CalculiX call;
both parser snapshots and the remediation evidence remain in the new bundle.
The source FEM-3A failed files and its parser/execution snapshots stay unchanged.

This correction concerns metadata arithmetic, not native equilibrium, continuous
PDE accuracy or a new tolerance for physical agreement. All prior preload,
strong/weak, safety and model-comparison definitions remain in force. New tests
use synthetic native-format examples and the saved successful linear STA;
no extra solver or integration is hidden in the reader regression.


### Native energy initialization qualification

The real linear job exposes a second, distinct native diagnostic limitation.
Its final STATIC internal energy is 3.736789e-8, consistent with independent
half displacement–body-load work 3.73678901575e-8. DYNAMIC stdout instead uses
initial energy 7.473578e-8, exactly twice the printed preload value. Normalizing
the unmodified native internal-plus-kinetic history by the physical preload
output gives a maximum relative change 1.00000013648. This approximately 100%
jump is retained explicitly; it is not hidden by replacing the denominator.

Read-only installed-source tracing localizes the bookkeeping path. On entry to
DYNAMIC, the routine allocates zero initial stress/strain baselines while copying
the previous internal-energy array. Initial-acceleration results then add the
preloaded elastic strain work again. The native reference is consequently
shifted. The mechanical force assembly uses the stress array, not this stored
energy value; actual preload U/S/E/RF reproduce the saved equilibrium exactly,
and recorded external and damping work after release are zero. This explains
the observed offset without changing the constitutive law, state, RHS or native
binary. It does not manufacture a corrected physical energy history.

A separately labelled drift relative to the solver's own DYNAMIC bookkeeping
reference is 6.82404e-8 for the linear job. This statistic is useful within that
native step; it is not substituted for the preload-referenced result or treated
as an independent physical conservation certificate. The complete native values,
source excerpts and independent work check are retained in
`energy_transfer_diagnostic.json`. Energy diagnostics remain PARTIAL even when
internal/kinetic outputs are present and the motion/release gates pass. ALPHA=0,
zero printed damping work and a small within-step bookkeeping drift do not
establish temporal accuracy. No solver energy correction or further test solve
is introduced in this continuation.


### Actual completed jobs, preload and release

Both corrected production jobs return zero and complete their own same-job static preload followed by .05T1 free motion. No second attempt or new fixture was used. Each preload reproduces FEM-2R medium U, support RF, S, E and recovered section fields exactly at their saved output precision. Zero differences are repetition evidence for the same finite mesh/output path, not a continuum certificate.

| Native case | Static increments | Dynamic increments / frames | First physical time | Final physical time | CCX seconds | Peak MiB |
|---|---:|---:|---:|---:|---:|---:|
| linear | 1 | 102 / 102 | 0.00259457 | 0.518914079527507 | 140.778 | 98.85 |
| nonlinear | 10 | 102 / 102 | 0.00259457 | 0.518914079527507 | 161.289 | 98.97 |

The native first dynamic frame is at positive time, not t=0. Both jobs have 102 exactly shared FRD time samples; final printed time is also retained separately from the input endpoint rounding. STATIC pseudo-time ends at 1.0 and is subtracted under explicit step metadata. There are no cutbacks or unexplained warnings. Fixed-face U/V stay zero; complete fields are finite and final deformation-gradient determinants are positive.

Independent support recovery subtracts consistent reference body-load contributions from RF. Linear force/moment imbalances are 4.84236e-8 / 3.12193e-9; nonlinear values are 1.51072e-8 / 7.98216e-9, below the unchanged 1e-5 gate. Applied force remains (0, -2.844950213371267e-5, 0). After release both native external work and damping work are zero. Initial w decreases and first w-velocity is negative, providing actual restoring free-motion evidence alongside the checked zero-velocity source/input contract.

**Nonlinear STATIC kinetic output qualification:** the native final STATIC ELKE is 9.631194e-8. It is not the physical kinetic energy of the released initial state. Installed checkconvergence.c computes predictor veold=(vold-vini)/dtime during STATIC; resultsmech can print kinetic energy from those pseudo-time velocities. nonlingeo.c explicitly zeros veold at the end of the static step before DYNAMIC. This is an output-semantics caveat, not an uninitialized-memory claim, a prescribed nonzero initial velocity or a transferred kinetic initial condition. Static ELKE is excluded from the free initial mechanical energy.

### Actual 1D references and short movement

The complete 252-coordinate linear K,M0 factorization gives the exact-in-time semidiscrete reference, with one factorization and no ODE/modal reduction or new physical root search. The single nonlinear p64/nq129 Radau integration reaches .518914079527507 from the frozen q0 and zero velocities. No new equilibrium or projection is used. The old callback authorization token is retained as execution provenance; separate metadata records FEM-3AR as the actual permission, without another integration.

| Movement | Initial midspan w | Final midspan w | Own initial-to-final change |
|---|---:|---:|---:|
| 1D linear | 0.005000000000 | 0.004808492610 | -1.9150739045e-04 |
| 1D nonlinear | 0.004991129921 | 0.004799622538 | -1.9150738291e-04 |
| 3D linear | 0.004821160850 | 0.004629628113 | -1.9153273657e-04 |
| 3D nonlinear | 0.004812087372 | 0.004620554660 | -1.9153271138e-04 |

Initial-to-final changes are observations of each defined IVP, not an amplitude/phase alignment replacing the full differences below. The inherited static 1D/3D offsets remain in every trajectory comparison. Final midspan NL-minus-L is -8.870071768e-6 in 1D and -9.073452872e-6 in 3D. Both remain negative; the difference includes each pair's initial static correction and is not isolated as a newly generated dynamic nonlinear effect.

Radau uses rtol=1e-10 and the same 504-entry componentwise atol (min 2.57172249937e-16, max 3.23507724041e-13), max_step=.00720939298770064. It accepts 811 internal steps, actual dt .000253098418–.000851629167, nfev5967, njev2, nlu248, mass factorizations5962; integration time6.009s. Total recorded numerical stage429.379s/3600s. One nonlinear ODE, two production CCX calls, zero fixtures/Gmsh/static-only/modal/physical-root/symbolic calls.

### Sampled physical differences and nonlinear correction

The tables use actual shared timestamps and the 41-point report grid in original material coordinates. Max is sampled max over x,t; L2 is max over sampled t of the physical length norm. Each relative max uses one fixed full-interval common characteristic scale, max of both compared profiles, not an instantaneous zero. These norms are model-to-model differences, not mesh/time error estimates or a continuum supremum.

| Pair / field | Absolute max | max_t L2 | Common characteristic scale | Relative max (%) |
|---|---:|---:|---:|---:|
| linear: u | 2.006578e-07 | 6.049128e-08 | 2.006578e-07 | 100.00000 |
| linear: w | 1.788702e-04 | 1.268681e-04 | 4.999995e-03 | 3.57741 |
| linear: theta | 5.601248e-04 | 3.764074e-04 | 1.365577e-02 | 4.10174 |
| linear: c_eff_diagnostic | 8.825888e-05 | 5.478406e-05 | 8.825888e-05 | 100.00000 |
| nonlinear: u | 1.078112e-06 | 5.254596e-07 | 5.530125e-06 | 19.49525 |
| nonlinear: w | 1.790736e-04 | 1.269854e-04 | 4.991125e-03 | 3.58784 |
| nonlinear: theta | 5.598079e-04 | 3.762804e-04 | 1.363078e-02 | 4.10694 |
| nonlinear: c_eff_diagnostic | 1.790806e-05 | 5.490617e-06 | 3.451091e-05 | 51.89102 |

| NL-minus-L comparison / field | Absolute max | max_t L2 | Common correction scale | Relative max (%) |
|---|---:|---:|---:|---:|
| u | 9.646108e-07 | 5.086260e-07 | 5.530125e-06 | 17.44284 |
| w | 2.034717e-07 | 1.210792e-07 | 9.073552e-06 | 2.24247 |
| theta | 3.364153e-07 | 2.161136e-07 | 2.562703e-05 | 1.31274 |
| c_eff_diagnostic | 8.773033e-05 | 5.614833e-05 | 1.058998e-04 | 82.84275 |

The w correction difference peaks at x=.5, t=.382699134 (tau=.03687), with max2.03472e-7 and L2 1.21079e-7: 2.24247% and 1.33442% of common correction scale9.07355245109e-6. Maximum sampled correction magnitudes are8.87008432162e-6 (1D) and9.07355245109e-6 (3D). Total w differences are about3.58%, dominated by the inherited static offset in this short interval. Correction agreement is less close for u (17.44%); c_eff (82.84%) remains a separate diagnostic and cannot certify the M-H coordinate. Linear u/c relative values near100% refer to very small or nonidentical quantities and must be read together with their absolute magnitudes.

Final-frame 41/81 recovery sensitivity for total w is max3.23915e-6 (linear) /3.23359e-6 (NL), L2 about1e-6; theta max1.74800e-5 /1.73431e-5. Small u and effective-c diagnostics are much more sensitive. These are observed recovery changes only at the final frame, not strict error bounds for the entire nonlinear correction history. DAT/FRD displacement agreement is5e-9; printed precision, recovery sensitivity and one retained 3D time policy limit interpretation of a signal around9e-6. No dynamic nonlinear signal certification or independent time/space uncertainty estimate is claimed.

### Energy and safety evidence

| Energy diagnostic | Linear | Nonlinear | Qualification |
|---|---:|---:|---|
| 1D max relative physical-energy drift | 4.10780e-13 | 3.61606e-12 | Own semidiscrete energy; removed gravity excluded |
| Native STATIC internal energy | 3.736789e-8 | 3.726886e-8 | Not augmented with static pseudo-kinetic output |
| Native DYNAMIC initial bookkeeping reference | 7.473578e-8 | 7.453772e-8 | Twice STATIC printed internal value |
| Max native change relative to STATIC | 1.00000013648 | 1.00000018782 | Raw approximately100% jump retained |
| Separate native DYNAMIC-reference drift | 6.82404e-8 | 9.39122e-8 | Bookkeeping diagnostic only |
| External / damping work after release | 0 /0 | 0 /0 | Native output, no fitted damping |

Energy status remains PARTIAL because of the established initialization offset and STATIC ELKE semantics, despite available histories. No energy is offset-corrected, projected or silently renormalized. Small 1D drift and small native within-step drift do not replace spatial/temporal accuracy checks.

Nonlinear 1D sampled min(1+c)=.999981554793, max|c|=1.84452e-5, mass Loewner lower bound .999963109927 and condition upper bound1.000036891434; all safety checks pass. Max|theta|=.0136590, max|u_s|=5.53134e-5, max|w_s|=.0150327 and L*max|theta_s|=.141992, below unchanged bounds. 3D fields/clamps/positive determinant checks pass on actual retained output; finite nodal/integration diagnostics are not exact maxima.

### Figures, reproducibility and bounded interpretation

Three actual-data figures are saved as PDF+PNG under the continuation bundle:

- `figures/linear_nonlinear_free_motion`: w(L/2), theta(L/4), u(L/4) for both models/cases.
- `figures/nonlinear_minus_linear_short_response`: midspan correction and final spatial profile.
- `figures/energy_and_release_diagnostics`: 1D physical drift, retained native +100% STATIC-reference jump and explicitly labelled native-DYNAMIC-reference drift.

```powershell
python scripts/analysis/resume_nlsp_nonlinear_dynamic_3d_fem.py --preflight
python scripts/analysis/resume_nlsp_nonlinear_dynamic_3d_fem.py --run-pilot
python scripts/analysis/resume_nlsp_nonlinear_dynamic_3d_fem.py --report-only results/nlsp_nonlinear_dynamic_3d_fem_resume/f6ac2f2f38510893
python scripts/analysis/resume_nlsp_nonlinear_dynamic_3d_fem.py --plot-only results/nlsp_nonlinear_dynamic_3d_fem_resume/f6ac2f2f38510893
```

The completed authorization is replayed from cache, without a second scientific execution. Bundle data retain immutable source identities, original and repaired parser snapshots, exact attempted inputs/native files, actual timestamps, two-attempt ledger, full1D coordinates/velocities, section histories, comparisons and energy qualifications. The corrected historical FEM-3A generator/failed guard and all source scientific bundles stay unchanged.

71 new focused tests (including the separate plotting-cache regression) and102 historical FEM-3A regressions pass, total173, using mocks/synthetic/saved evidence without extra scientific jobs. JUnit records are `targeted_tests.xml`, `plot_cache_regression.xml` and `historical_targeted_tests.xml`. The renderer validates the existing cache before figure mutation and saves its updated manifest afterward; it performs zero scientific calculations.

| FEM-3AR status | Outcome |
|---|---|
| NLSP_FEM3AR_SOURCE_PRESERVATION | PASS |
| NLSP_FEM3AR_INPUT_REMEDIATION | PASS |
| NLSP_FEM3AR_LINEAR_PRELOAD | PASS |
| NLSP_FEM3AR_LINEAR_RELEASE | PASS |
| NLSP_FEM3AR_LINEAR_DYNAMIC | PASS |
| NLSP_FEM3AR_NONLINEAR_PRELOAD | PASS |
| NLSP_FEM3AR_NONLINEAR_RELEASE | PASS |
| NLSP_FEM3AR_NONLINEAR_DYNAMIC | PASS |
| NLSP_FEM3AR_1D_REFERENCE | PASS |
| NLSP_FEM3AR_TRANSIENT_RECOVERY | PASS |
| NLSP_FEM3AR_ENERGY_DIAGNOSTICS | PARTIAL |
| NLSP_FEM3AR_SHORT_RESPONSE_COMPARISON | PASS |
| Overall | PILOT_COMPLETE_WITH_QUALIFICATIONS |

The pilot now demonstrates practical short 3D free motion and a meaningful unaligned comparison to the retained 1D model, including matching signs and similar bending corrections for this load/mesh/horizon. It does not establish PHYSICAL_DYNAMIC_VALIDATION_PASS, accuracy of every V0 coefficient/variable inertia, continuum truth of p64, temporal/spatial convergence of the 3D correction, physical stability, a nonlinear periodic orbit or critical amplitude. No unique causal allocation of discrepancies is established.

FEM-1R linear PASS; FEM-2R and planar physical sanity DIAGNOSTIC_COMPLETE_WITH_QUALIFICATIONS; strict float64/prepared and historical zero-u/c PARTIAL; LONG CLOSED; EB/RLB-KV PAUSED_FOR_SUPERVISOR_DIRECTION; angular same-clamp UNAVAILABLE remain unchanged. After this bounded report stop. No full T1, fine/refined dynamics, extra time level/amplitude, angular joint, out-of-plane perturbation or Floquet study is run or selected automatically.


### Additional read-only final sensitivity and cache checks

Paired final-frame recovery uses the same saved quadrature for linear and
nonlinear fields before subtracting them, rather than using the total-field
sensitivity as a bound on the small correction. Delta w from 41 versus 81
material sections differs by max7.78842e-9, L2 2.48752e-9, .0858312% of its
common correction scale. Final midspan values are -9.0734528721e-6 (41) and
-9.0741124007e-6 (81), difference6.59529e-10. The same-quadrature primary41
reproduction agrees with the stored w history within1.12e-16. Evidence is retained
in `paired_final_recovery_sensitivity.json/.npz`. This is a final-frame diagnostic,
not a whole-history recovery bound or temporal/spatial certification.

The paired nodal DAT decimal-rounding allowance is4.6989525e-9; it is explicitly
not propagated as a recovered-section or full solver error bound. These observed
measures support numerical visibility of the retained final correction. The
single selected medium mesh/time level still leaves independent dynamic
accuracy unestablished; `dynamic_nonlinear_signal_certified=false` is preserved.

Actual cached preflight, completed run-pilot, report-only and plot-only were
checked with scientific entry points forbidden: all four routes perform zero
new CCX/Gmsh/1D ODE/static/physical-root/BVP/symbolic calls. Native inputs/outputs
remain identical. `cache_readonly_verification.json` retains the audit. Plotting
updates figures and their new bundle manifest only, without a scientific solve.
The parent FEM-3A failed ledger/outputs and frozen physical helpers remain intact.
[Completed scoped memory K15](../memory/knowledge.md#nlsp-k15) records the result;
[decision D14](../memory/decisions.md#nlsp-d14) supplies its separate authorization.


**Dynamic-correction resolution qualification:** `midspan_summary.json` retains
the initial offsets separately. Maximum midspan Delta-w change from its own
preload is only7.92010e-12 in1D and1.66047e-10 in3D, below paired native nodal
DAT format allowance4.69895e-9. The observed2.24247% correction difference is
therefore largely inherited static difference over this short horizon, not
independent resolution of its nonlinear dynamic evolution. Main norms/profile
comparisons are unchanged, with no phase/amplitude alignment. Actual linear/NL
trajectories are available; accuracy of a newly evolving dynamic correction
remains unestablished.


<a id="fem-3b-nonlinear-correction-evolution-and-longer-horizon"></a>

## FEM-3B — эволюция нелинейной поправки на более длинном интервале

2026-10-09. Новое явное задание [NLSP-D15](../memory/decisions.md#nlsp-d15)
продолжает ту же задачу свободного движения после снятия GRAV. Прежний
FEM-3AR остаётся PILOT_COMPLETE_WITH_QUALIFICATIONS; его manifest, native
outputs, траектории и записи памяти не изменяются. Цель нового этапа —
измерить изменение поправки NL-minus-L относительно начального preload,
а затем отдельно сопоставить это изменение в 1D и 3D.

### Сохранённая задача и последовательность

Геометрия L=1, b=.20, h=.10 и материал E=rho=1, nu=.3, kappa=5/6
сохранены. Начальная нагрузка g=.0014224751066856333,
q=2.844950213371267e-5 не подбирается. Linear и nonlinear движения начинаются
из собственных сохранённых статических равновесий с нулевыми скоростями.
Положительное w соответствует отрицательной глобальной Y. После release
GRAV равна нулю; закреплённые торцы сохраняются.

Сначала read-only анализируется прежний интервал .05T1. Затем одна p64
nonlinear интеграция проходит до .5T1; при наличии корректных сохранённых
p48 координат выполняется один p48 control. Базис, все четыре независимых
поля, V0/V4, variable mass, inertial terms, RHS/Jacobian и tight time policy
не меняются. Полные exact-time семидискретные линейные решения используют
сохранённые q_linear и q_nonlinear без нового static solve или проекции.

До новых 3D-результатов выбирается один горизонт H=.25T1 либо H=.5T1,
где T1=10.37828159055014. Предварительно зафиксирована planning heuristic:
предпочесть .25T1, если midspan-сигнал эволюции в 1D превосходит десятикратно
каждый из исторических индикаторов 5e-9 и 7.78842e-9, а также наблюдаемую
p48/p64-разность. Эти индикаторы не являются строгими границами 3D-ошибки,
а правило не является критерием физической валидации. Иначе допустим .5T1
с явной qualification и проверкой safety/resources.

После фиксации решения разрешены два sequential CalculiX jobs на прежней
medium mesh: 5649 nodes, 3120 C3D10 elements. Каждый job содержит STATIC
preload и DYNAMIC в одной execution. Исправленный output routing, OP=NEW,
zero GRAV, STEP amplitude, ALPHA=0 и нулевые физические начальные скорости
сохраняются. Initial increment=T1/4000, maximum=T1/2000, minimum=1e-4 initial,
maximum increments=2000. После первого solver failure вся новая серия
останавливается без retry. Новые meshes, modal и static-only jobs запрещены.

### Определения и пределы сравнения

Для каждой модели сохраняются полная поправка Delta w=w_NL-w_L и её
эволюционная часть delta_evol w=Delta w(t)-Delta w(0). В 3D значение при
нулевом физическом времени берётся из подтверждённого static preload.
Оно не переименовывается в native dynamic t=0 frame. Знаки, пространственные
профили и собственные начальные состояния сохраняются без phase/amplitude
alignment.

Дополнительная 1D-диагностика раскладывает полную поправку на линейную
эволюцию разности начальных состояний и nonlinear-minus-linear движение
из одного nonlinear initial state. Тождество проверяется в сохранённых
координатах и физических полях. Такое разложение не переносится автоматически
на 3D: отдельное линейное 3D движение из nonlinear preload здесь не задано.

Используются actual native timestamps. При необходимости явно маркированная
postprocessing interpolation допускается только внутри фактического overlap;
сохраняется разность linear/PCHIP. Interpolated values не считаются native
samples. Все L2/max показатели относятся к сохранённой сетке по x,t,
а не к доказанному непрерывному supremum. Effective c_eff остаётся отдельной
3D-диагностикой и не тождественно координате M-H c.

Native STATIC/DYNAMIC energy discrepancy сохраняется без вычитания константы.
Отдельно показываются raw jump, drift внутри DYNAMIC, доступная независимая
kinetic-energy диагностика и собственная mechanical energy в 1D. Одна medium
mesh и один native time-step level не дают независимой 3D temporal/spatial
сертификации. Strict float64 PARTIAL и EXPLORATORY_NOT_CERTIFIED сохраняются;
PHYSICAL_DYNAMIC_VALIDATION_PASS этим этапом не назначается.

### Организация и воспроизводимость

Новый [config](../../data/input/nlsp_nonlinear_dynamic_long_horizon.json) имеет
отдельную authorization `explicit_user_FEM3B_2026_10_09`, immutable parent
manifest и namespace `results/nlsp_nonlinear_dynamic_long_horizon/`.
Существующий [continuation CLI](../../scripts/analysis/resume_nlsp_nonlinear_dynamic_3d_fem.py)
получает явный `--long-horizon` режим. Его прежний FEM-3AR default, authorization
ledger и cache semantics сохраняются. Scoped
[helper](../../scripts/lib/nlsp_fem3b_continuation.py) переиспользует прежние
generator, monitored runner, streaming readers, section recovery и 1D solver;
новый runnable FEM solver не создаётся. Matching cache/report/plot не выполняют
новых scientific calls.

### Stage A — что показали сохранённые .05T1 данные

Собственный manifest FEM-3AR имеет SHA256
`187870ece572dfe1b999b83d99016d01c1fdbe10e8835ed45543d1119d6f36ea`;
проверены все 959 записанных artifact hashes. Read-only анализ не выполнял
CCX, Gmsh, Radau, static equilibrium, eigenanalysis или BVP.

| Сохранённый показатель | 1D | 3D |
|---|---:|---:|
| Initial midspan Delta w | -8.8700793068e-6 | -9.0734780570e-6 |
| Final midspan Delta w | -8.8700717676e-6 | -9.0734528721e-6 |
| Max midspan delta_evol w | 7.9200994721e-12 | 1.6604665646e-10 |
| Max по x,t delta_evol w | 1.2448739465e-8 | 1.7346753932e-8 |
| max_t physical L2 evolution | 6.0718314705e-9 | 8.1957546499e-9 |
| Место spatial maximum x/L | .200 | .175 |
| Spatial evolution / initial correction max | .0014034530 | .0019118087 |

Малость midspan-поправки не является следствием ошибочного вычитания или
повторённого snapshot. Использованы реальные timestamps, независимый
static initial offset и все сохранённые формы. Однако spatial maximum
эволюции уже достигает порядка 1e-8 вне середины стержня. Поэтому вывод
«эволюция ниже output indicators» относится к midspan-наблюдению;
его нельзя автоматически переносить на весь профиль.

Сохранённые начальные 1D ускорения w совпадают: разность в середине
3.51282e-17, maximum по профилю 4.40837e-16. Оба сохранённых equilibrium
имеют одну и ту же transverse variational force, а принятая kinetic energy
содержит одинаковый constant translational mass block для w. После снятия
нагрузки это объясняет отсутствие существенного t²-вклада в Delta w.
Вывод относится к сохранённой semidiscrete задаче и не определяет все
последующие члены временной эволюции.

Приближение w(0)+w_tt(0)t²/2 воспроизводит основное раннее изменение
midspan. К концу .05T1 его разность с 1D составляет около 5.68e-7,
то есть .296% начального-to-final изменения. Для 3D ускорение оценено
из первой positive-time velocity, поэтому это finite-time estimate,
а не фактически выведенный native acceleration при t=0.

Из сохранённых коэффициентов линейной теории Тимошенко получены характерные
скорости sqrt(S/m)=.5661385171 и sqrt(Bp/jp)=1. Соответствующие half-length
времена .8831760866 и .5 дают ориентир ранней стадии движения. Coupling,
dispersion и глобальный Galerkin basis не позволяют считать эти оценки
доказанным временем прихода волнового фронта или единственной причиной
малости поправки.

Независимая kinetic energy вычислена из реальных FRD velocities существующей
положительной 14-point degree-5 C3D10 квадратурой по reference volume.
Проверены affine midside geometry, положительные determinants и total mass=.02.
Maximum разности с native kinetic energy составляет 2.13295e-14 (linear)
и 2.13225e-14 (nonlinear), около 5.31e-6 собственного фиксированного K-scale.
Это отдельная mass-weighted диагностика; она не исправляет внутреннюю энергию
и не делает native energy accounting физически сертифицированным.

Raw STATIC-to-DYNAMIC reference jump остаётся +100% в обоих старых jobs.
Первый positive-time independent K равен около 1.36591e-13; external/damping
work равны нулю. Drift внутри native DYNAMIC reference равен 6.82404e-8
и 9.39122e-8. Независимое восстановление StVK internal energy имеет NOT_RUN:
проверенного общего constitutive/element evaluator здесь нет и новый FE
implementation ради этой диагностики не создаётся.

### Stage B — 1D сигнал и зафиксированный выбор горизонта

Обе разрешённые nonlinear интеграции из сохранённых p64/p48 координат
достигли .5T1=5.189140795275070. Initial coefficients не перепроецировались;
начальные скорости нулевые. Сохранены полные q/velocity и accepted Radau
cubic polynomials. Их повторная оценка на выбранных times воспроизводит
сохранённые состояния с maximum 6.62e-24 (p64) и 2.65e-23 (p48), без нового
ODE solve. Все прежние динамические DOF остаются независимыми.

| Предварительный 1D prefix | S(H), midspan max evolution | Spatial max evolution | max_t L2 evolution | p48/p64 max evolution difference |
|---|---:|---:|---:|---:|
| 0.25T1 | 3.5856198942e-06 | 3.5856198942e-06 | 2.0833400275e-06 | 1.0244278486e-10 |
| 0.50T1 | 1.8004816179e-05 | 1.8004816179e-05 | 1.1045204350e-05 | 2.1270618039e-10 |

Выбран .25T1=2.594570397637535. Его ожидаемый 1D сигнал 3.58562e-6
превышает заранее установленную planning requirement 7.78842e-8.
p48/p64 sensitivity определена консервативно как maximum по всем 41 сечениям
на том же prefix, а не только в середине. Safety пройдена; свободного места
и numerical budget достаточно. Решение сохранено до первого нового native
attempt в `horizon_decision.json` вместе с frozen config. Оно не выбиралось
по будущему совпадению 1D/3D.

Стоимость двух CCX jobs оценена по прежним counters: около 704s linear
и 806s nonlinear для .25T1; примерно 510 accepted increments на case.
Это planning estimate, не гарантия линейного масштабирования времени.
Индивидуальный timeout сохранён 1200s, память 4GiB, один поток.

Точность именно evolving w не следует смешивать с общим four-field PASS.
На полном предварительном .5T1 прежние all-eight spatial gates имеют PARTIAL:
проходят u,w displacement; остальные six компоненты превышают хотя бы один
прежний L2/max criterion. Пороги, numerical floor и знаменатели не изменены.

| Компонента, p48/p64 до .5T1 | Absolute max | max_t L2 | Relative max | Relative L2 | Preset tolerance | Status |
|---|---:|---:|---:|---:|---:|---|
| u | 1.4072070e-09 | 5.0990206e-10 | 2.5400601e-04 | 1.5322242e-04 | 1e-03 | PASS |
| w | 1.5171190e-07 | 4.4015085e-08 | 3.0035103e-05 | 1.3640162e-05 | 1e-04 | PASS |
| theta | 1.7939244e-06 | 6.0934193e-07 | 1.2845682e-04 | 6.1094283e-05 | 1e-04 | FAIL |
| c | 2.5638735e-08 | 8.6010249e-09 | 1.3031392e-03 | 4.8658428e-04 | 1e-03 | FAIL |
| u_t | 8.8528615e-08 | 3.1893158e-08 | 1.1678783e-02 | 7.3896397e-03 | 1e-03 | FAIL |
| w_t | 1.1061976e-05 | 2.7548773e-06 | 3.6515637e-03 | 1.4041049e-03 | 1e-04 | FAIL |
| theta_t | 1.4576711e-04 | 3.9189251e-05 | 1.5761507e-02 | 6.1941158e-03 | 1e-04 | FAIL |
| c_t | 1.8998311e-06 | 5.9135270e-07 | 8.9261161e-02 | 4.4114501e-02 | 1e-03 | FAIL |

Малую p48/p64-разность evolving w около 1.0e-10 на .25T1 нельзя использовать,
чтобы снять эту qualification или объявить p64 точным continuum solution.
Отдельного 1D time refinement в FEM-3B не выполнялось: остаётся прежняя tight
policy с rtol=1e-10 и componentwise atol, analytic Jacobian и max_step
.0072093929877006385.

| Nonlinear 1D case, .5T1 | nDOF / nq | Solver seconds | Accepted steps | RHS / Jacobian / LU | Own energy drift |
|---|---:|---:|---:|---:|---:|
| p64 | 252 / 129 | 55.6446 | 8024 | 58828 / 2 / 2296 | 3.62142e-11 |
| p48 | 188 / 97 | 27.8147 | 6276 | 45704 / 2 / 1422 | 2.40217e-11 |

На .5T1 p64 min(1+c)=.9999803258; relative mass lower bound=.9999606520.
Поля, slopes и curvature остаются в прежних safety bounds. Полный exact-time
линейный reference использует все Shen coordinates; это не modal truncation.
Дополнительное linear движение из nonlinear initial state требует только
того же сохранённого полного eigensystem, без ещё одной nonlinear интеграции.

### Stage C — два фактически завершённых native jobs

Оба новых CalculiX 2.22 jobs завершились с returncode=0 и реальным
JOB FINISHED. Каждый содержит собственный STATIC preload и DYNAMIC в одной
execution. Preload U/S/E/RF и section profiles воспроизводят сохранённые
FEM-2R medium данные точно в пределах native output precision. Applied force
остаётся (0,-2.844950213371267e-5,0). Независимые force/moment imbalances
равны 4.84236e-8/3.12193e-9 для linear и 1.51072e-8/7.98216e-9 для nonlinear;
прежний equilibrium gate 1e-5 сохранён.

| Case | STATIC increments | DYNAMIC increments / frames | CCX seconds | Peak MiB | Reached physical horizon | Cutbacks |
|---|---:|---:|---:|---:|---:|---:|
| linear | 1 | 502 / 502 | 691.863903 | 99.4492 | 2.594570397637535 | 0 |
| nonlinear | 10 | 502 / 502 | 729.338809 | 99.4492 | 2.594570397637535 | 0 |

First native dynamic time=.002594570000000074; оба набора имеют 502 точно
совпадающих timestamps. Interpolation не потребовалась. Нулевые GRAV,
external work и damping work подтверждены; ALPHA остаётся нулём. Нет
необъяснённых warnings. Это actual displacement trajectories, а не вывод
только по exit code или заранее предполагаемому solver capability.

Сохранённые поля конечны, fixed-end U/V равны нулю. Проверка deformation
Jacobian выполнена в существующих 14 volume quadrature points каждого C3D10
на всех 502 сохранённых dynamic frames, также на сохранённых STATIC frames.
Observed min detF=.9936290399 (linear) и .9936553037 (nonlinear); maximum
Green-strain magnitude=.00704544/.00707669. Это sampled extrema на указанных
точках и временах, не доказанные непрерывные extrema в элементе или времени.

На общем участке FEM-3AR совпадают 101 native dynamic frames до
.517616794: сохранённые U/V/S/E/RF/ENER воспроизводятся без разности.
Собственный последний старый кадр .518914079527507 не является точно общим
после изменения t_bound. Он не интерполировался и не заменялся новым кадром.
1D prefix оценивается в точных старых times сохранённым dense output и
линейным operator; source q0 остаётся идентичным. Разница few-ulp initial
back-transform отдельно сохраняется как arithmetic qualification.

### Stage D — движение и изменение поправки

Оба 1D сравнения вычислены в фактические 3D times из сохранённых accepted
Radau polynomials и полного linear eigensystem. Новых ODE/eigen solves
для postprocessing нет. Нормы ниже относятся к 41 material sections и
502 native timestamps. У каждого показателя один фиксированный characteristic
scale на всём горизонте. Начальные amplitudes, phases и time scale не
выравниваются; 3D не объявляется точным continuum reference.

| Midspan movement | Initial physical w | Final w at .25T1 |
|---|---:|---:|
| 1D linear | 5.000000000000e-03 | -6.440990783203e-05 |
| 1D nonlinear | 4.991129920693e-03 | -6.969436724462e-05 |
| 3D linear | 4.821160849899e-03 | -1.971551818700e-04 |
| 3D nonlinear | 4.812087371842e-03 | -2.023605683438e-04 |

Отрицательный final w означает, что наблюдение уже перешло через нуль в
данной свободной IVP. Это не solver failure и не условие периодичности.
T1 остаётся периодом первой линейной 1D формы; nonlinear period не искался.

| Model difference / field | Absolute max | max_t L2 | Common characteristic scale | Relative max (%) |
|---|---:|---:|---:|---:|
| linear: u | 2.0065778e-07 | 6.0491281e-08 | 2.0065778e-07 | 100.00000 |
| linear: w | 1.8770425e-04 | 1.2686814e-04 | 4.9999952e-03 | 3.75409 |
| linear: theta | 5.6012480e-04 | 3.7640736e-04 | 1.3655773e-02 | 4.10174 |
| linear: c_eff_diagnostic | 8.8258882e-05 | 5.4784065e-05 | 8.8258882e-05 | 100.00000 |
| nonlinear: u | 1.0781118e-06 | 5.2545964e-07 | 5.5301247e-06 | 19.49525 |
| nonlinear: w | 1.8788768e-04 | 1.2698542e-04 | 4.9911251e-03 | 3.76444 |
| nonlinear: theta | 5.5980792e-04 | 3.7628037e-04 | 1.3630776e-02 | 4.10694 |
| nonlinear: c_eff_diagnostic | 1.8866759e-05 | 5.4906172e-06 | 3.4510909e-05 | 54.66897 |
| Delta: u | 9.6461081e-07 | 5.0862598e-07 | 5.5301247e-06 | 17.44284 |
| Delta: w | 2.0448465e-07 | 1.2107918e-07 | 9.0735839e-06 | 2.25363 |
| Delta: theta | 4.6185220e-07 | 2.6843553e-07 | 2.5627029e-05 | 1.80221 |
| Delta: c_eff_diagnostic | 8.7730328e-05 | 5.6148325e-05 | 1.0589982e-04 | 82.84275 |
| evolving Delta: u | 9.5705118e-07 | 5.0134172e-07 | 5.6777568e-06 | 16.85615 |
| evolving Delta: w | 2.8247169e-07 | 1.8211832e-07 | 3.8680916e-06 | 7.30261 |
| evolving Delta: theta | 5.2836824e-07 | 2.8481797e-07 | 1.1716720e-05 | 4.50952 |
| evolving Delta: c_eff_diagnostic | 8.7690580e-05 | 5.6104037e-05 | 1.0605722e-04 | 82.68233 |

Величины u/w дополнительно нормированы h=.1, а dimensionless theta/c —
фиксированным unit scale; эти данные сохранены в JSON/CSV. Effective c_eff
остаётся диагностикой поперечных 3D-деформаций. Его большое различие не
представляется как совпадение отдельной M-H DOF или незаметное малое поле.

Полная Delta-w разность имеет maximum 2.04485e-7, или 2.25363% своего общего
масштаба. Для evolving части разность увеличивается до 2.82472e-7:
7.30261% по max и 4.70822% по L2 при общем фиксированном масштабе 3.86809e-6.
Maximum достигается в середине при t=.25T1, signed 1D-minus-3D=-2.82472e-7.
Таким образом, близость полной поправки около 2.25% не заменяет самостоятельное
сравнение её развивающейся части.

| Evolving w on 0… .25T1 | 1D | 3D |
|---|---:|---:|
| Max abs | 3.5856198942e-06 | 3.8680915832e-06 |
| max_t L2 | 2.0833400275e-06 | 2.2622255288e-06 |
| Final signed midspan | 3.5856198942e-06 | 3.8680915832e-06 |
| Evolution / max initial correction | 4.0423763646e-01 | 4.2630748197e-01 |

Initial midspan Delta w=-8.870079307e-6/-9.073478057e-6 в 1D/3D.
Final значения -5.284459413e-6/-5.205386474e-6 сохраняют отрицательный знак,
но существенно меняются относительно preload. Эволюция составляет около
40.42% и 42.63% initial correction scale. На старом .05T1 midspan-изменение
было ниже наблюдаемого разрешения; новый горизонт даёт самостоятельный
динамический сигнал.

### Дополнительное 1D разложение

| Signed midspan component at .25T1 | Value | Over fixed full-horizon total-correction scale |
|---|---:|---:|
| Total NL-L | -5.2844594126e-06 | -0.5957620219 |
| Linear evolution of initial-state difference | -2.1126310406e-07 | -0.0238174852 |
| Nonlinear-minus-linear evolution from the same NL initial state | -5.0731963085e-06 | -0.5719445366 |

На выбранном конце остаточная linear initial-state component мала сравнительно
с same-IC nonlinear component. Относительно конкретного ненулевого final
Delta w их additive signed contributions составляют примерно 4.00% и 96.00%.
Это описание одной фиксированной конечной точки, не instantaneous percentage
на всей траектории и не energy/modal fractions. Сумма компонентов проверена
до arithmetic precision. Такое разложение отдельно не рассчитывалось для 3D;
поэтому указанные доли нельзя автоматически приписать 3D движению.

### Разрешённость и оставшиеся численные ограничения

Observed 3D evolving-w signal=3.868091583e-6. DAT/FRD difference=5.0e-9;
paired initial/final 41/81 evolving-correction difference=9.386996962e-9.
Signal превосходит наибольшее из этих измеренных различий примерно в 412 раз.
Native pairing не потребовала interpolation, её diagnostic difference равна
нулю. В этом ограниченном смысле нелинейная динамическая эволюция различима
на сохранённой сетке наблюдений. Эти различия не являются строгими верхними
границами continuum error.

Раздельные qualifications сохраняются: 1D evolving w мало чувствительна к p,
но общий full-.5T1 all-eight spatial check PARTIAL; 1D использует один tight
level без нового temporal control; 3D имеет одну medium mesh и одну исходную
time-step policy. Independent 3D temporal/spatial certification отсутствует.
Для количественного physical-accuracy утверждения нужны отдельно разрешённые
time/mesh controls именно Delta/evolving Delta, а также решение оставшихся
four-field 1D spatial limitations. Эти действия в FEM-3B не запускались.

### Энергия на новом горизонте

| Raw/native energy diagnostic | Linear | Nonlinear |
|---|---:|---:|
| STATIC internal energy | 3.7367890000e-08 | 3.7268860000e-08 |
| DYNAMIC bookkeeping initial reference | 7.4735780000e-08 | 7.4537720000e-08 |
| STATIC-to-DYNAMIC reference jump | 1.0000000000e+00 | 1.0000000000e+00 |
| Drift within native DYNAMIC reference | 1.3380471864e-07 | 1.3416026134e-07 |

Старый +100% reference discrepancy сохраняется без вычитания константы.
External/damping work в DYNAMIC равны нулю. Независимая kinetic energy
из фактических nodal velocities даёт maximum difference 6.91922e-14 linear
и 6.45969e-14 nonlinear, около 1.86e-6/1.74e-6 fixed kinetic scale.
Она использует проверенную существующую C3D10 mass/quadrature и не заменяет
недостающую independent internal-energy реконструкцию; последняя NOT_RUN.
Native ENERGY_DIAGNOSTICS остаётся PARTIAL. Нулевой ALPHA и малый within-step
drift не доказывают точность 3D интегрирования.

### Сохранённый результат и остановка

Bundle `results/nlsp_nonlinear_dynamic_long_horizon/7d2b499e6a1eb990/`
содержит source/provenance hashes, old-data diagnostic, обе p64/p48 nonlinear
траектории, полные linear references и same-NL-initial-state auxiliary reference,
accepted dense output, frozen pre-FEM decision/config, два actual native jobs,
DAT/FRD/STA/logs, section histories, all-eight spatial table, correction/evolution
NPZ/CSV, recovery/kinetic diagnostics и три PDF/PNG figures. Numerical stage
2367.679s из 6000s; число production CCX calls=2, nonlinear 1D calls=2.
Gmsh/modal/static-only/physical-root/symbolic calls=0.

| FEM-3B status | Result |
|---|---|
| `NLSP_FEM3B_SOURCE_PRESERVATION` | PASS |
| `NLSP_FEM3B_OLD_SIGNAL_DIAGNOSTIC` | PASS |
| `NLSP_FEM3B_1D_PRELIMINARY` | PASS |
| `NLSP_FEM3B_HORIZON_SELECTION` | PASS |
| `NLSP_FEM3B_LINEAR_3D` | PASS |
| `NLSP_FEM3B_NONLINEAR_3D` | PASS |
| `NLSP_FEM3B_PREFIX_REPRODUCTION` | PASS |
| `NLSP_FEM3B_RESPONSE_COMPARISON` | PASS |
| `NLSP_FEM3B_DYNAMIC_NONLINEAR_SIGNAL` | PASS |
| `NLSP_FEM3B_ENERGY_DIAGNOSTICS` | PARTIAL |

Overall FEM3B_DIAGNOSTIC_COMPLETE_WITH_QUALIFICATIONS. Новый динамический
сигнал измерен и сопоставлен; PHYSICAL_DYNAMIC_VALIDATION_PASS не присваивается.
Физика V0/V4, variable mass, RHS/Jacobian, coefficients, материал, геометрия,
нагрузка, заделки и source static states сохранены. Старые bundles и ledgers
не перезаписаны. Новых meshes/modal/static-only jobs, fullT1, других amplitudes,
fine/refined dynamics, automatic timestep refinement, angular joints, periodic
orbits, Floquet и critical-amplitude search нет. После этого отчёта остановка.

Воспроизведение существующего завершённого cache:

```powershell
python scripts/analysis/resume_nlsp_nonlinear_dynamic_3d_fem.py --long-horizon --compute
python scripts/analysis/resume_nlsp_nonlinear_dynamic_3d_fem.py --long-horizon --report-only results/nlsp_nonlinear_dynamic_long_horizon/7d2b499e6a1eb990
python scripts/analysis/resume_nlsp_nonlinear_dynamic_3d_fem.py --long-horizon --plot-only results/nlsp_nonlinear_dynamic_long_horizon/7d2b499e6a1eb990
```

Три figures: `linear_nonlinear_free_motion`, `nonlinear_correction_and_evolution`,
`one_d_initial_state_evolution_decomposition`. Underlying data находятся в
`dynamic_comparison.npz`, `comparison_metrics.csv`, `midspan_response.csv`,
source trajectories и native frame files.
[Результат памяти NLSP-K16](../memory/knowledge.md#nlsp-k16) сохраняет выводы
и ограничения; исторические D14/K15 остаются неизменными.


### Почему малость старой полной поправки не означает отсутствия nonlinear effect

Дополнительный 1D same-initial-state reference показывает существенное взаимное
сокращение двух additive компонентов уже на .05T1. Полная разность исходных
IVP мала относительно Delta w(0), хотя nonlinear-minus-linear эволюция из
одного nonlinear initial state уже не мала в том же смысле.

| 1D midspan at .05T1 | Signed value |
|---|---:|
| Initial correction Delta w(0) | -8.870079306756e-6 |
| Linear initial-state component I(t) | -8.344349028554e-6 |
| Change I(t)-I(0) | +5.257302782020e-7 |
| Same-IC nonlinear component N(t) | -5.257227391623e-7 |
| Total evolving correction | +7.539039705018e-12 |

Таким образом, исчезновение заметного изменения полной Delta w в середине
на коротком интервале подтверждено как cancellation в принятой 1D декомпозиции.
Это не утверждение, что intrinsic nonlinear evolution ещё отсутствовала.
Равенство начальных ускорений объясняет начальное отсутствие t²-вклада
в полную исходную разность; оно не зануляет каждый компонент по отдельности.
Wave-time estimates остаются дополнительной ограниченной интерпретацией.
Для 3D такая same-initial-state декомпозиция не выполнялась, поэтому её
причинное разделение здесь не объявляется установленным.


### Targeted verification

98 focused new/selected-historical tests PASS. Они используют synthetic fixtures
и сохранённые outputs; дополнительных FEM/ODE jobs в tests нет. Проверены
algebraic decomposition, frozen-source reuse, horizon/authorization/ledger,
dense accepted-output equivalence, explicit native/interpolated pairing,
signal-resolution arithmetic, fresh-compute finalization/call accounting и
старые cache/failure/output-safety contracts. При завершённом cache scientific
entrypoints запрещены: compute/report/plot не повторяют solves.

Postcompletion исправлены только routing/finalization и read-only completion
недостающих linear/decomposition artifacts из сохранённых полных factors.
При отсутствии factors eigenanalysis retry запрещён. Actual execution-code
snapshots сохранены отдельно от final reporting-code revisions. Default
FEM-3AR, его старый failed guard, физический helper и scientific thresholds
не изменены этим оформлением результата.


<a id="fem-3c-numerical-robustness-and-dissertation-verification"></a>

## FEM-3C — устойчивость динамического сопоставления и диссертационный итог

### Отдельное разрешение, frozen sources и заранее заданные gates

Задание FEM-3C продолжает [FEM-3B](#fem-3b-nonlinear-correction-evolution-and-longer-horizon)
с отдельной authorization `explicit_user_FEM3C_2026_10_09`, описанной в
[NLSP-D16](../memory/decisions.md#nlsp-d16). Новый
[config](../../data/input/nlsp_nonlinear_dynamic_validation.json) и namespace
`results/nlsp_nonlinear_dynamic_validation/c6256269eb8143ef/` сохраняют
immutable parent FEM-3B и его manifest SHA256
`e4bcc291fab5f04a2fb73103c35a2fa6481aecf5f058b8ce6ae4a5661efa0c4e`.
Старые bundles, failure ledgers, научные статусы и записи D/K не переписываются.

Физическая задача сохраняет L=1, b=.20, h=.10; E=rho=1, nu=.3, kappa=5/6,
g=.0014224751066856333 и q=2.844950213371267e-5. Сохраняются собственные
linear/NL static equilibria, нулевые скорости, fixed faces, instantaneous
OP=NEW/zero-GRAV/STEP release, ALPHA=0 и same-job STATIC/DYNAMIC transfer.
V0/V4, variable mass, RHS/Jacobian, Shen–Legendre basis и все независимые
координаты остаются неизменными. Исторический strict float64 PARTIAL,
EXPLORATORY_NOT_CERTIFIED и admitted=False сохраняются.

FEM-3C1 заранее ограничен четырьмя sequential production jobs на .25T1:
linear/NL на medium с уменьшенным вдвое шагом, затем linear/NL на сохранённой
fine mesh при той же временной policy. Medium содержит 5649 nodes/3120 C3D10,
fine — 11553 nodes/6670 C3D10; новый Gmsh mesh не создаётся. Каждый preload
сравнивается со своим mesh-level FEM-2R reference, а не с medium reference
для обеих сеток.

| Stage | Mesh | Horizon | Initial / maximum increment | Output cadence | Administrative increment limit |
|---|---|---|---|---:|---:|
| C1 temporal | medium | .25T1 | T1/8000 / T1/4000 | every 2 accepted increments | 2000 |
| C1 spatial | fine | .25T1 | T1/8000 / T1/4000 | every 2 accepted increments | 2000 |
| C2, conditional illustration | medium | T1 | T1/4000 / T1/2000 | every 5 accepted increments | 3000 |

T1=10.37828159055014; контрольный конец H=2.594570397637535. Minimum increment
равен 1e-4 initial, ALPHA=0. STATIC output остаётся frequency1; необходимые
DYNAMIC U/V/RF/S/E/ENER и internal/kinetic energy outputs явно следуют выбранной
cadence. Фактический final accepted frame обязателен даже при прореженном output.
Это изменение вывода, а не физики, интегратора или шага. Материал, нагрузка,
заделки и физические input cards сохранены; INCLUDE сравнивается после проверки
разрешённого абсолютного пути и SHA256 одного frozen mesh.

До первого нового результата зафиксированы 201 common physical times на
[0,.25T1], включая концы, primary linear и diagnostic PCHIP interpolation.
Точка t=0 берётся из подтверждённого static preload и не выдаётся за native
DYNAMIC frame. Все Delta w=NL-minus-L и delta_evol w=Delta w(t)-Delta w(0)
рассчитываются после согласованного переноса обеих траекторий. Phase/amplitude
fitting и extrapolation отсутствуют.

Baseline E_model=2.824717e-7 сохраняет историческое абсолютное расхождение.
Заранее объявлено Rt=Dt/E_model<=.25, Rh=Dh/E_model<=.25, где Dt и Dh — sampled
max differences эволюционной поправки при temporal и spatial refinement.
Interpolation difference дополнительно должна быть <=.25 min(Dt,Dh) и
<=.25 E_model. Это практические диагностические ориентиры, не теорема о FEM,
не строгие error bounds и не физический порог универсальной применимости.

Full-period 3D разрешён только после всех четырёх actual completed/gated C1 jobs,
численной robustness и нового resource preflight. Он использует только medium
mesh и заранее заданную original time policy. Результаты на T1 являются
иллюстративными: сходимость на .25T1 не переносится автоматически на весь период.
T1 — период первой ЛИНЕЙНОЙ формы принятой 1D модели, не найденная nonlinear
periodic orbit. Предусмотрены один p64 и максимум один p48 nonlinear 1D solve
до T1 из saved static coordinates плюс full exact-time linear references из
сохранённых факторов; eigenanalysis и новый static solve не выполняются.

Существующий [continuation CLI](../../scripts/analysis/resume_nlsp_nonlinear_dynamic_3d_fem.py)
получил явный `--validation` preset с scoped
[orchestration](../../scripts/lib/nlsp_fem3c_validation.py),
[read-only diagnostics](../../scripts/lib/nlsp_fem3c_diagnostics.py) и
[saved-state 1D helper](../../scripts/lib/nlsp_fem3c_1d.py).
Нового runnable FEM solver нет; исторические default/--long-horizon cache и
failure guards сохраняются. Frozen code и actual prelaunch/postprocessing
snapshots различимы. Обнаруженные до первого native call расхождения relative
INCLUDE из директорий разной глубины и output-cadence inheritance устранены
проверяемым input/output adapter; это не solver retries. Legacy generator,
математические helpers и научные thresholds не изменяются.

До запуска оценён bounded numerical cost 18522.628s с postprocessing allowance,
по фактическим FEM-3B counters. Это planning estimate, не гарантия. Per-job timeout
5400s, total numerical budget 24000s, память 4GiB и один поток фиксированы.
Максимум 6 новых CCX jobs и 2 nonlinear 1D integrations. Первый реальный solver
failure останавливает последовательность; code-hash change не даёт retry.
Новые meshes/modal/static-only jobs и дальнейшие автоматические исследования
не разрешены.

### FEM-3C1A — фактически завершённый temporal control

Оба medium jobs прошли completion, preload, release, finite-field, fixed-U/V и
sampled positive-detF gates. Preload U/S/E/RF/section profiles повторяют собственные
FEM-2R medium states с нулевой измеренной разностью. Каждый job имеет 1002 accepted
DYNAMIC increments и 501 actual saved frames, достигая H. Warnings и cutbacks
отсутствуют; external/damping work после release равна нулю.

| Medium refined-time job | CCX seconds | Peak MiB | Initial midspan w | Final midspan w | Min sampled detF |
|---|---:|---:|---:|---:|---:|
| linear | 1121.826110 | 99.5703 | .004821160849899094 | -.000197142378497060 | .993629043 |
| nonlinear | 1185.491228 | 99.4063 | .004812087371842098 | -.000202347611824860 | .993655306 |

Read-only preliminary 201-grid comparison сохранён в `temporal_preliminary.json`.
Он не разрешает full-period 3D до получения spatial control.

| Temporal preliminary quantity | Existing medium | Medium refined-time |
|---|---:|---:|
| Max evolving w correction | 3.868091583184e-6 | 3.868244729196e-6 |
| Evolving correction at H, midspan | +3.868091583184e-6 | +3.868244729196e-6 |
| 1D/3D max evolving-w difference | 2.824716890162e-7 | 2.826248350274e-7 |
| 1D/3D max-time spatial L2 difference | 1.821183200207e-7 | 1.821173843239e-7 |
| Relative max on each own common evolution scale, % | 7.30261120 | 7.30628114 |
| Relative max on fixed historical scale3.868091583184e-6, % | 7.30261120 | 7.30657041 |

Dt=5.400904445896e-10, Rt=.001912016123, то есть изменение равно .191202% baseline
model discrepancy и проходит заранее заданный .25 ratio. Максимальная linear/PCHIP
разность evolving w на двух medium resolutions — 9.17254e-11; совместный comparison
с min(Dt,Dh) требует ещё spatial результата. Проценты с обновлённым и историческим
знаменателями сохранены отдельно; исторические 7.3026% не переписаны.

Dt сопоставимо с half-quantum 7-significant-digit nodal DAT output: при U около .005
оно составляет 5e-10. Поэтому Dt характеризует наблюдаемое совокупное изменение
при уменьшении шага вместе с округлением вывода и postprocessing interpolation.
Оно не является отдельно разрешённой чистой temporal error или её строгой верхней
границей. Сам evolving-w signal около 3.868e-6 существенно больше наблюдаемых
output/recovery differences. Эта qualification не меняет заранее выбранные
Rt/Rh thresholds и не отменяет необходимость fine comparison.

Raw native STATIC/DYNAMIC energy reference jump +100% сохранён. Within-DYNAMIC
reference drift medium L/NL — 1.338047e-7/1.341603e-7; native data не перенормируются.
Independent kinetic energy рассчитывается прежней C3D10 mass/quadrature policy;
independent internal StVK reconstruction NOT_RUN. ENERGY_DIAGNOSTICS остаётся PARTIAL.

На момент данного промежуточного раздела fine spatial pair ещё не завершена,
full-period 1D/3D results NOT_RUN. Общий verification result пока не присваивается.
Последующий фактический результат дополняет этот раздел, сохраняя temporal evidence.


### FEM-3C1 — completed temporal and spatial robustness, 2026-10-10

Все четыре разрешённых C1 production jobs завершены. Новая fine pair имеет
11553 nodes/6670 C3D10, 1002 accepted DYNAMIC increments и 501 actual output frames
на каждый case при той же refined-time policy, что medium. Fine preload повторяет
свой собственный FEM-2R U/S/E/RF/section reference с нулевой разностью.
Все completion/release/finite/clamp/positive sampled-detF/equilibrium gates PASS;
warnings и cutbacks отсутствуют. Исторические FEM-3B jobs не повторялись.

| Actual C1 job | CCX s | Recovery s | Peak MiB | Dynamic increments / saved frames | Final midspan w |
|---|---:|---:|---:|---:|---:|
| medium refined-time linear | 1121.826110 | 312.008715 | 99.5703 | 1002 /501 | -.000197142378497060 |
| medium refined-time nonlinear | 1185.491228 | 314.805116 | 99.4063 | 1002 /501 | -.000202347611824860 |
| fine refined-time linear | 2921.618411 | 634.882226 | 210.2109 | 1002 /501 | -.000189931647633494 |
| fine refined-time nonlinear | 3070.336089 | 642.462548 | 210.2031 | 1002 /501 | -.000195192612941096 |

Совместный read-only 201-grid результат сохранён в `robustness_comparison.json`,
`robustness_comparison.npz` и `robustness_summary.csv`. Дополнительных scientific
calls при обработке нет. Пространственное сравнение использует одинаковые
refined time settings и исходные материальные сечения. Каждый собственный
STATIC t=0 profile явно отделён от positive-time native samples.

| Quantity | Existing medium | Medium refined-time | Fine refined-time | Observed change |
|---|---:|---:|---:|---|
| Max evolving NL w correction | 3.868091583184e-6 | 3.868244729196e-6 | 3.887471189460e-6 | max values change +1.5314601e-10, then +1.9226460e-8 |
| Midspan evolving correction at .25T1 | +3.868091583184e-6 | +3.868244729196e-6 | +3.887471189460e-6 | signed increases; same location/endtime |
| 1D/3D max evolving-w difference | 2.824716890162e-7 | 2.826248350274e-7 | 3.018512952918e-7 | +1.5314601e-10, then +1.9226460e-8 |
| 1D/3D max-time spatial L2 difference | 1.821183200207e-7 | 1.821173843239e-7 | 1.941159534063e-7 | -9.3569682e-13, then +1.1998569e-8 |
| Temporal/spatial change divided by fixed baseline 2.824717e-7 | baseline | Rt=.001912016123 | Rh=.074609339668 | both <=.25 predeclared guide |

Изменения максимальных значений в первых строках не подменяют максимальную
разность пространственно-временных функций. Главные observed refinement
metrics — Dt=5.400904445896e-10 и Dh=2.107502701193e-8. Temporal maximum signed
change=-5.400904445896e-10, spatial=+2.107502701193e-8; положения/времена и
cumulative curves сохранены. Rt соответствует .191202% baseline model difference,
Rh — 7.460934%. Это observed changes, не строгие temporal/spatial error bounds.

| Evolving-w model normalization | Existing medium | Medium refined-time | Fine refined-time |
|---|---:|---:|---:|
| Own pair characteristic max scale | 3.868091583184e-6 | 3.868244729196e-6 | 3.887471189460e-6 |
| Relative max on own pair scale, % | 7.30261120 | 7.30628114 | 7.76472109 |
| Relative L2 on own pair scale, % | 4.70822151 | 4.70801092 | 4.99337343 |
| Relative max on fixed historical scale 3.868091583184e-6, % | 7.30261120 | 7.30657041 | 7.80362328 |
| Relative max on fixed all-resolution scale 3.887471189460e-6, % | 7.26620662 | 7.27014610 | 7.76472109 |

Обновлённая модельная разность после spatial refinement слегка увеличивается;
исходные 7.3026% не сохраняются искусственно. Таблица отдельно показывает
старый знаменатель, pair-specific denominators и общий denominator для
всех resolutions. Absolute/L2 differences относятся к sampled 41x201 domain,
без удаления boundary/initial regions и без phase fitting.

Linear/PCHIP evolution-w differences old/mediumRT/fineRT равны 9.1725445e-11,
6.3546902e-11 и 8.5550607e-11. Их максимум составляет .16983349 min(Dt,Dh) и
.00032472437 baseline; interpolation-comparability gate PASS. Dt близок
к исходному DAT nodal rounding indicator 5e-10 при U~.005. Поэтому результат
поддерживает малое совокупное изменение при time refinement и сохранение
эволюционного сигнала, но не выделяет отдельно чистую temporal error из
округления и postprocessing. Новый порог по масштабу округления после просмотра
данных не вводится; ratio/interpolation policies сохранены.

`NLSP_FEM3C_TEMPORAL_CONTROL`, `NLSP_FEM3C_SPATIAL_CONTROL` и
`NLSP_FEM3C_ROBUSTNESS` = PASS в объявленном диагностическом смысле.
Результат достаточно устойчив для ограниченного quantitative bending comparison
на .25T1. Это не EXACT_3D_CONTINUUM_REFERENCE и не physical validation всех V0
coefficients/семиполевых nonlinear couplings. Научная формулировка приведена
в [dissertation summary](nlsp_straight_rod_3d_fem_verification_summary.md).
Full-period program остаётся отдельным conditional этапом; actual results
и окончательный статус добавляются после завершения. Quarter-period robustness
не переносится на T1 автоматически.


### Full-period 1D trajectories and pre-FEM3C2 decision, 2026-10-10

После completed C1 получены обе разрешённые nonlinear 1D trajectories до
T1=10.37828159055014. p64 является основным comparator, p48 — spatial diagnostic.
Saved static q_nonlinear каждого p переиспользованы точно, v0=0; нового static
solve, initial projection, BVP, eigenanalysis или modal truncation нет. Complete
linear factors из FEM-3B дают exact-time linear references из собственных q_linear
и из q_nonlinear для additive decomposition. Все 252/188 координаты независимы,
новых derivative BC нет; сохранены четыре поля и четыре скорости.

| Nonlinear 1D case | p64 | p48 |
|---|---:|---:|
| Spatial quadrature / coordinate DOF | 129 /252 | 97 /188 |
| Actual end | 10.37828159055014 | 10.37828159055014 |
| Output samples | 28499 | 28499 |
| Accepted Radau steps | 16006 | 12524 |
| RHS / Jacobian / LU calls | 117588 /2 /4748 | 91330 /2 /2926 |
| Mass factorizations | 117583 | 91324 |
| Integration seconds | 115.168661 | 56.960507 |
| rtol | 1e-10 | 1e-10 |
| Componentwise atol entries | 504 | 376 |
| Min / max componentwise atol | 2.571722499368e-16 /3.235077240413e-13 | 2.977456670877e-16 /3.745467216092e-13 |
| max_step | .0072093929877006385 | .0072093929877006385 |
| Own initial-energy relative drift | 7.272492585099e-11 | 4.853427737969e-11 |

Output sampling preserves 12 samples per highest retained linear frequency,
старые common times и required 0,.05,.25,.5,.75,1 times/T1. Оно не заменяет
временной convergence test. Actual accepted Radau dense polynomials сохраняются
для read-only evaluation; воспроизведение sampled states отличается максимум
на 5.29396e-23/6.61744e-24. Новые nonlinear integrations — ровно 2, linear
eigendecompositions — 0; linear dynamics сохраняет все прежние DOF.

p64 safety samples дают min(1+c)=.999978734453, relative mass lower bound
.999957469358, upper bound 1.000004149342 и condition bound1.000042532451.
max|c|=2.12655472e-5, max|theta|=.01434791965, max|w_s|=.01576015180,
max L|theta_s|=.1419920812; retained axial/shear strains 7.11122e-5/.00221867645.
Все прежние safety gates пройдены. E=0.5*v.T*M(q)*v+V4(q) использует собственную
initial energy; GRAV potential после release отсутствует. Никакой energy
projection, field filtering или RHS/Jacobian correction не применялись.

Full-period all8 spatial comparison использует прежние physical norms,
characteristic scales и numerical-floor policy без изменения thresholds.

| Component | Absolute max-time L2 difference | Absolute max-space-time difference | Relative L2 | Relative max | Preset tolerance | Status |
|---|---:|---:|---:|---:|---:|---|
| u | 5.22770648e-10 | 1.44021181e-9 | .000136329671 | .000222657100 | 1e-3 | PASS |
| w | 4.40150853e-8 | 1.60100073e-7 | .000013640162 | .000031234777 | 1e-4 | PASS |
| theta | 6.25307263e-7 | 1.79392442e-6 | .000061295523 | .000125073325 | 1e-4 | FAIL |
| c | 1.37172108e-8 | 4.11639954e-8 | .000743918344 | .001935861925 | 1e-3 | FAIL |
| u_t | 3.41589028e-8 | 9.38491103e-8 | .004279005212 | .007579325821 | 1e-3 | FAIL |
| w_t | 2.75487733e-6 | 1.39780804e-5 | .001404104891 | .004297964295 | 1e-4 | FAIL |
| theta_t | 4.17178354e-5 | 1.45767110e-4 | .006346214069 | .015011550142 | 1e-4 | FAIL |
| c_t | 1.00383735e-6 | 3.25824158e-6 | .062713041475 | .113646555013 | 1e-3 | FAIL |

Только u,w displacement проходят оба прежних criteria; overall all8 spatial
qualification PARTIAL. Непройденные criteria не скрываются малым energy drift,
близостью w или диаграммой семи полей. p64 не является exact continuum truth;
полная numerical certification four-field dynamics не заявляется. Отдельного
нового 1D temporal refinement в FEM-3C нет, сохраняется прежняя tight policy.

После actual C1 checks заморожен `full_period_decision.json`: numerical robustness,
preload/release/recovery и resource preflight PASS, поэтому условные 2 medium
CCX full-period jobs разрешены. На момент решения remaining numerical budget
13577.022172s; estimate для native pair 5684.810850s плюс 1800s postprocessing
вмещается. Original initial/max=T1/4000,T1/2000, output every 5, INC 3000 выбраны
в исходном pre-result config, без подбора по будущему 1D/3D совпадению. Размер шага
не укрупняется ; 3000 — только administrative accepted-increment allowance.

На момент этой записи full-period 3D cases ещё не завершены. Иллюстративный
one-period comparison, его displacement/correction/energy tables и окончательный
научный вывод добавляются только после фактических результатов. Сходимость C1
на .25T1 не объявляется доказанной на полном T1; no full-period fine/time study.


### Read-only 1D correction sensitivity on the full period

`one_d_correction_spatial_diagnostic.json/.npz` отдельно сравнивает сохранённые
p48/p64 Delta w, delta_evol w и additive decomposition. Максимальная разность
полной поправки на T1 равна 4.368619881129e-10, spatial L2 — 1.710369822038e-10;
для её эволюции 4.368619945097e-10 и 1.710369823462e-10. Последняя разность равна
.0015465690705 baseline C1 model discrepancy и .0014472755338 updated fine C1
model discrepancy. На quarter-prefix evolving-w difference 1.024427848594e-10.
Таким образом, main bending-correction comparator существенно менее чувствителен
к выбранным p, чем некоторые полные поля и скорости. Correlated linear/NL
cancellation может уменьшать эту разность, поэтому она не подменяет all8 PARTIAL
и не даёт strict continuum error bound.

На полном saved-grid characteristic scales Delta w=1.785704042240e-5 и
delta_evol w=2.672711972923e-5; это maxima over 41 material sections/28499 times,
не continuum supremum. Scales для quarter-window и полного T1 сохраняются отдельно.
1D correction growth и изменение знаменателя не смешиваются с обновлённой 3D
C1 модельной разностью.

Удлинённый 1D run сохраняет 14301 exact old FEM-3B times на [0,.5T1], bitwise q0/v0
и исходные physical initial fields. Prefix differences составляют примерно 3e-15
для nonlinear w и не более 5.77e-12 для relevant velocities; last samples могут
отличаться вследствие удлинённого t_bound и последнего accepted dense polynomial.
Линейная reference использует те же saved factors, без нового eigenanalysis.
Это reproduction diagnostic, а не ещё одна ODE integration или независимый
time refinement. Старые FEM-3B файлы остаются неизменными.


### FEM-3C2 — actual full-period illustration completed, 2026-10-10

Оба conditional medium preload+dynamic jobs достигают T1. Каждый содержит
2002 accepted DYNAMIC increments и 401 actual positive-time frames при заранее
выбранном output frequency 5. Final accepted state выведен; missing prefix не
заполняется последним значением. Original initial/max=T1/4000,T1/2000, ALPHA=0,
zero-GRAV/STEP release и все material/load/clamp settings сохранены. STATIC
preload U/S/E/RF/section profiles повторяют собственные medium FEM-2R equilibria
с нулевой разностью. Completion, finite fields, fixed-end U/V, independent
reaction balance и positive sampled-detF PASS; warnings/cutbacks отсутствуют.

| Actual full-period job | CCX s | Recovery s | Peak MiB | Dynamic increments / frames | Actual end |
|---|---:|---:|---:|---:|---:|
| linear | 1961.864073 | 254.045876 | 99.5508 | 2002 /401 | 10.37828159055014 |
| nonlinear | 2154.299542 | 256.693777 | 99.4805 | 2002 /401 | 10.37828159055014 |

Сохраняются 100 exact shared native timestamps с FEM-3B quarter-prefix каждого
case; STATIC datum совпадает отдельно. W difference linear=0, nonlinear около
1.2e-21; остальные recovered fields отличаются максимум на arithmetic scale
около 2.1e-14. Полное bitwise coincidence recovery не заявляется; unmatched
samples не интерполируются ради этой reproduction check. Native output every5
объясняет иное число frames; никаких дополнительных CCX solves для prefix нет.

| Model / case | Initial midspan w | Final midspan w at T1 | Initial Delta w | Final Delta w | Final delta_evol w |
|---|---:|---:|---:|---:|---:|
| 1D linear p64 | .005000000000000 | .005134569565480 | — | — | — |
| 1D nonlinear p64 | .004991129920693 | .005125309947300 | -8.870079306756e-6 | -9.259618180290e-6 | -3.895388735325e-7 |
| 3D medium linear | .004821160849899 | .004889861345824 | — | — | — |
| 3D medium nonlinear | .004812087371842 | .004877646304410 | -9.073478056996e-6 | -1.221504141406e-5 | -3.141563357066e-6 |

T1 задаётся первой линейной 1D частотой. Ни nonlinear periodic orbit, ни q(T1)=q0
не предполагались; linear loaded static profile также содержит другие retained
modes. Поэтому отличие конечного состояния от начального само по себе не является
ошибкой интегратора. На comparison grid midspan w maximum 1D достигается при
10.32639018, 3D nonlinear при 10.11882455. Наблюдается различие времени extrema,
согласующееся с сохранявшейся разностью линейных частот; его единственная причина
не выделена, phase/amplitude/time fitting и новый period search не выполнялись.

### Full-horizon quantitative differences and seven-field recovery

`full_period_comparison.json/.npz` использует 401 common physical times и 41
исходное материальное сечение. Primary linear/diagnostic PCHIP transfer сохраняет
actual native data отдельно; t=0 — confirmed STATIC state. Все absolute max/L2
и normalization scales относятся к полному интервалу, sampled maxima only.

| Full-period w comparison | Absolute max | Max-time spatial L2 | Declared common characteristic scale | Relative max /L2, % |
|---|---:|---:|---:|---:|
| Linear 1D minus 3D | 4.843367702484e-4 | 2.840388216388e-4 | 5.136938014182e-3 | 9.428511 /5.529341 |
| Nonlinear 1D minus 3D | 4.845878707917e-4 | 2.841297344404e-4 | 5.128394451779e-3 | 9.449115 /5.540325 |
| Total Delta w difference | 2.955423233772e-6 | 1.775486564906e-6 | 1.824704270381e-5 | 16.196724 /9.730270 |
| Evolving delta_evol w difference | 2.752024483533e-6 | 1.654889705716e-6 | 2.732052076081e-5 | 10.073104 /6.057314 |

Evolving max difference достигается при x=.5,t=T1. На прежнем fixed historical
scale 3.868091583184e-6 это 71.14683%, а абсолютная разность составляет 9.74266
старого E_model baseline. Увеличение full-horizon denominator не объявляется
улучшением: абсолютная разность выросла по сравнению с четвертью периода.
Полная Delta w и её эволюция используют разные denominators и не смешиваются.
Full-period time/mesh convergence не сертифицирована результатом C1; это
иллюстративный longer-response comparison.

Four-active full nonlinear comparisons не сводятся к хорошо разрешённому w:

| Active field | Absolute 1D/3D max difference | Own characteristic scale | Relative max, % | Meaning |
|---|---:|---:|---:|---|
| u | 2.062367484812e-6 | 6.465992182219e-6 | 31.89561 | actual section-coordinate comparison |
| w | 4.845878707917e-4 | 5.128394451779e-3 | 9.44911 | bending displacement |
| theta | 1.373552156838e-3 | 1.431652417134e-2 | 9.59417 | finite section orientation |
| c vs c_eff | 2.342770126204e-5 | 3.959002504924e-5 | 59.17577 | proxy comparison only |

Следовательно, statement об agreement всех четырёх nonlinear fields не делается.
3D c_eff включает finite director-stretch/quadratic measurement даже для
линейного Cartesian displacement field; linear1D M-H c почти нулевое при этой
нагрузке. Это не ошибка отдельной FEM c-DOF, поскольку такой DOF нет, и не
тождественный physical coordinate. Независимая 1D c и все coupling coefficients
не меняются для согласования proxy.

Canonical seven-field ordering (u,w,v,Phi,psi,theta,c) сохранён. Заранее выбранные
observations: u/theta/c при x=.25, w/v при x=.5, Phi/psi при x=.25. Three 1D zeros
v=Phi=psi=0 отражают выбранную planar invariant subspace. Actual 3D out-of-plane
remainders не принудительно зануляются: nonlinear max|v|=1.36429e-6,
max|Phi|=1.56939e-5,max|psi|=1.64020e-6 на sampled domain. Они показаны с declared
scale relative to active displacement/rotation. Ни planarity stability, ни
Floquet вывод из их малости не делается.

Full-horizon linear/PCHIP w difference максимум 2.08656e-7 (linear), 2.10182e-7 (NL),
а correction/evolution difference 5.5562571392e-8. Это примерно 2.02% full evolving
model gap и сохраняется как observed interpolation sensitivity, не error bound.
Output cadence 5 не объявляется самостоятельной temporal accuracy проверкой.

### Actual energy, resource totals, figures and bounded conclusion

| Full-period energy evidence | 3D linear | 3D nonlinear |
|---|---:|---:|
| Final STATIC internal energy | 3.736789e-8 | 3.726886e-8 |
| Native initial DYNAMIC energy reference | 7.473578e-8 | 7.453772e-8 |
| Raw STATIC/DYNAMIC reference jump, % | +100 | +100 |
| Within-DYNAMIC native-reference relative drift | 1.3380471842e-7 | 1.3416026134e-7 |
| Max absolute independent/native kinetic difference | 6.8563192387e-14 | 6.8143461588e-14 |
| Difference divided by fixed own kinetic scale | 1.8410783640e-6 | 1.8346962269e-6 |
| Native external /damping work after release | 0 /0 | 0 /0 |

Raw energies не сдвигаются и не перенормируются. Independent mass-weighted K
использует прежнюю positive 14-point C3D10 quadrature и actual nodal velocities;
StVK internal reconstruction NOT_RUN, ENERGY_DIAGNOSTICS PARTIAL. Native small
drift — bookkeeping diagnostic, не certificate of physical energy continuity.
Own 1D energy drift/safety и all8 PARTIAL остаются отдельными результатами.

В full cases minimum sampled detF L/NL=.9936316326/.9936579599; max sampled
Green strain=.0070452460/.0070764944. Independent preload force imbalance
4.84236e-8/1.51072e-8, moment imbalance 3.12193e-9/7.98216e-9 проходят прежний 1e-5
gate. Это 14 volume points/elements в 401 native frames, не непрерывный supremum.

Всего ровно 6 sequential CCX jobs, 2 nonlinear 1D integrations ; 0 new Gmsh/modal/
static-only/eigen/BVP/root/symbolic jobs. Actual native time 12415.43545s,
recorded numerical program 15059.22507s <24000s. Per-case memory меньше 211MiB
при 4GiB ceiling; максимум individual native runtime 3070.336089s <5400s.
Полные input/output, accepted timestamps/increments, original и actual execution
snapshots, source/config/binary/runtime hashes, attempt ledger, histories,
NPZ/CSV и manifest сохранены в отдельном bundle c6256269eb8143ef.

Техническая [robustness figure](../../results/nlsp_nonlinear_dynamic_validation/c6256269eb8143ef/numerical_robustness_controls.pdf)
показывает C1 controls. Dissertational figures используют только actual data:
[seven-field observations](../../results/nlsp_nonlinear_dynamic_validation/c6256269eb8143ef/full_period_seven_fields.pdf),
[linear/NL movement and corrections](../../results/nlsp_nonlinear_dynamic_validation/c6256269eb8143ef/full_period_free_motion_and_correction.pdf),
[representative profiles](../../results/nlsp_nonlinear_dynamic_validation/c6256269eb8143ef/full_period_representative_profiles.pdf).
PNG versions и underlying arrays/`full_period_native_observations.csv`,
`full_period_observations.csv`, `robustness_summary.csv` сохранены.
Отдельный [scientific summary](nlsp_straight_rod_3d_fem_verification_summary.md)
даёт формулировку для диссертации, сохраняя различие numerical convergence,
physical model differences и отсутствующей experimental validation.

| Final FEM-3C status | Outcome |
|---|---|
| NLSP_FEM3C_SOURCE_PRESERVATION | PASS |
| NLSP_FEM3C_INPUT_PROTOCOL | PASS |
| NLSP_FEM3C_TEMPORAL_CONTROL | PASS |
| NLSP_FEM3C_SPATIAL_CONTROL | PASS |
| NLSP_FEM3C_ROBUSTNESS | PASS |
| NLSP_FEM3C_FULL_PERIOD_1D | PASS |
| NLSP_FEM3C_FULL_PERIOD_3D | PASS |
| NLSP_FEM3C_SEVEN_FIELD_RECOVERY | PASS |
| NLSP_FEM3C_ENERGY_DIAGNOSTICS | PARTIAL |
| NLSP_FEM3C_VERIFICATION_SUMMARY | PASS, scoped evidence complete |
| Full 1D all8 spatial /strict float64 | PARTIAL /PARTIAL |

Overall STRAIGHT_ROD_NONLINEAR_3D_FEM_VERIFICATION_COMPLETE_WITH_QUALIFICATIONS.
Это завершённое независимое limited bending verification, не
UNIVERSAL_NONLINEAR_MODEL_VALIDATION_PASS или agreement всех active fields.
В рассмотренном случае сохраняются близкие linear frequencies, sign/scale
static bending correction и numerically robust quarter-period evolving-w
comparison. Full-period results дают meaningful illustration с видимыми
model differences; их spatial/time convergence отдельно не доказана.
[NLSP-K17](../memory/knowledge.md#nlsp-k17) сохраняет итог, D16 — authorization.

```powershell
python scripts/analysis/resume_nlsp_nonlinear_dynamic_3d_fem.py --validation --compute
python scripts/analysis/resume_nlsp_nonlinear_dynamic_3d_fem.py --validation --report-only results/nlsp_nonlinear_dynamic_validation/c6256269eb8143ef
python scripts/analysis/resume_nlsp_nonlinear_dynamic_3d_fem.py --validation --plot-only results/nlsp_nonlinear_dynamic_validation/c6256269eb8143ef
```

Matching completed compute/report/plot не выполняют новых CCX/Gmsh/Radau/static/
eigen/BVP calls. V0, quartic action, variable mass, RHS/Jacobian, coefficients,
материал/геометрия/нагрузка/BC и saved IC не менялись; старые bundles/failed
ledgers не переписаны. Ни full-period fine/refined dynamics, ни другие amplitudes,
angular joints, out-of-plane perturbations/Floquet, nonlinear periodic orbits
или critical-amplitude search не выполнялись. LONG CLOSED, EB/RLB-KV
PAUSED_FOR_SUPERVISOR_DIRECTION, angular same-clamp UNAVAILABLE и все исторические
qualified/PARTIAL statuses сохраняются. После отчёта остановка.
