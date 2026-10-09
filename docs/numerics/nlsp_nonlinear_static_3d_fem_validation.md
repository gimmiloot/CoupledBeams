# NLSP FEM-2: nonlinear static 1D / 3D comparison

2026-10-09. Initial checkout: main, HEAD
`e10b23d7f6c27b4dae6cdf869704f8c58b2ada44`; initially clean working tree and index.
**PARTIAL: the 1D static preflight completed, but the first 3D linear job failed
while reading its input deck. No 3D static equilibrium was obtained.**

The concrete failure is a numeric serialization defect in the new static deck,
not an observed disagreement between physical models or nonlinear divergence.
The prescribed medium-level hard gate stopped the remaining five jobs. The
failed input, original execution code, stdout/stderr and manifest are retained;
no hidden retry or replacement of the failed evidence occurred.

## 1. Purpose of FEM-2

The authorized comparison concerns one monolithic fixed-fixed rectangular rod
under a uniform dead transverse body load. It compares linear and nonlinear
statics of the adopted four-field planar restriction with independent 3D
finite-strain elasticity. The quantities of interest are complete profiles and
`Delta w = w_NL - w_linear`, each relative to its own linear baseline. It is not
a modal calculation, nonlinear dynamics, inertia audit or universal validation
of the seven-field law. [User scope D11](../memory/decisions.md#nlsp-d11).

## 2. Frozen FEM-1R baseline

The [linear FEM-1/FEM-1R report](nlsp_linear_rectangular_3d_fem_validation.md)
and [readiness audit](nlsp_3d_fem_environment_readiness.md) remain historical.
The two source bundles and their complete artifact manifests are verified:

| Source | Bundle | Pinned manifest SHA256 |
|---|---|---|
| FEM-1 | `results/nlsp_linear_rectangular_3d_fem/4262efa427b03dad/` | `9f2d5139b84b2aa133b20d9a7cae806bac085c178fba506e087cae2941da1d90` |
| FEM-1R | `results/nlsp_linear_rectangular_3d_fem_refinement/63d44daae533389c/` | `bf42a0bbb1bbee8c64a18068848377bea477d8decd2990cf7f36bbedc832eb24` |

Saved medium/fine/refined meshes are loaded directly; no Gmsh call, new mesh,
repeated modal job or new torsional section reduction occurs. The existing
coefficients are read from the frozen FEM-1 reference. The existing action
archive `results/weakly_nonlinear_spatial_rod/a9cedd4b6de99295/result.json`, SHA256
`9a71ae8ba7f24c377a5f87d07c13389048b75f75d4906c659943952c4bc308df`, is reused
without a symbolic derivation.

## 3. Physics and constitutive comparison

| Quantity | Fixed value / contract |
|---|---|
| Geometry | L=1, b=0.20, h=0.10; A0=0.02 |
| Material | E=rho=1, nu=0.3; project kappa=5/6 |
| Physical axes | X along rod; thickness along Y, width along Z; positive w is negative global Y |
| 1D fields | independent `(u,w,theta,c)`; all four endpoint values zero |
| 3D support | all three translations zero on both complete end faces; free lateral faces |
| Internal connection | none: no joint, rigid insert, spring or additional MPC |

The [accepted reduced law](../theory/weakly_nonlinear_spatial_rod.md) supplies
`V_le4 = V2 + V3 + V4`; its coordinate gradient is cubic. It is not replaced by
a fitted stretching term, exact untruncated energy or von Karman functional.
The material, kappa, H and all independent fields remain unchanged. No slope
constraint or prescribed `c=-nu*u_s` is introduced.

Local CalculiX 2.22 documentation establishes a different 3D constitutive
contract: `*ELASTIC + NLGEOM` is St-Venant–Kirchhoff, Green–Lagrange strain with
internal second Piola–Kirchhoff stress (manual printed p255); printed stress S
is Cauchy (p600). FRD strain E is Green–Lagrange. These are not asserted identical
to V0. Nodal stress/strain output is extrapolated/averaged and is not an exact
pointwise maximum. Normalized E/rho/length units are retained; there is no Hz or
Lambda conversion and no specified laboratory material.

## 4. Load selection before FEM

The load was frozen from the independent linear fixed-fixed Timoshenko
compliance, before the first static job:

| Quantity | Selected value |
|---|---:|
| Target max w_linear / h | 0.05 |
| Target max w_linear | 0.005 |
| Acceleration g, global direction `(0,-1,0)` | 0.0014224751066856333 |
| q = rho A0 g | 2.844950213371267e-5 |
| Total load q L | 2.844950213371267e-5 |
| Linear max slope | 0.0149886793741 |
| Linear max rotation | 0.0136877470839 |
| Linear max shear strain | 0.00221906116643 |
| Linear bending surface strain | 0.00711237553343 |

The surface strain is below the predeclared diagnostic ceiling 0.01. The backup
w/h=0.03 was not needed or evaluated. This ceiling is a load-selection guide,
not a universal material-validity threshold. The 1D line load is positive along
w. CalculiX `*DLOAD, GRAV` uses the same fixed global acceleration and reference
mass (manual pp470–475 and local `rhs.f` / `e_c3d_rhs.f`); it is not follower
pressure. The load was not selected from FEM agreement and was not changed after
the failure.

## 5. Linear 1D equilibrium

The unchanged `PlanarGalerkin.linear_stiffness` is the Hessian of V2 at zero.
The load vector is the integral of the w basis times q. The assembled
`K0 a = f_q` solution is checked against the independent fixed-fixed Timoshenko
uniform-load profile, preserving independent section rotation and shear.
The MH block stays unforced in the linear problem: u=c=0. Essential endpoint
residuals are zero; no derivative condition is imposed.

## 6. Nonlinear 1D equilibrium

The same `PlanarGalerkin.potential(..., gradient=True, hessian=True)` supplies
energy, gradient and analytic Hessian for stationary `V_le4-q*integral(w)`.
Ten load increments continue the unloaded branch; each converges in four Newton
iterations. No subdivision, damping, branch search or changed physical energy
was needed. Quadrature is nq=2p+1; all four coefficient blocks remain independent.

| p / nq / DOF | Linear midspan w | Nonlinear midspan w | Delta w at midspan | Relative NL residual | Min tangent eigenvalue |
|---|---:|---:|---:|---:|---:|
| 48 / 97 / 188 | 0.004999999999999941 | 0.004991129920693116 | -8.870079306825e-6 | 5.08049e-14 | 0.368447216 |
| 64 / 129 / 252 | 0.004999999999999884 | 0.004991129920693126 | -8.870079306758e-6 | 7.13923e-14 | 0.368447216 |

The tangent is positive in the finite coefficient space; its condition number
is about 1.82e6 / 5.57e6. Its eigenvalues are static Hessian diagnostics, not new
eigenfrequencies. The nonlinear correction is negative: the midspan deflection
is reduced by about 0.177402% relative to the model's own linear value.

| p48 to p64 quantity | Absolute max difference | Relative max difference |
|---|---:|---:|
| Linear w | 6.07153e-17 | 1.21431e-14 |
| Nonlinear w | 1.30104e-17 | 2.60671e-15 |
| Delta w | 7.11237e-17 | 8.01838e-12 |
| Nonlinear u | 8.17916e-20 | 1.47521e-14 |
| Nonlinear theta | 2.35706e-16 | 1.72521e-14 |
| Nonlinear c | 1.63033e-18 | 8.83878e-14 |
| Endpoint reactions, largest relative entry | 1.14864e-16 | 8.07497e-12 |

Relative profile values use the saved p64 characteristic maximum, not local
zeros. L2 differences and full profiles are retained. These sampled comparisons
establish agreement of the two finite-dimensional static solutions, not exact
continuum error bounds. The strict independent float64 action/strong comparison
remains PARTIAL: nonlinear relative differences 1.16441e-11 / 4.41817e-11 exceed
the unchanged 2e-12 threshold. It is reported separately from the successfully
converged variational Newton residual; no old strict qualification is removed.

## 7. 3D static input and actual solver stop

One new focused [CLI](../../scripts/analysis/verify_nlsp_nonlinear_static_3d_fem.py)
reuses the existing action/helper, saved C3D10 meshes/audit and CalculiX runner.
This new static load/equilibrium/result contract is the Script Proliferation
Control reason for a distinct entry point rather than changing historical modal
commands. The [config](../../data/input/nlsp_nonlinear_static_3d_fem.json) fixes
all load, increment, output, recovery, source and resource policies.

The intended paired decks differ only by `NLGEOM`: same mesh/material/density,
full-face support, dead GRAV and static output. Static pseudo-time is load factor,
not a physical nonlinear trajectory. Increments are 0.1, total step1, minimum
1e-6, maximum0.1; field controls1e-8 are recorded. The attempted medium linear
job used the verified CalculiX 2.22 absolute executable and one-thread rule.
No PATH, installation, binary location or runtime cleanup changed.

The first job returned201 after **0.205941 s**, peak working set9,056,256 bytes.
The stdout explicitly reports `*ERROR reading *STATIC` and fatal `calinput`
termination. The new deck serialized the minimum increment with `.17g` as
`9.9999999999999995E-07` (22 characters), but local `statics.f` reads a 20-character
numeric field. The truncated exponent caused failure before stiffness assembly
or equilibrium. The actual bad card and source evidence are preserved. This is
a confirmed new I/O defect, not a nonlinear constitutive or convergence failure.
The local source reads each numeric entry using a 20-character field
(`statics.f`, lines181–199); a short source member and the actual failing card
are archived as independent evidence of the cause.

The current CLI now uses the existing bounded `ccx_float` representation for
all new static numeric cards. The minimum increment is serialized as a short
scientific value rather than an overlength `.17g` decimal. A regression checks
native field width, parseability and numerical round-trip error. The original
execution-code copy and actual failed input/log are untouched. A corrected
medium linear input and current-code snapshot are retained under
`remediation_preview/`, explicitly **NOT_RUN**. Its g serialization changes the
stored decimal by only +2.57794e-13 relative; the selected physical load remains
unchanged. No new preflight, corrected solver job or replacement bundle was
created. The formatting correction is verified locally; real execution of the
corrected medium pair has not been tested.

| Level | Linear static | Nonlinear static | Final load / displacement / reactions |
|---|---|---|---|
| medium | FAILED_PREPROCESSING; one actual job | NOT_RUN | not obtained |
| fine | NOT_RUN | NOT_RUN | not obtained |
| refined | NOT_RUN | NOT_RUN | not obtained |

The medium hard gate was enforced. No nonlinear, fine or refined job was
launched and no automatic retry occurred. Correcting numeric serialization
requires a separate execution test; passing a synthetic formatting test is not
successful 3D equilibrium evidence.

## 8. Source mesh quality

| Saved level | Target size | Nodes | C3D10 | End nodes left/right |
|---|---:|---:|---:|---:|
| medium | 0.0333333333333 | 5649 | 3120 | 119 / 119 |
| fine | 0.025 | 11553 | 6670 | 193 / 193 |
| refined | 0.020 | 20752 | 12687 | 279 / 279 |

These are the original verified meshes, not newly generated meshes. Source
hashes, geometry, material and restraints pass. The medium pre-job audit confirms
one connected C3D10 volume, correct bbox X=[0,1],Y=[-.05,.05],Z=[-.10,.10],
volume/mass0.02 and positive Jacobians. There was no deformed 3D mesh to assess.
Historical modal convergence PASS does not imply static-correction convergence.

## 9. Static section recovery contract

The implemented static reader distinguishes actual final STEP/increment/time
from intermediate frames and requires complete U/RF/S/E blocks. DAT U/RF
(seven significant digits) is preferred to FRD's float32/E12.5 vector precision;
FRD provides an independent consistency check. This parser has only synthetic
verification in this stage: the failed job produced no usable final 3D state.

The saved policy is quadrature-weighted original material slabs, local cubic
variation along X and affine transverse variation. Centroid displacement is
recovered relative to the original material section. A finite polar section
orientation is used for both linear and nonlinear states, with 41 versus81
sections as sensitivity control. Effective c is a thickness-strain diagnostic,
not a Cartesian FEM DOF. Synthetic rigid motion/affine examples check signs;
real section/recovery sensitivity remains NOT_RUN.

## 10. Linear deflections

| Model | Midspan w | Profile comparison |
|---|---:|---|
| 1D p64 | 0.004999999999999884 | p48/p64 checked |
| 3D medium/fine/refined | unavailable | NOT_RUN |

No physical 1D/3D linear discrepancy is reported without a 3D equilibrium.

## 11. Nonlinear deflections

| Model | Midspan w | Max u | Max c / effective c |
|---|---:|---:|---|
| 1D p64 quartic | 0.004991129920693126 | 5.544408787929e-6 | 1.844516529676e-5 |
| 3D NLGEOM | unavailable | unavailable | unavailable |

No nonexistent final fields, reactions or stresses are reconstructed from the
failed output. Missing intervals/increments are not filled or extrapolated.

## 12. Nonlinear corrections

| Quantity | 1D p64 | 3D | 1D/3D comparison |
|---|---:|---|---|
| Midspan Delta w | -8.870079306758e-6 | unavailable | NOT_RUN |
| Max abs Delta w | 8.870079306758e-6 | unavailable | NOT_RUN |
| Delta w / own linear midspan w | -0.177401586% | unavailable | NOT_RUN |
| Correction sign agreement | negative in 1D | unknown | cannot establish |
| Correction magnitude agreement | finite resolved1D signal | unknown | cannot establish |

The intended comparison retains signed/max/L2 differences and fixed displacement
scale h=0.10. The common correction denominator is the maximum absolute signal
across p64 and actual 3D results, never a local zero. With no 3D result this scale
is not used to invent a model-agreement percentage or physical validation.

## 13. Static mesh convergence

| Transition | Linear response change | Nonlinear response change | Delta w change |
|---|---|---|---|
| medium to fine | NOT_RUN | NOT_RUN | NOT_RUN |
| fine to refined | NOT_RUN | NOT_RUN | NOT_RUN |

The predeclared diagnostics are absolute/L2/max Delta w changes, a common signal
scale, sign stability, recovery/rounding uncertainty and equilibrium checks.
The old frequency threshold0.1% is not applied to static Delta w. No static
signal-to-mesh-change ratio, continuum correction or extrapolation can be inferred
from the incomplete series. No fourth static grid is added.

## 14. Reactions, strains and equilibrium

1D support forces/moments are obtained independently from endpoint derivatives
of V_le4, with endpoint signs -1 left/+1 right; they are not inferred from the
balance subsequently checked. Global Y reactions are positive against negative
Y gravity. The p64 nonlinear force/moment balances are at about1e-11 relative.

| 1D p64 observable | Linear | Nonlinear |
|---|---:|---:|
| Max abs u_s | 0 | 5.52995905e-5 |
| Max abs w_s | 0.0149886794 | 0.0149614200 |
| Max abs theta | 0.0136877471 | 0.0136623964 |
| Max abs c | 0 | 1.84451653e-5 |
| Max abs Gamma1 | 0 | 6.15870518e-5 |
| Max abs Gamma2 | 0.00221906117 | 0.00221906117 |
| Max bending surface strain | 0.00711237553 | 0.00710328220 |
| Support Y force, each end | 1.42247511e-5 | 1.42247511e-5 |
| Support X force, left/right | 0 / 0 | -1.10589612e-6 / +1.10589612e-6 |
| Global Z support couple, left/right | +2.37079184e-6 / -2.37079184e-6 | +2.36776073e-6 / -2.36776073e-6 |

Virtual-work/energy identity residuals are separately saved. min(1+c) in 1D is
0.9999815548; no NaN/Inf is observed. These are static checks, not dynamic mass
or state-dependent rotational inertia validation.

For future actual 3D output, local documentation pp558/562 and solver source show
`RF = support reaction + consistent applied body nodal load`. The implemented
recovery uses `R_support = DAT_RF_support - f_body_support`, not an unqualified
sum of RF. Reference C3D10 consistent weights are -V/20 for each vertex and V/5
for each midside; their independent14-point quadrature check on the saved medium
mesh agrees to2.86e-18 absolute. Full external resultant and moments are checked
independently; no reaction is constructed from balance. Actual3D reaction,
symmetry, out-of-plane, strain and deformation checks remain NOT_RUN.

## 15. Interpretation and limitations

The frozen load and adopted1D quartic static solution are reproducible and well
resolved in the p48/p64 comparison. The independent 3D part is blocked by one
localized deck-formatting defect. It supplies no evidence for agreement or
disagreement of V0 with 3D elasticity. The correction sign and magnitude cannot
yet be compared, and the absence of3D correction data is not
NONLINEAR_SIGNAL_UNRESOLVED: the static solve was never reached.

There is no basis from this stage to proceed directly to FEM-3. The next required
validation, if separately authorized, is execution of the corrected medium
linear/nonlinear pair and its parser/reaction gate before considering the
remaining saved grids. No future job is authorized by this report.

Solid-face restraint and1D c=0/section rotation clamp describe local deformation
differently. 3D section/end effects and StVK versus reduced quartic conventions
must remain explicit in any eventual discrepancy interpretation. No causal
mechanism is established here. The stage does not test inertia, dynamics of c,
frequency shifts, periodic orbits, angular joints, out-of-plane stability/Floquet,
damping or critical amplitude.

## 16. Reproducibility, statuses and stop

Primary attempted-job bundle:
`results/nlsp_nonlinear_static_3d_fem/6714bd9f2778e6d7/`.
Its own manifest and execution-code copy preserve the original serialized deck
and the observed failure. Historical sources are referenced by pinned hashes;
they are not regenerated because the new CLI hash changes. Matching cache/report/
plot execute zero new CCX/Gmsh/static/BVP/ODE/eigen/symbolic calls, as verified by
forbidden-call guards. Preflight, compute, report-only and plot-only replay the
same failed attempt even after the localized current-code fix; the attempt
ledger prevents an automatic retry from being triggered by a changed code hash.
A repaired current code path cannot retroactively certify or overwrite the
historical failure. Figure hashes are unchanged by cache replay.

The recorded primary numerical stage is **1.09647 s / 3600 s**: 0.88298 s
for the actual 1D preflight and 0.21349 s for the failed 3D stage. Native CCX time
is 0.20594 s, peak working set **8.63672 MiB**. The separate development 1D
fragment/check runs used 0.87045+0.69032 s and are recorded explicitly;
these together with the primary stage total about **2.65724 s**, before the
lightweight tests and cached postprocessing. No additional real FEM job is
hidden in tests or report commands.

Two PDF+PNG figures show only obtained 1D evidence:
[static profiles](../../results/nlsp_nonlinear_static_3d_fem/6714bd9f2778e6d7/figures/one_d_static_preflight_profiles.pdf)
and [nonlinear correction](../../results/nlsp_nonlinear_static_3d_fem/6714bd9f2778e6d7/figures/one_d_static_preflight_correction.pdf).
They identify p48/p64 and the missing 3D equilibrium explicitly. No 3D curve,
mesh-convergence plot or extrapolated comparison is fabricated.

```powershell
python scripts/analysis/verify_nlsp_nonlinear_static_3d_fem.py --preflight
python scripts/analysis/verify_nlsp_nonlinear_static_3d_fem.py --run-fem
python scripts/analysis/verify_nlsp_nonlinear_static_3d_fem.py --report-only results/nlsp_nonlinear_static_3d_fem/6714bd9f2778e6d7
python scripts/analysis/verify_nlsp_nonlinear_static_3d_fem.py --plot-only results/nlsp_nonlinear_static_3d_fem/6714bd9f2778e6d7
```

The run-fem command is the execution entry point, not permission for another
attempt after this hard-gate stop. Report/plot are read-only with respect to
science. No plotted3D curve or fake full-comparison figure is produced.

| Status | Outcome |
|---|---|
| NLSP_FEM2_LOAD_PREFLIGHT | PASS |
| NLSP_FEM2_1D_LINEAR_STATIC | PASS |
| NLSP_FEM2_1D_NONLINEAR_STATIC | PASS |
| NLSP_FEM2_3D_LINEAR_STATIC | FAIL |
| NLSP_FEM2_3D_NONLINEAR_STATIC | NOT_RUN |
| NLSP_FEM2_STATIC_SECTION_RECOVERY | NOT_RUN |
| NLSP_FEM2_STATIC_EQUILIBRIUM | NOT_RUN |
| NLSP_FEM2_NONLINEAR_CORRECTION_MESH_CHECK | NOT_RUN |
| NLSP_FEM2_1D_3D_COMPARISON | NOT_RUN |
| Overall | PARTIAL |

[Tests](../../tests/test_nlsp_nonlinear_static_3d_fem.py) use synthetic static
inputs/parsers/affine recovery/reaction fixtures and saved evidence, never extra
FEM jobs. **53 targeted tests PASS (2.79 s)**, including the confirmed native
numeric-field regression and failed-attempt cache behavior. Verification/counters
are retained separately in the new evidence. **16 selected historical
regressions PASS (1.29 s)**, without fresh roots or real FEM jobs; historical tests
are unchanged. JUnit reports are retained as `targeted_tests.xml` and
`historical_regressions.xml`. Source/history/link preservation and whitespace
checks are recorded separately in final verification.
[NLSP-K12](../memory/knowledge.md#nlsp-k12) records this partial result.
FEM-1R PASS remains confined to its linear scope; physical sanity stays
DIAGNOSTIC_COMPLETE_WITH_QUALIFICATIONS, prepared planar strict PARTIAL,
LONG CLOSED, EB/RLB-KV PAUSED_FOR_SUPERVISOR_DIRECTION and angular same-clamp
reference UNAVAILABLE. V0/cubic equations, physical coefficients, basis, old
solvers/results/manifests and historical D/K are unchanged. One static CCX
attempt, zero new Gmsh/modal/ODE calls; no nonlinear 3D job or automatic further
study. Stop after the report.


<a id="fem-2r-controlled-continuation-after-input-serialization-failure"></a>

## FEM-2R - controlled continuation after input serialization failure

2026-10-09. Initial checkout: main, HEAD
`a650b77db3ebcb73bce1c57bce9bc98472950c9b`, clean working tree and index.
This section extends the historical FEM-2 result above; it does not replace its
PARTIAL status or its failed medium job. [Authorization D12](../memory/decisions.md#nlsp-d12),
[result K13](../memory/knowledge.md#nlsp-k13).

**FEM2_DIAGNOSTIC_COMPLETE_WITH_QUALIFICATIONS.** All six newly authorized
static jobs completed and passed their output/equilibrium gates. The nonlinear
bending correction has the same negative sign in both models. On refined 3D,
Delta w at midspan is -9.1626824679e-6 versus saved 1D -8.8700793068e-6;
the full-profile max difference is 3.19353% of the common correction scale.
The signal exceeds the observed mesh/recovery/printed-rounding measures, with
no claim of an exact 3D reference or universal physical validation of V0.

### Frozen sources and serialization gate

The separate bundle is `results/nlsp_nonlinear_static_3d_fem_resume/210b74b8b166997c/`.
Its identity pins the historical failed parent6714bd9f2778e6d7 manifest SHA256
`e0a972bf852c64658202bc5e385fabbede8ea5f5396248603e0ed6fe76d3d62c`, all31 parent
artifact hashes, the two FEM-1/FEM-1R source manifests and their94/41 artifacts,
saved p48/p64 NPZ profiles, frozen load/material/geometry, source mesh hashes,
corrected generator, solver/runtime DLLs, code/environment and explicit
`explicit_user_FEM2R_2026_10_09` authorization. Sources are referenced directly;
large historical arrays are not copied or recomputed. The old failed input,
code/logs/manifest and original failed-attempt guard are unchanged.

The preserved 1D profiles independently reproduce w_linear=.005,
w_NL=.004991129920693126 and Delta w=-8.8700793068e-6 at x=.5. Existing static
residuals, energies, reactions and p48/p64 comparison are reused; this stage
performs zero new 1D equilibria, BVP/eigen/section reduction or symbolic work.
The strict float64 strong/weak relative threshold2e-12 remains PARTIAL.

Before CCX, original failed input, corrected preview and freshly serialized
inputs were compared. All six linear/NL cards pass finite parseability and
native20-character checks for each numeric field in STATIC/CONTROLS/DLOAD/
ELASTIC/DENSITY. Maximum new token width is18; the failed minimum-increment
token was22. New medium input matches the preview except its mesh include path;
linear/NL pairs differ only by NLGEOM. The same frozen g is represented within
+2.57794e-13 relative rounding, not replaced by a fitted load.
No additional generator or parser correction was needed in this continuation.

### Actual static execution and output recovery

Same L=1,b=.20,h=.10; E=rho=1,nu=.3,kappa=5/6; g=.0014224751066856333,
q=F_total=2.844950213371267e-5, global load direction(0,-1,0). Full end faces
remain fixed, lateral faces free, with no internal joint/MPC/contact/spring.
No Gmsh or modal job ran; the saved medium/fine/refined C3D10 meshes were reused.
All source mesh quality/material/clamp contracts passed before execution.

CalculiX2.22 binary is unchanged:
`results/_smoke/3d_fem_environment_check/calculix_2p22/CalculiX-2.22.0-win-x64/bin/ccx.exe`.
One thread, job timeout1200s and4GiB ceiling; jobs execute sequentially in the
explicit ledger order. Medium linear and NL both passed before fine/refined.
Each linear solve has one accepted full-load increment; each NLGEOM solve has
10 accepted increments of0.1, two iterations each, no cutback/retry/warning.
STA final step time is1, i.e. full load factor1 under the preserved ramp.
This load-step pseudo-time is not a time-dynamic trajectory.

| New case | Nodes / C3D10 | CCX seconds | Peak working set, MiB | Final load | Output / equilibrium |
|---|---:|---:|---:|---:|---|
| medium linear | 5649 / 3120 | 1.2113 | 69.70 | 1 | PASS |
| medium NLGEOM | 5649 / 3120 | 10.1626 | 74.21 | 1 | PASS |
| fine linear | 11553 / 6670 | 2.8220 | 150.36 | 1 | PASS |
| fine NLGEOM | 11553 / 6670 | 27.0664 | 159.46 | 1 | PASS |
| refined linear | 20752 / 12687 | 6.5445 | 331.57 | 1 | PASS |
| refined NLGEOM | 20752 / 12687 | 80.6070 | 347.02 | 1 | PASS |

The numerical stage, including mesh audits and recovery, is138.31963s/3600s.
All return codes0 are supported independently by completion markers, complete
final DAT/FRD U/RF/S/E blocks, STA increments, finite nodal fields, zero fixed-face
displacements and checked force/moment balances. No parser repair or solver
repetition occurred. Exact decks, source include hashes, native command, logs,
DAT/FRD/STA, nodal arrays, sections and per-case diagnostics are retained at once.

An independent read-only audit parses INP/DAT/FRD/STA and integrates the body
load without importing the production static parser. It confirms six complete
cases,33 accepted increments (3linear+30NL),60 nonlinear Newton iterations and
zero cutbacks. Native stdout prints average/residual forces with only six
fixed decimal places; printed zero is not proof of zero numerical residual.
Accepted correction/increment ratios range from about1.2034e-4 at the first
increment to1.35-1.37e-7 at the last, so nominal control can=1e-8 cannot be
reported as a demonstrated bound for every accepted state. Local2.22 convergence
source permits alternative acceptance branches; coarse stdout cannot identify
the exact branch. Controls are unchanged. Independent final RF balances pass,
but no complete native Newton-error bound or tighter-control test is claimed.
Small recovered asymmetry/transverse fields persist on the unstructured mesh;
perfect symmetry/planarity is not asserted.

Recovery uses the unchanged weighted original material-section policy, cubic X
variation/transverse affine fit, finite polar orientation for both linear/NL
states and41-versus81-section check. Positive w is negative global Y; the same
undeformed material coordinate x is used, with no surface-node substitution,
profile/phase/amplitude alignment or altered bounding box. DAT displacement
(seven significant digits) is primary; FRD is independent rounded evidence.

### Deflection and correction comparison

| Model / mesh | Linear midspan w | NL midspan w | Signed Delta w | Relative NL effect, % | Status |
|---|---:|---:|---:|---:|---|
| saved 1D p64 | .004999999999999884 | .004991129920693126 | -8.8700793068e-6 | -.177401586 | existing equilibrium PASS; strict PARTIAL |
| 3D medium | .004821160849899094 | .004812087371842098 | -9.0734780570e-6 | -.188201106 | PASS |
| 3D fine | .004831278709155392 | .004822130272658330 | -9.1484364971e-6 | -.189358492 | PASS |
| 3D refined | .004837259676401260 | .004828096993933395 | -9.1626824679e-6 | -.189418867 | PASS |

Relative nonlinear effects use each model's own linear midspan w. Profile
comparisons use the predeclared common801-point original-x grid, composite
trapezoid L2 and sampled maxima; no continuous supremum is proved. The shared
correction scale is max absolute Delta w across saved p64/all three3D results,
9.1626824679e-6. Translation differences also retain the fixed scale h=.10;
angles/strain diagnostics use fixed scale1. No denominator is a local zero.

| Refined 1D-minus-3D w quantity | Midspan signed difference | Absolute profile max | L2 | Relative max, % |
|---|---:|---:|---:|---:|
| Linear displacement | +1.6274032360e-4 | 1.6287136025e-4 | 1.1526475242e-4 | 3.25743 of common linear scale |
| Total nonlinear displacement | +1.6303292676e-4 | 1.6316227131e-4 | 1.1543881626e-4 | 3.26904 of common NL scale |
| Nonlinear correction | +2.9260316111e-7 | 2.9261288903e-7 | 1.7764915170e-7 | 3.19353 of common correction scale |

The refined3D correction max/L2 are9.1626824679e-6 /5.6997874068e-6;
saved1D values are8.8700793068e-6 /5.5222537356e-6. Both are negative in the
interior and reduce deflection relative to their own linear response. The
small nonlinear signal is not mistaken for the much larger linear offset.

For context, saved full-profile comparisons include all recovered fields.
Refined nonlinear-correction max differences are9.86842e-7 for u (17.7989% of
the common u-correction scale) and4.47796e-7 for theta (1.73663%). Effective c
is only a finite-section thickness-contraction diagnostic: its comparison is
not identity of independent M-H and solid DOFs, and its finite polar extraction
has a nonzero linear baseline. Its correction discrepancy is82.9202% of the
diagnostic common scale; close Delta w must not be reported as agreement of all
four fields or validation of c/V0 coefficients. Solid-face/clamp and section
reduction distinctions are retained. Small v/Phi/psi diagnostics are saved,
without claiming new spatial dynamics or classifying modes.

### Static mesh and small-signal qualification

The frequency criterion0.1% is not applied to Delta w. Instead the original
signal rule compares signal with observed successive-grid, recovery and printed
rounding measures and retains sign/decreasing-change checks.

| Transition | Max change w_linear | Max change w_NL | Max change Delta w | L2 change Delta w | Delta w change / common signal, % |
|---|---:|---:|---:|---:|---:|
| medium to fine | 1.1469171206e-5 | 1.1391400023e-5 | 7.7790574643e-8 | 4.9220784763e-8 | .848993 |
| fine to refined | 5.9866240008e-6 | 5.9723940965e-6 | 1.4245970802e-8 | 7.7699733409e-9 | .155478 |

Absolute and L2 correction changes decrease; correction sign is stable. No
extrapolated value is substituted for the actual refined result and no static
fourth grid is introduced. The 1D/3D max correction differences are2.03416e-7,
2.80317e-7,2.92613e-7 for medium/fine/refined, respectively; refinement does not
artificially drive them to agreement.

| Refined w-correction evidence | Absolute magnitude |
|---|---:|
| Nonlinear signal | 9.1626824679e-6 |
| Last fine/refined profile change | 1.4245970802e-8 |
| Paired correction recovery change,41 vs81sections | 3.1405842805e-9 |
| Conservative summed DAT nodal printed-rounding allowance | 4.906427e-9 |
| Signal / largest observed measure | 643.177 |

The saved outcome is SIGNAL_EXCEEDS_OBSERVED_UNCERTAINTY. The last change is
not a rigorous continuum-error bound; the rounding allowance is a pointwise
formatting estimate, not a full propagated recovery/solver error theorem.
DAT/FRD nodal displacement differences are at most5e-9; equilibrium convergence
is checked separately. No numerical green threshold was invented after seeing
the result. The observation supports a resolved negative bending correction
and a quantified finite-grid comparison, not exact continuum truth.

### Reactions, equilibrium and strains

Actual support forces are independently recovered as
`R_support=DAT_RF_support-consistent_reference_bodyload_support`.
Bodyload is integrated from saved reference C3D10 geometry, not deduced from
reaction balance. Nonlinear moment balance uses current nodal coordinates and
the corresponding dead-load moment. Raw RF is not treated as pure support
reaction; the free-node RF/bodyload residual is also retained.

| Refined support observable | Linear | NLGEOM |
|---|---:|---:|
| Left R_X | +5.9510850e-11 | -1.1373996719e-6 |
| Right R_X | -5.9422700e-11 | +1.1373997021e-6 |
| Left R_Y | 1.4224656856e-5 | 1.4224656517e-5 |
| Right R_Y | 1.4224844976e-5 | 1.4224845426e-5 |
| Left M_Z about face centroid | +2.3882004421e-6 | +2.3851623695e-6 |
| Right M_Z about face centroid | -2.3882944827e-6 | -2.3852567457e-6 |
| Force imbalance / total force | 1.18159e-8 | 7.04492e-9 |
| Moment imbalance / total force times L | 9.86754e-9 | 2.15632e-9 |
| Maximum fixed-face displacement | 0 | 0 |
| Maximum recovered nodal strain component | .00798158 | .00790246 |
| Maximum recovered nodal stress component | .0107444 | .0105592 |

Across all six cases the maximum relative force/moment imbalance is
4.84236e-8 /1.36216e-8, below the unchanged1e-5 equilibrium gate. All fields are
finite. Strain is infinitesimal for linear and Green-Lagrange for NLGEOM; printed
stress is Cauchy. FRD S/E are extrapolated/averaged nodal diagnostics, not exact
local maxima; their mesh variation near full-face clamps remains qualified.
The unchanged1D reaction/strain diagnostics above remain a separate reduction.

### Interpretation, replay and stop

For this one load/geometry the adopted reduced quartic model and StVK3D both
predict a small reduction of transverse deflection. The resolved correction
profiles differ by about3.19% on the declared common max scale, while total
linear/NL deflections retain about3.26% offsets. This is useful independent
physical evidence for the sign and order of the static bending response;
it is not PHYSICAL_NONLINEAR_VALIDATION_PASS for all V0 coefficients.
Section/shear/clamp/constitutive differences are possible contributors, without
unique causal separation. Hardening is not attributed solely to axial stretching.
Variable inertia, nonlinear dynamics/Floquet, damping, angular joints and
amplitude dependence are not checked by this static task.

[Continuation CLI](../../scripts/analysis/resume_nlsp_nonlinear_static_3d_fem.py)
and [config](../../data/input/nlsp_nonlinear_static_3d_fem_resume.json) implement
only the explicit immutable-parent/authorization/ledger orchestration contract,
reusing original FEM-2 functions. The old failed guard remains in place. Cache
identity includes parent manifest/source profiles/meshes, corrected generator,
code/solver/runtime/environment and frozen numerical/recovery settings. A new
solver failure cannot acquire a hidden retry merely by changing code hash.

```powershell
python scripts/analysis/resume_nlsp_nonlinear_static_3d_fem.py --check-source
python scripts/analysis/resume_nlsp_nonlinear_static_3d_fem.py --run-fem
python scripts/analysis/resume_nlsp_nonlinear_static_3d_fem.py --report-only results/nlsp_nonlinear_static_3d_fem_resume/210b74b8b166997c
python scripts/analysis/resume_nlsp_nonlinear_static_3d_fem.py --plot-only results/nlsp_nonlinear_static_3d_fem_resume/210b74b8b166997c
```

Three PDF+PNG figures preserve full profiles and separate correction/numerical
mesh diagnostics: [linear/NL profiles](../../results/nlsp_nonlinear_static_3d_fem_resume/210b74b8b166997c/figures/linear_and_nonlinear_static_profiles.pdf),
[correction profiles](../../results/nlsp_nonlinear_static_3d_fem_resume/210b74b8b166997c/figures/static_nonlinear_corrections.pdf),
[mesh changes/reactions](../../results/nlsp_nonlinear_static_3d_fem_resume/210b74b8b166997c/figures/static_mesh_convergence_and_reactions.pdf).

| FEM-2R status | Outcome |
|---|---|
| NLSP_FEM2R_SOURCE_PRESERVATION | PASS |
| NLSP_FEM2R_INPUT_SERIALIZATION | PASS |
| NLSP_FEM2R_MEDIUM_LINEAR | PASS |
| NLSP_FEM2R_MEDIUM_NONLINEAR | PASS |
| NLSP_FEM2R_FINE_PAIR | PASS |
| NLSP_FEM2R_REFINED_PAIR | PASS |
| NLSP_FEM2R_STATIC_OUTPUT_RECOVERY | PASS |
| NLSP_FEM2R_EQUILIBRIUM | PASS |
| NLSP_FEM2R_NONLINEAR_SIGNAL_RESOLUTION | PASS |
| NLSP_FEM2R_1D_3D_COMPARISON | PASS |
| Overall | FEM2_DIAGNOSTIC_COMPLETE_WITH_QUALIFICATIONS |

A PASS for comparison means actual complete, numerically qualified evidence,
not automatic exact physical agreement or universal applicability. Historical
FEM-2 PARTIAL/failed attempt and strict float64 PARTIAL are not raised.

V0/cubic equations, linear coefficients, material/load, Shen basis and BC are
unchanged. Exactly six new static CCX jobs; zero new meshes, modal/1D/BVP/ODE or
symbolic computations. LONG CLOSED, EB/RLB-KV PAUSED_FOR_SUPERVISOR_DIRECTION,
angular same-clamp reference UNAVAILABLE, prepared strict PARTIAL and physical
sanity DIAGNOSTIC_COMPLETE_WITH_QUALIFICATIONS remain. Stop after the report;
FEM-3, other amplitudes/geometries and angular/Floquet studies are not authorized
or launched automatically. Test/cache/link/preservation evidence is recorded in
the continuation bundle and the final verification paragraph below.


Final targeted verification:75 new continuation tests PASS and42 selected
historical serialization/I/O/recovery/failure-cache regressions PASS (117total),
using synthetic fixtures and saved evidence with no fresh scientific solver.
The first historical-test invocation had32pass/10setup errors solely from pytest
default temporary-directory ACL; both that XML and the42/42 successful fresh-
basetemp rerun are retained, with no code/threshold changes to make tests pass.
Artifacts: `targeted_tests.xml`, `historical_targeted_tests.xml` and
`historical_targeted_tests_initial_env_error.xml` in the continuation bundle.
`cache_checks.json` confirms all four source/cached-compute/report/plot routes
with numerical entry points forbidden: zero new CCX/Gmsh/1D/eigen/BVP/ODE/symbolic
calls and identical figure hashes. Historical report/D/K exact byte-prefix
preservation and unique new anchors pass; `git diff --check` passes. Final
source/hash/link/index/HEAD verification is recorded separately in the bundle.
