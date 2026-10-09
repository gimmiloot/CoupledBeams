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
