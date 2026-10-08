# Weakly nonlinear planar physical sanity checks

Date: 2026-10-08. Initial main HEAD: `5878f5b34d8ff04aae926c7a3e4be959cc78a4b6`.

**DIAGNOSTIC_COMPLETE_WITH_QUALIFICATIONS.** The saved large-amplitude motion
and one new half-amplitude motion exhibit the expected small-amplitude patterns
on one fixed linear period. These are internal consistency checks of the adopted
reduced constitutive law V0. Its physical accuracy against experiment or 2D/3D
elasticity has not been established.

## Scope and immutable sources

The model/action/kinetic energy, cubic residuals, coefficients, Radau RHS and
analytic Jacobian are unchanged. Four independent fields `(u,w,theta,c)` use
Shen-Legendre functions with essential values zero at both ends; no slope BC
or dynamic endpoint-jet constraint is introduced. The G20 geometry is
E=rho=1, nu=.3, b=.20, h0=.05, L=1, kappa=5/6. Absolute numbers below are
in this normalized system, not predictions for a specified laboratory material.

Canonical sources are the [frozen spatial model](weakly_nonlinear_spatial_rod.md),
[generated expansion](weakly_nonlinear_spatial_rod_expansion_generated.md),
[second-order response](planar_second_order_axial_response.md),
[prepared-state and one-period report](planar_prepared_initial_state.md#prepared-one-period-feasibility),
and [historical time pilot](weakly_nonlinear_planar_time_pilot.md).
The older linear equations in `equations.tex` and analytic production branches
are untouched; this diagnostic uses the separate audited NLSP action archive.
No literature search or new material assumption is made.

Historical bundles were read under their OWN manifests, not the new CLI hash:

- `results/planar_prepared_one_T1/795dcb14d3cd3a55/`: existing p64 tight .05 history,
  q/v, actual95057 timestamps, snapshots and existing p/time uncertainty.
- `results/planar_prepared_initial_state/5ea8d41faf8ede54/`: frozen common evaluator,
  analytical bending pair, quintic Theta3 and saved p96 Legendre stat/harm profiles.
  Their sum is byte-identical to the common U_star/C_star coefficients.
- `results/planar_second_order_axial_response/b3ea4eb6ac95d6e1/`: provenance of
  those leading profiles; no old BVP/eigensystem/million-point history recomputed.
- `results/planar_prepared_feasibility/284a4039177391d1/`: preserved strict
  qualification and numerical policy. Float64 relative strong/weak gate2e-12
  remains unpassed; its previous independent precision evidence is reused.

The action archive `weakly_nonlinear_spatial_rod/a9cedd4b6de99295/result.json`
has SHA256 `9a71ae8ba7f24c377a5f87d07c13389048b75f75d4906c659943952c4bc308df`.
It is checked against the historical provenance before deserializing the frozen
polynomials. Frozen model/helper/runner hashes also match the one-T1 manifest.

## Two amplitudes and initial state

With epsilon=A/h0, both motions use the SAME fixed functions:

\[
 u_0=\epsilon^2U_*,\quad c_0=\epsilon^2C_*,\quad
 w_0=\epsilon W,\quad \theta_0=\epsilon\Theta+\epsilon^3\Theta_3,
 \qquad \dot q_0=0.
\]

The new case has epsilon=.025, A=.00125, p64/nq129/252 independent coordinates;
the source has epsilon=.05, A=.0025. The common projection rule remains
`common_endpoint_constrained_L2`, only for INITIAL coefficients. It is not a
uniform half-scaling of the old q0: axial/contraction terms scale by1/4 and
the known cubic rotation correction by1/8. Theta3 and U_star/C_star are not refit.

The new preparation passes the original1e-6 initial policy: worst profile
L2/max relative error3.95125e-9 and fixed-scaled formal endpoint compatibility
coefficient2.93377e-8. This does not remove the historical strict strong/weak
qualification. Execution is **EXPLORATORY_NOT_CERTIFIED**, `admitted=False`.

omega1=.3174742907880648 and T1=19.791162590151373 are fixed by the original
linear bending pair. Both histories use the same95057 actual timestamps over
0...T1. Sampling maxima are not certified continuous suprema. Radau remains tight:
rtol1e-10, max_step=.0036046964938503193, and the full previous componentwise
atol prescription evaluated at the new amplitude (half of the old dimensional
vector). No interpolation, phase alignment, period or amplitude fitting is used.

## Comparison with leading response

The saved stat/harm profiles are evaluated directly at each actual timestamp:

\[
 (u,c)_{\rm lead}=\epsilon^2[(U,C)_{\rm stat}+
                              (U,C)_{\rm harm}\cos(2\omega_1t)],
 \quad (w,\theta)_{\rm lead}=\epsilon(W,\Theta)\cos(\omega_1t).
\]

Time derivatives use -2omega1 sin(2omega1t) and -omega1 sin(omega1t), respectively.
Theta's known epsilon^3 Theta3 at t=0 is not included in the linear reference;
its difference is a prescribed higher-order contribution, not integrator error.

Physical L2 norms use Gauss100 on [0,L]; maxima use the same spatial samples.
One fixed diagnostic scale per amplitude and component is used over ALL time:
q scales `(epsilon^2 h0,epsilon h0,epsilon h0/L,epsilon^2)` and velocity scales
multiply these by `(2omega1,omega1,omega1,2omega1)`. These asymptotic reporting
scales do not replace the historical trajectory acceptance denominators/gates.
The own-characteristic percentage uses the full-history maximum of that field,
never an instantaneous zero. The following L2 values use L=1.

| Component | .05 absolute L2 | .05 absolute max | .05 own max % | .05 fixed-scaled max | .025 own max % | Reduction of fixed-scaled deviation |
|---|---:|---:|---:|---:|---:|---:|
| q_u | 3.192898e-09 | 5.292747e-09 | 0.38615 | 4.234198e-05 | 0.09657 | 3.998929 |
| q_w | 5.047084e-06 | 8.130011e-06 | 0.32537 | 3.252004e-03 | 0.08140 | 3.996957 |
| q_theta | 1.781235e-05 | 2.497886e-05 | 0.33383 | 9.991545e-03 | 0.08352 | 3.997152 |
| q_c | 1.574177e-08 | 1.637445e-08 | 0.35580 | 6.549780e-06 | 0.08900 | 3.997664 |
| velocity_u | 2.623489e-09 | 4.436355e-09 | 1.01314 | 5.589562e-05 | 0.25439 | 3.994869 |
| velocity_w | 2.129728e-06 | 3.604140e-06 | 0.45398 | 4.541017e-03 | 0.11352 | 4.001373 |
| velocity_theta | 2.089334e-05 | 2.925334e-05 | 1.22111 | 3.685759e-02 | 0.30718 | 4.000731 |
| velocity_c | 1.591959e-08 | 3.012138e-08 | 2.05230 | 1.897564e-05 | 0.51577 | 3.996342 |

The leading approximation describes displacement u/c to about0.36-0.39% of their
large-case maxima; velocity deviations are about1.01% and2.05%. Their normalized
remainders decrease approximately fourfold when epsilon is halved. This is
consistent with higher-order terms, not a fitted exponent or exact-error claim.
The leading p96 profiles are numerical reference profiles, not continuum truth.

For context, unchanged large-case spatial and temporal uncertainties are:

| Component | saved p48-p64 absolute max | saved tight-extra absolute max | spatial difference / leading deviation |
|---|---:|---:|---:|
| q_u | 1.233158e-12 | 7.520141e-16 | 2.329901e-04 |
| q_w | 2.605483e-10 | 6.071137e-14 | 3.204772e-05 |
| q_theta | 7.636980e-09 | 3.228225e-12 | 3.057377e-04 |
| q_c | 1.419387e-11 | 3.668439e-14 | 8.668300e-04 |
| velocity_u | 6.248847e-11 | 2.101604e-13 | 1.408554e-02 |
| velocity_w | 1.313503e-08 | 1.507356e-11 | 3.644428e-03 |
| velocity_theta | 3.809423e-07 | 4.053737e-10 | 1.302218e-02 |
| velocity_c | 1.028686e-09 | 3.673507e-12 | 3.415137e-02 |

Large-case spatial comparison remains7/8 PARTIAL (theta_t max criterion failed),
temporal8/8 PASS. The new half case has NO independent p/time control. The tables
show measured differences; they cannot fully separate higher-order physics from
all numerical uncertainty or certify V0's physical accuracy.

## Amplitude structure and sign symmetry

Normalize u,c and their velocities by epsilon^2, w,theta and their velocities by
epsilon. The comparison uses the SAME physical times and profiles. Fixed scales
for these normalized quantities are `(h0,h0,h0/L,1)` with the above frequency
multipliers for velocities. Ratios below are full-history characteristic maxima,
not ratios taken near instantaneous zeros.

| Component | expected large/half ratio | measured ratio | normalized difference L2 | normalized difference max | fixed-scaled normalized max |
|---|---:|---:|---:|---:|---:|
| q_u | 4 | 4.000223 | 9.577672e-07 | 1.587683e-06 | 3.175365e-05 |
| q_w | 2 | 2.000013 | 7.568759e-05 | 1.219192e-04 | 2.438385e-03 |
| q_theta | 2 | 2.000084 | 2.671229e-04 | 3.745940e-04 | 7.491880e-03 |
| q_c | 4 | 4.000000 | 4.721191e-06 | 4.911379e-06 | 4.911379e-06 |
| velocity_u | 4 | 4.012291 | 7.874225e-07 | 1.330337e-06 | 4.190377e-05 |
| velocity_w | 2 | 2.001175 | 3.195044e-05 | 5.406829e-05 | 3.406152e-03 |
| velocity_theta | 2 | 2.012859 | 3.134159e-04 | 4.388269e-04 | 2.764488e-02 |
| velocity_c | 4 | 4.017350 | 4.776177e-06 | 9.033659e-06 | 1.422739e-05 |

Several displacement maxima occur at t=0, so their2/4 ratios partly express the
prescribed IC. Full normalized histories, velocity ratios and leading-deviation
reductions provide additional dynamic evidence. No new physical ratio tolerance
or claim of an exact general amplitude law is introduced.

For sign reversal of the selected bending amplitude, the local parity is
P=diag(+1,-1,-1,+1). The quartic kinetic/potential actions are invariant;
coordinate residuals transform with P. This includes all inertial terms.
Exact Fraction polynomial checks and nine saved-state RHS checks give
`RHS(Py)=P RHS(y)`; the largest observed relative/absolute discrepancy is0.
The common IC evaluator has the same parity, including the odd cubic Theta3;
the linear constrained projection preserves it. No negative-amplitude ODE ran.
This is bending reflection/sign symmetry, distinct from global rotational
invariance of a Taylor-truncated action.

## Formal classical bulk stretching benchmark

This is a controlled reduction for a benchmark, not a new production closure.
First impose no shear IN the reduced action and choose positive director
alignment: Gamma2=0 implies Gamma1=|r_s|-1. A constrained shear reaction must not
be discarded by blindly substituting Gamma2=0 into every independent theta PDE.

Formally neglect contraction gradient/inertia in the BULK normal energy
C/2(Gamma1^2+2nu Gamma1 c+c^2). Local stationarity gives c=-nu Gamma1 and hence
C(1-nu^2)Gamma1^2/2=EA Gamma1^2/2. For von Karman ordering u_s=O(delta^2),
w_s=O(delta), expand sqrt((1+u_s)^2+w_s^2)-1:

\[
 \Gamma_1=u_s+\tfrac12w_s^2+O(\delta^4).
\]

Quasistatic axial equilibrium makes N=EA(u_s+w_s^2/2) constant.
Integrating u_s with u(0)=u(L)=0 fixes

\[
 N=\frac{EA}{2L}\int_0^L w_s^2\,ds,\qquad
 V_{\rm stretch}=\frac{L N^2}{2EA}
   =\frac{EA}{8L}\left(\int_0^L w_s^2\,ds\right)^2.
\]

Its variation is N integral(w_s delta w_s ds), giving the potential gradient
-N w_ss when delta w vanishes at the endpoints. The coefficient1/8 and positive
sign follow from elimination, and the contribution is even in w and O(A^4).
With w dimension length, integral(w_s^2ds) has dimension length, N is force
and Vstretch is energy. Exact rational checks verify these operations.

In the finite M-H problem c=0 at both clamps. The relaxed relation c=-nu Gamma1
need not satisfy that BC. Therefore this is NOT a uniform exact limit/elimination
of the same finite four-field clamp problem; possible boundary layers are not
analyzed. No giant shear penalty/tiny H calculation or new EB time solver ran.
One positive stretching mechanism does not establish universal hardening of the
complete M-H/Timoshenko nonlinear frequency response.

## Retained strains and resultants

The retained cubic kinematic measures are

\[
 \Gamma_1^{[3]}=u_s+\theta w_s-\tfrac12\theta^2-\tfrac12u_s\theta^2,
 \quad
 \Gamma_2^{[3]}=w_s-\theta-u_s\theta-\tfrac12w_s\theta^2+\tfrac16\theta^3.
\]

They are evaluated on saved nonlinear coordinates; the trajectory was not
replaced by untruncated trigonometric dynamics. Sampled maxima over the entire
history are below (locations/time refer to the large case):

| Measure | .05 max | .025 max | large s/L | large t/T1 |
|---|---:|---:|---:|---:|
| u_s | 1.495435e-05 | 3.738588e-06 | 0.500000 | 0.000000 |
| w_s | 7.657787e-03 | 3.828893e-03 | 0.774940 | 0.000000 |
| theta | 7.483300e-03 | 3.741493e-03 | 0.774940 | 0.000000 |
| c | 4.602236e-06 | 1.150559e-06 | 0.185892 | 0.000000 |
| L_c_s | 5.436427e-04 | 1.359029e-04 | 0.000086 | 0.999547 |
| L_theta_s | 6.845121e-02 | 3.421748e-02 | 0.999914 | 0.498763 |
| Gamma1 | 1.534195e-05 | 3.835413e-06 | 0.185892 | 0.000000 |
| Gamma2 | 2.071786e-04 | 1.034390e-04 | 0.000086 | 0.499293 |
| surface_bending_strain | 1.711280e-03 | 8.554370e-04 | 0.999914 | 0.498763 |

These magnitudes are consistent with the declared small-motion regime: rotation
and slope below.008, surface bending strain about.00171. This is descriptive,
not a new universal validity threshold. Six stored snapshots comparing exact
reduced geometry with Gamma[3] give maximum large/half differences
Gamma1:4.04239e-10/2.52623e-11; Gamma2:1.71744e-13/5.36604e-15.
Exact reduced geometry here is only an auxiliary kinematic diagnostic.

With C=EA/(1-nu^2), S=kappa GA, Bp=EI, H=kappa GI, material-frame quantities
are N[3]=C(Gamma1[3]+nu c), Q[3]=S Gamma2[3], M=Bp theta_s, Rc=H c_s:

| Resultant | .05 sampled max | .025 sampled max |
|---|---:|---:|
| N3 | 1.534207e-07 | 3.835434e-08 |
| Q3 | 6.640340e-07 | 3.315353e-07 |
| M | 1.426067e-07 | 7.128642e-08 |
| Rc | 3.630093e-10 | 9.074716e-11 |

## Boundary forces, signs and momentum

Local t=+EX, n=-EY, k=t cross n=-EZ; theta is signed about k. Define r=(s+u,w)
in the t,n basis. Canonical force fluxes are F=(Fu,Fw)=partial V4/partial(u_s,w_s).
They equal the CUBIC RETAINED rotation of material N/Q:
Fu=trunc3(N cos theta-Q sin theta), Fw=trunc3(N sin theta+Q cos theta).
Multiplying truncated measures without truncating the product would introduce
orders absent from V4. Moment and contraction fluxes are M and Rc above.

For outward sigma=-1 left, +1 right, support reactions ON the rod are sigma F,
sigma M k and sigma Rc. Physical global forces are `(sigma Fu,-sigma Fw)`;
the physical global couple is `-sigma M` about EZ. Forces exerted by the rod on
the support have opposite signs. Rc is the scalar c-conjugate reaction, not a
Cartesian force or physical bending couple despite its force-length units.

Representative LEFT reactions on the rod (local signs) are:

| epsilon | actual t/T1 | Fu | Fw | M | Rc |
|---|---:|---:|---:|---:|---:|
| 0.050 | 0.00000000 | -1.492765e-07 | -6.627533e-07 | -1.426246e-07 | 3.668213e-10 |
| 0.050 | 0.25000000 | 5.166469e-11 | 7.793786e-10 | 1.622460e-10 | -1.301004e-13 |
| 0.050 | 0.50000000 | -1.492829e-07 | 6.640058e-07 | 1.426587e-07 | 3.668369e-10 |
| 0.050 | 1.00000000 | -1.492854e-07 | -6.636021e-07 | -1.426462e-07 | 3.668464e-10 |
| 0.025 | 0.00000000 | -3.731912e-08 | -3.313767e-07 | -7.131054e-08 | 9.170533e-11 |
| 0.025 | 0.25000000 | 1.153834e-11 | 9.755424e-11 | 2.029801e-11 | -2.926793e-14 |
| 0.025 | 0.50000000 | -3.731955e-08 | 3.315333e-07 | 7.131483e-08 | 9.170638e-11 |
| 0.025 | 1.00000000 | -3.731980e-08 | -3.314833e-07 | -7.131335e-08 | 9.170721e-11 |

Both-end forces/moments and weak-lift reactions at all nine actual snapshots
are retained in CSV/JSON. At the initial positive transverse displacement,
negative local transverse support force is restoring. Clamp reactions are
nonzero and accelerate total momentum; zero total force/moment was not imposed.

Let E_U,E_theta denote the independently evaluated STRONG interior coordinate
residuals, and ell=jp(1+c)^2 theta_t the spin density. In local signed orientation,

\[
 \dot P-[F]_0^L=\int E_U ds,\qquad
 \dot J-[r\times F+M]_0^L
   =\int(r\times E_U+E_\theta-K_4)ds.
\]

Here P=integral(m r_t ds), J=integral(r cross m r_t+ell ds). Constant translations
and the global-rotation generator do not belong to the homogeneous essential
Galerkin test space; finite-p strong residuals need not vanish under these tests.
Their contributions are explicitly integrated rather than forced to zero.
The scalar c is not cyclic, so R balance is not a standalone conserved momentum:

\[
 \int jp c_{tt}ds=[R_c]_0^L-
             \int(V_{4,c}-jp(1+c)\theta_t^2)ds.
\]

A necessary truncation qualification follows directly from the frozen V4:

\[
 K_4=V_{4,\theta}+(1+u_s)F_w-w_sF_u
 =\theta^3u_s(-C/2+2S/3)
  +(\nu C/2)c w_s\theta^2+2(C-S)u_sw_s\theta^2.
\]

It is degree4, vanishes through degree3, and is O(epsilon^5) for the prepared
scaling. The original untruncated reduced potential is rotationally objective;
its finite Taylor polynomial need not be exactly invariant under global rotation
at discarded orders. K4 is retained in the diagnostic identity, not corrected
in the RHS. In this symmetric history its integral nearly cancels; that does
not prove the generic local expression zero or exact angular conservation.
The local snapshot maxima of |K4| are7.69468e-14 and2.40441e-15, while the
maximum symmetric integral magnitudes are8.97093e-26 and2.80261e-27.
The saved auxiliary snapshot audit uses the displayed K4 expression.

Maximum discrepancies over nine snapshots:

| Diagnostic | .05 | .025 |
|---|---:|---:|
| translation_difference | 4.496410e-12 | 5.623214e-13 |
| angular_difference | 2.248205e-12 | 2.811605e-13 |
| translation_decomposition_error | 2.431479e-20 | 9.573148e-21 |
| angular_decomposition_error | 1.353720e-20 | 5.765514e-21 |
| fixed_scaled_translation_discrepancy | 1.784468e-06 | 4.463314e-07 |
| fixed_scaled_angular_discrepancy | 8.922339e-07 | 2.231656e-07 |

The difference minus its independently integrated residual/truncation expression
is near arithmetic precision. No unexplained sign or momentum identity failure
was found. This verifies accounting, not exact continuum balance of the p64
approximation. Fixed force scales are m(2omega1)^2 epsilon^2 h0 L for axial,
m omega1^2 epsilon h0 L for transverse; moment scale is transverse scale times L.
No division by an instantaneous zero or new balance acceptance threshold occurs.
Weak lift reactions are assembled separately from endpoint fluxes; agreement of
their summed balance is variational bookkeeping, not independent physical validation.

## Energy, safety and cost

Energy is the unchanged semidiscrete quartic-action total
E_h=1/2 v^T M(q)v+V_h(q), using EACH initial energy, no energy projection.

| Quantity | saved .05 tight | new .025 tight |
|---|---:|---:|
| initial_energy | 1.259738347e-09 | 3.147228160e-10 |
| relative_energy_drift_max | 2.132074057e-12 | 1.383133905e-13 |
| mass_lower_bound_min | 9.999907955e-01 | 9.999976989e-01 |
| mass_condition_bound_max | 1.000009205e+00 | 1.000002301e+00 |

Safety/energy/mass checks pass under the unchanged policy. The initial-energy
ratio is4.00269152, consistent with a leading quadratic energy plus higher orders,
not required to equal4 exactly. Relative mass bounds are the existing weighted
Gram/Loewner bounds, not new eigensolves. Min(1+c) is.9999953978 in the large case
and.9999988494 in the new case. Small energy drift does not prove spatial accuracy.
The forced leading longitudinal subsystem separately obeys
E2_t=x_t^T[f0+f2 cos(2omega1t)] and is not individually energy-conserving.
No energy fraction is used to classify a mode.

Exactly ONE new ODE reachedT1 in118.547837s, with5491 accepted steps,
nfev38439, njev1, nlu10982 and38435 state-dependent mass factorizations.
The old large run remains immutable (48.004957s, nfev38445,njev2,nlu2566).
The greater LU count makes a linear cost forecast unreliable; rejected-step
counts are not inferred. Charged primary numerical work is237.72s, including
19.94s supplemental postprocessing and a conservative35s allowance for a failed
report-only list conversion. This error was corrected before final reporting;
it never affected the integrator, RHS, coefficients or historical sources.
The initial manifest path-separator preflight also stopped before any ODE.
All numerical work, tests/cache verification remain under600s.

## Reproduction, artifacts and verification

New result: `results/nlsp_planar_physical_sanity_checks/888042e17cfc315a/`.
It contains manifests/source hashes, the new q/v history, actual times,
initial coefficients, all8 comparison curves/tables, physical profiles at
six timestamps, resultants/strains, nine reaction/balance snapshots, energy,
mass/safety and original execution provenance. The one-run execution snapshot
is preserved separately; final code additions only strengthened provenance and
saved-data reporting, and the half-history SHA is unchanged. Old bundles were
not overwritten.

The new focused CLI has a genuinely different diagnostic I/O contract from
initial preparation and one-period convergence. It reuses `PlanarGalerkin`,
the existing Radau runner and frozen prepared evaluator; there is no second
physics solver. Config is
[nlsp_planar_physical_sanity_checks.json](../../data/input/nlsp_planar_physical_sanity_checks.json).

```powershell
python scripts/analysis/check_weakly_nonlinear_planar_physics.py --compute
python scripts/analysis/check_weakly_nonlinear_planar_physics.py --report-only results/nlsp_planar_physical_sanity_checks/888042e17cfc315a
python scripts/analysis/check_weakly_nonlinear_planar_physics.py --plot-only results/nlsp_planar_physical_sanity_checks/888042e17cfc315a
```

Matching compute/report/plot read saved evidence with ZERO new ODE/BVP/eigensolves
or symbolic model derivations. Cache identity includes source manifests and action
hash, amplitudes, common projection policy, p/quadrature, exact time prescription,
horizon/sampling contract, code, dependency versions and BLAS thread environment.
A missing historical source is an explicit blocker, not permission to recompute it.
Tests prohibit hidden integration/eigensolver calls.

Three figures, PDF+PNG, are saved in `figures/`:
`amplitude_normalized_motion`, `leading_axial_response`,
`support_reaction_diagnostics`. Sorted/modal/energy tracking is not performed.
Full numerical tables are CSV and JSON; no figure conceals a strict FAIL.
37 independent exact Fraction checks verify limited identities/dimensions;
62 distinct targeted tests ultimately PASS (60 new tests and2 unchanged
one-T1 regressions). Two exact-equality assertions initially failed on arithmetic
representation alone; the corrected assertions use the original amplitude
formula and a4-epsilon arithmetic tolerance and both pass on targeted recheck.
No numerical/physical acceptance criterion was changed. Tests cover source
reuse, IC powers/parity, derivatives,
canonical flux signs, integrated balance accounting, safety/cache and frozen hashes.
Historical strict XFAIL markers and thresholds are unchanged.

## Status and stopping boundary

| Status | Outcome | Scope |
|---|---|---|
| NLSP_SANITY_SECOND_ORDER_COMPARISON | PARTIAL | measured leading consistency; finite-p/strict and half-case uncertainty |
| NLSP_SANITY_AMPLITUDE_SCALING | PARTIAL | expected pattern observed; one half run, no independent p/time certificate |
| NLSP_SANITY_REFLECTION_SYMMETRY | PASS | exact local parity and numerical RHS spots, no negative ODE |
| NLSP_SANITY_CLASSICAL_STRETCHING_LIMIT | PASS | formal bulk benchmark, not uniform finite c-clamp elimination |
| NLSP_SANITY_STRAINS | PASS | retained diagnostic complete and existing safety passes, not applicability proof |
| NLSP_SANITY_REACTIONS_AND_MOMENTUM | PARTIAL | residual/truncation accounting closes; no exact continuum certification |
| NLSP_PLANAR_PHYSICAL_SANITY_CHECKS | DIAGNOSTIC_COMPLETE_WITH_QUALIFICATIONS | all permitted checks completed, physical V0 accuracy remains unestablished |

No unexplained action/sign/assembly defect was found. The numerical motion is
consistent with the expected weakly nonlinear patterns in this bounded test,
which supports using the method for the declared task with visible qualifications.
It does not establish continuum convergence, experiment/3D truth, a periodic
orbit, universal amplitudes/geometries or out-of-plane stability.

The existing LONG line remains CLOSED in its adopted scope; EB/RLB spring/KV
remains PAUSED; the angular same-clamp out-of-plane reference is UNAVAILABLE.
Historical zero-u/c PARTIAL and prepared strict PARTIAL are not promoted.
[NLSP-D08](../memory/decisions.md#nlsp-d08) and
[NLSP-K08](../memory/knowledge.md#nlsp-k08) record this bounded result and stop.
No new physical coefficient, angular/joint model, FEM, Floquet, critical-amplitude
search or nonlinear periodic-orbit calculation was introduced. No next stage is
selected or automatically authorized by this result.
