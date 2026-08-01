# Chapter-2 unequal-thickness validation

Overall unequal-thickness validation: **IN_PROGRESS**

UT-0 section-routing and beta=0 stepped-reference smoke:
**PASS**

UT-1 full beta=0 unequal-thickness validation:
**PARTIAL_PASS**

UT-1a beta=0 EB FEM matrix-level exchange audit:
**INCONCLUSIVE_NUMERICAL_AUDIT**

## 1. UT-0 scope

UT-0 is a finite validation gate for the existing Chapter-2 Timoshenko and
rectangular Euler--Bernoulli models. It is not a new production model and does
not complete the unequal-thickness validation program.

The only calculated case uses the real-elastic HMS/DX-209 material,
`theta_1=theta_2=0 deg`, `L_1=L_2=0.2 m`, `b_1=b_2=0.020 m`, the existing
shear factor `5/6`, and `beta=0 deg`. The rectangular dimension `a_i` is
parallel to `e_z`. Geometry is prescribed directly; total mass is not held
fixed.

Three independent `RodPoint` objects are constructed with separate `Geometry`
objects:

| point | `a` (m) | role |
|---|---:|---|
| `a4` | 0.004 | spectral arm 1 |
| `a5` | 0.005 | existing-baseline regression |
| `a6` | 0.006 | spectral arm 2 |

The spectral smoke is only `(point_1,point_2)=(a4,a6)`.

## 2. Section-property routing

Each point independently uses

\[
A=ab,\qquad I_y=\frac{a^3b}{12},\qquad
I_p=\frac{ab(a^2+b^2)}{12}.
\]

The transverse-shear stiffness remains the existing quantity

\[
K_s=\kappa G_{xz}A,
\]

and the generalized rectangular torsion series is unchanged. At `theta=0`,
`Sbar16=0` and `C_T=Cbar=C_SV`.

| point | `A` (m^2) | `I_y` (m^4) | `I_p` (m^4) | `K_s` (N) | `Cbar=C_T` (N m^2) | series terms | estimated relative tail |
|---|---:|---:|---:|---:|---:|---:|---:|
| `a4` | `8.000000000000e-5` | `1.066666666667e-10` | `2.773333333333e-9` | `2.000000000000e5` | `1.118656336082` | 298 | `9.927360166223e-13` |
| `a5` | `1.000000000000e-4` | `2.083333333333e-10` | `3.541666666667e-9` | `2.500000000000e5` | `2.106097187308` | 298 | `9.927492364894e-13` |
| `a6` | `1.200000000000e-4` | `3.600000000000e-10` | `4.360000000000e-9` | `3.000000000000e5` | `3.503243508855` | 298 | `9.928476259971e-13` |

For all three points, `Ex=1.91e11 Pa` and `Gxz=3.0e9 Pa`; both maximum
cross-point relative residuals are zero. `Cbar` is recomputed from each
geometry and differs across the three sections. The explicit `a5` construction
matches a point made through `cantilever_geometry(0.2)` with zero maximum
relative residual over the compared geometry, material, and torsion fields.

The maximum section-identity, scaled-`Sbar16`, `C_T=Cbar`, and direct-`Cbar`
recomputation residuals are all zero. The largest estimated relative
torsion-series tail is `9.928476259971e-13`, below the unchanged `1e-12`
tolerance. All physically positive section quantities are finite and positive.

## 3. Arm-specific state coefficients

The diagnostic frequency is fixed before evaluation at `f_diag=500 Hz`.
For the accepted state order `[w,psi,Phi,Q,M,M_T]`, each Timoshenko point is
checked against

\[
A_{w,\psi}=1,\quad A_{w,Q}=\frac{1}{K_s},\quad
A_{\psi,M}=\frac{1}{E_xI_y},\quad A_{\Phi,M_T}=\frac{1}{C_T},
\]

\[
A_{Q,w}=-\rho A\omega^2,\quad
A_{M,\psi}=-\rho I_y\omega^2,\quad A_{M,Q}=-1,\quad
A_{M_T,\Phi}=-\rho I_p\omega^2.
\]

The two material bending--torsion entries vanish at `theta=0`. The EB check
uses the unchanged equations

\[
w'=\psi,\quad \psi'=\frac{M}{E_xI_y},\quad
\Phi'=\frac{M_T}{C_{SV}},\quad Q'=-\rho A\omega^2w,
\quad M'=-Q,\quad M_T'=-\rho I_p\omega^2\Phi.
\]

The maximum Timoshenko and EB coefficient residuals are both zero. Because
`A`, `I_y`, `I_p`, `K_s`, and torsional stiffness come from each point, the
checked coefficients change with the selected arm geometry.

## 4. Existing state scaling

For each point, the physical transfer recovered by the existing scaled
`state_corrected` implementation is compared with the direct physical matrix

\[
\exp\!\left(A(\omega_{diag})L\right).
\]

The relative Frobenius residuals are `2.299066698625e-17` for `a4`,
`1.634715368057e-15` for `a5`, and `5.875495157454e-16` for `a6`. The maximum
is `1.634715368057e-15`, below the fixed `1e-10` gate. No private scaling
helper is exposed or copied.

## 5. Independent beta=0 stepped reference

For the accepted book state, the local joint-state map is introduced only in
the diagnostic script as

\[
y_{2,J}=R_0y_{1,J},\qquad
R_0=\operatorname{diag}(1,-1,-1,-1,1,1).
\]

Thus `w_2=w_1`, `psi_2=-psi_1`, `Phi_2=-Phi_1`, `Q_2=-Q_1`,
`M_2=M_1`, and `M_T,2=M_T,1`. Deterministic states confirm that this map
satisfies the independently evaluated `J_book(0)` conditions.

The direct stepped Timoshenko boundary matrix is

\[
D_{step}^{T}=C_2T_2^{-1}R_0T_1B_1,
\]

implemented as `C2 @ solve(T2, R0 @ T1 @ B1)`. `B1` is the existing physical
`book_slope_clamp` map and `C2` is the existing straight right clamp selector
built from arm 2, including the `a6` shear stiffness. The EB reference is

\[
D_{step}^{EB}=C_2^{EB}(T_2^{EB})^{-1}R_0T_1^{EB}B_1^{EB},
\]

with the corresponding existing EB transfer and clamp maps. Both stepped
matrices are `3 x 3`; neither calls `joint_matrix_book`, a coupled boundary
builder, coupled block assembly, or a copied coupled determinant.

## 6. First-three sorted spectra

The unchanged `find_elastic_roots` is used with `state_corrected`, three
positive roots, `10 Hz` scan step, `5000 Hz` initial maximum, and `100000 Hz`
hard maximum. These are sorted spectra only; no branch identity is asserted.

| model | mode | coupled (Hz) | direct stepped (Hz) | relative difference |
|---|---:|---:|---:|---:|
| Timoshenko | 1 | `332.915895232818` | `332.915895232821` | `1.109836531947e-14` |
| Timoshenko | 2 | `787.843604846504` | `787.843604846582` | `9.956788066229e-14` |
| Timoshenko | 3 | `951.334941111405` | `951.334941111270` | `1.420883899135e-13` |
| EB | 1 | `334.385945984019` | `334.385945984032` | `4.079842675938e-14` |
| EB | 2 | `787.843604846502` | `787.843604846559` | `7.258354199005e-14` |
| EB | 3 | `974.239141599259` | `974.239141599289` | `3.057355243823e-14` |

The Timoshenko maximum coupled--stepped relative difference is
`1.420883899135e-13`; the EB maximum is `7.258354199005e-14`. Both pass the
fixed `1e-8` gate. The Timoshenko--EB difference is not an acceptance
criterion and no continuum-model accuracy ranking is inferred.

## 7. Root quality and thresholds

Every root is finite, positive, accepted, and retains scaled plus
physical-raw determinant/SVD diagnostics. The unchanged rectangular-gate rule
passes a root when either scaled quality or normalized physical-raw quality
passes, provided the root status is not rejected.

Across all 12 root/formulation evaluations, the maxima are:

| residual | maximum | threshold |
|---|---:|---:|
| scaled determinant | `2.779978222317e-4` | `1e-8` |
| scaled relative singular | `5.556871988090e-4` | `1e-8` |
| physical-raw normalized determinant | `4.401217951480e-18` | `1e-8` |
| physical-raw relative singular | `4.096193230446e-14` | `1e-8` |

Two stepped torsion-root rows do not pass through their scaled representation,
but both pass the pre-existing alternative physical-raw normalized quality
gate. No threshold is waived or changed.

All fixed UT-0 thresholds are:

- section identities, scaled zero quantities, baseline, common material
  values, and state coefficients: `1e-12`;
- estimated rectangular torsion-series tail: `1e-12`;
- physical-transfer relative Frobenius residual: `1e-10`;
- determinant and relative singular root residuals: `1e-8`;
- coupled--stepped sorted-spectrum relative difference: `1e-8`.

## 8. Status and exclusions

All UT-0 scientific gates and the targeted/regression tests pass. Therefore:

```text
Overall unequal-thickness validation: IN_PROGRESS
UT-0 section-routing and beta=0 stepped-reference smoke: PASS
```

UT-0 explicitly excludes `beta=30 deg`, `beta=90 deg`, the swapped `(6,4)`
spectrum, exchange symmetry, root 7, EB/Timoshenko/3D FEM, mesh refinement,
`theta != 0`, complex roots, damping, mode shapes, MAC, mode tracking,
parameter maps, mass-preserving parametrization, plots, and any new production
anisotropic API. These calculations were not run.

Generated text-only evidence is written under
`results/anisotropic_rods/yartsev_ch2_unequal_thickness_validation/ut0_smoke/`.

## 9. UT-1 — full beta=0 unequal-thickness validation

UT-1 extends only the `beta=0` gate. It retains elastic HMS/DX-209,
`theta_1=theta_2=0`, `L_1=L_2=0.2 m`, `b_1=b_2=0.020 m`, the existing
`5/6` shear factor, and `a_i || e_z`. It evaluates the three prescribed
section orders `(5,5)`, `(4,6)`, and `(6,4) mm`. No mass-preserving
parameterization is used.

The existing UT-0 direct stepped transfer builders remain script-local. The
existing root finder is used without changes with seven positive roots,
`10 Hz` scan steps, `5000 Hz` initial maximum, and `100000 Hz` hard maximum.
Every spectrum completed inside the initial `5000 Hz` scan range.

## 10. Seven-root continuum spectra

All 12 main spectra contain exactly seven finite, positive, accepted roots.
The following values are sorted frequencies; no branch identity is claimed.

| case | model | formulation | f1 | f2 | f3 | f4 | f5 | f6 | f7 | units |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| `(5,5)` | Timoshenko | coupled | 351.012188841 | 766.861317902 | 951.727469927 | 1533.722635804 | 1820.109925368 | 2300.583953706 | 2913.573068905 | Hz |
| `(5,5)` | Timoshenko | stepped | 351.012188841 | 766.861317902 | 951.727469927 | 1533.722635804 | 1820.109925368 | 2300.583953704 | 2913.573068898 | Hz |
| `(5,5)` | EB | coupled | 353.181083065 | 766.861317902 | 973.557255447 | 1533.722635804 | 1908.561482197 | 2300.583953706 | 3067.445271607 | Hz |
| `(5,5)` | EB | stepped | 353.181083065 | 766.861317903 | 973.557255447 | 1533.722635804 | 1908.561482197 | 2300.583953705 | 3067.445271606 | Hz |
| `(4,6)` | Timoshenko | coupled | 332.915895233 | 787.843604847 | 951.334941111 | 1401.537220814 | 1734.935407839 | 2308.867763347 | 2867.089958267 | Hz |
| `(4,6)` | Timoshenko | stepped | 332.915895233 | 787.843604847 | 951.334941111 | 1401.537220814 | 1734.935407839 | 2308.867763348 | 2867.089958260 | Hz |
| `(4,6)` | EB | coupled | 334.385945984 | 787.843604847 | 974.239141599 | 1401.537220814 | 1809.159163981 | 2308.867763347 | 2878.317774375 | Hz |
| `(4,6)` | EB | stepped | 334.385945984 | 787.843604847 | 974.239141599 | 1401.537220815 | 1809.159163983 | 2308.867763348 | 2878.317774374 | Hz |
| `(6,4)` | Timoshenko | coupled | 332.915895233 | 787.843604847 | 951.334941111 | 1401.537220814 | 1734.935407839 | 2308.867763347 | 2867.089958267 | Hz |
| `(6,4)` | Timoshenko | stepped | 332.915895233 | 787.843604847 | 951.334941111 | 1401.537220815 | 1734.935407839 | 2308.867763347 | 2867.089958265 | Hz |
| `(6,4)` | EB | coupled | 334.385945984 | 787.843604847 | 974.239141599 | 1401.537220814 | 1809.159163981 | 2308.867763347 | 2878.317774375 | Hz |
| `(6,4)` | EB | stepped | 334.385945984 | 787.843604847 | 974.239141599 | 1401.537220815 | 1809.159163982 | 2308.867763348 | 2878.317774376 | Hz |

The maximum coupled--stepped relative differences are:

| case | Timoshenko | EB | threshold |
|---|---:|---:|---:|
| `(5,5)` | `2.633987928940e-12` | `4.490478570501e-13` | `1e-8` |
| `(4,6)` | `2.278899522742e-12` | `7.356003544981e-13` | `1e-8` |
| `(6,4)` | `1.073601495836e-12` | `3.580675762534e-13` | `1e-8` |

No near-degenerate cluster below the fixed `1e-3` relative-gap diagnostic
required reassignment in these comparisons.

## 11. Root quality

The two existing acceptance branches are retained. Across the 84 roots in the
12 main spectra plus the 14 homogeneous-reference roots, 63 pass through the
scaled branch and 35 through the normalized physical-raw branch. No root is
rejected.

| accepted quality basis | maximum determinant residual | maximum relative singular residual | threshold |
|---|---:|---:|---:|
| scaled | `1.868018030350e-15` | `4.890652111179e-14` | `1e-8` |
| physical raw | `5.731647065818e-18` | `6.054486459828e-13` | `1e-8` |

For completeness, maxima across representations, irrespective of which
branch accepted a root, are `3.748369515937e-4` for the scaled determinant,
`1.843259330288e-2` for the scaled relative singular residual,
`5.731647065818e-18` for the normalized raw determinant, and
`6.054486459828e-13` for the raw relative singular residual. Large scaled
values do not waive or replace the predeclared alternative raw gate.

## 12. Equal-thickness and exchange regressions

For `(5,5)`, the coupled and direct stepped spectra are additionally compared
with independent homogeneous fixed--fixed rods of length `0.4 m` using the
existing straight Timoshenko and EB builders. The maximum pairwise relative
differences are `4.716236859975e-12` for Timoshenko and
`7.418406988057e-13` for EB, both below `1e-8`.

The `beta=0` stepped-order exchange/reflection maxima for `(4,6)` versus
`(6,4)` are:

| model | formulation | maximum relative difference | threshold |
|---|---|---:|---:|
| Timoshenko | coupled | `2.014339082597e-14` | `1e-8` |
| Timoshenko | direct stepped | `1.818139283769e-12` | `1e-8` |
| EB | coupled | `3.141979986753e-15` | `1e-8` |
| EB | direct stepped | `6.208552453820e-13` | `1e-8` |

This is only a straight stepped-order result at `beta=0`; it is not
generalized to angular joints.

## 13. Independent EB FEM

UT-1 reuses the existing independent Hermite EB bending and linear
Saint-Venant torsion elements, their consistent translational and torsional
mass matrices, common-joint congruence map, fixed outer nodes, and dense
symmetric generalized eigensolver. The analytic determinant is not used to
obtain FEM eigenvalues, and `J_book` is not imposed as a FEM constraint.

All three `N_1=N_2=64` assemblies have 381 reduced DOFs, zero stiffness and
mass symmetry residuals at saved precision, positive reduced mass, zero
spurious modes, seven finite positive roots, zero joint-kinematic residual,
and maximum joint-equilibrium residual below the existing `1e-7` gate.

| case | E3(64) | E6(64) | root-7 error (diagnostic) | E3 threshold | E6 threshold |
|---|---:|---:|---:|---:|---:|
| `(5,5)` | `2.510018347396e-5` | `2.259131248401e-4` | `4.016439332810e-4` | `1e-4` | `5e-4` |
| `(4,6)` | `2.920357232077e-5` | `2.773278662256e-4` | `3.791829271970e-4` | `1e-4` | `5e-4` |
| `(6,4)` | `2.918306050920e-5` | `2.773256322238e-4` | `3.791829001104e-4` | `1e-4` | `5e-4` |

All single-mesh structure and accuracy gates pass. Root 7 is present but has
no accuracy threshold.

## 14. Representative FEM convergence

Only `(4,6)` is evaluated on the fixed `(16,16)`, `(32,32)`, `(64,64)`
sequence. Acceptance uses raw frequencies and aggregate `E6`; no extrapolated
frequency or per-mode monotonicity rule is introduced.

| elements per arm | E3 | E6 | aggregate order from previous |
|---:|---:|---:|---:|
| 16 | `4.671130953825e-4` | `4.438659031444e-3` | n/a |
| 32 | `1.167706851710e-4` | `1.109380926989e-3` | `2.000369060806` |
| 64 | `2.920357232077e-5` | `2.773278662256e-4` | `2.000090334837` |

Both `E6` reductions are strict and both empirical orders exceed the fixed
`1.5` minimum. The representative convergence gate passes.

## 15. FEM exchange result and UT-1 status

The independent mesh-64 FEM `(4,6)`--`(6,4)` maximum relative spectral
difference is `5.762284502891e-8`, above the fixed `1e-8` exchange threshold.
Both FEM assemblies independently pass their structural and analytic-accuracy
gates; their matrix symmetry and joint-kinematic residual differences are
zero at saved precision. This isolated FEM exchange failure is not treated as
evidence of a continuum-model error, and no element, eigensolver, matching
rule, or threshold is changed.

The status rules therefore give:

```text
Overall unequal-thickness validation: IN_PROGRESS

UT-0 section-routing and beta=0 stepped-reference smoke:
PASS

UT-1 full beta=0 unequal-thickness validation:
PARTIAL_PASS
```

## 16. Timoshenko--EB diagnostic

The maximum first-six relative differences are `4.859682132184e-2` for
`(5,5)`, `4.278185562840e-2` for `(4,6)`, and `4.278185562839e-2` for
`(6,4)`. These values have no PASS threshold and do not rank the accuracy of
the continuum models or establish branch identity.

## 17. UT-1 exclusions and evidence

UT-1 did not add or run any nonzero `beta`, `theta != 0`, complex roots,
damping, Timoshenko FEM, 3D FEM, mode-shape study, MAC, branch tracking,
parameter map, thickness grid, mass-preserving parameterization, plot, PDF,
production anisotropic API, new root finder, or slender-limit study. No
`J_book`, continuum equation, scaling formula, FEM element, joint map,
historical threshold, or release version is changed.

Generated text-only evidence is under
`results/anisotropic_rods/yartsev_ch2_unequal_thickness_validation/ut1_beta0/`.

## 18. UT-1a — beta=0 EB FEM matrix-level exchange audit

UT-1a isolates the only failed UT-1 criterion: the mesh-64 EB FEM exchange
difference `5.762284502891e-8 > 1e-8` between `(a_1,a_2)=(4,6)` and
`(6,4) mm`. It asks whether the two FEM assemblies/reductions cease to be
exchange-equivalent, or whether the difference first appears when two
congruent generalized eigenproblems are solved in different DOF orderings.
The audit does not change the UT-1 status.

Only elastic HMS/DX-209, `theta_1=theta_2=0`, `beta=0`,
`L_1=L_2=0.2 m`, `b_1=b_2=0.020 m`, `N_1=N_2=64`, and seven native FEM
roots are used. Exactly two production FEM assemblies are constructed. No
continuum roots or new mesh sequence are evaluated.

## 19. Full and reduced DOF orderings

The production assembly metadata and reduction entries were checked before
constructing either permutation. Each full arm block contains all `N+1`
local node triples in the order

```text
outer clamp -> internal nodes -> joint endpoint,
local node DOFs = [w, psi, Phi].
```

Thus the full vector and its a-priori arm-swap permutation are

```text
x = [x_1, x_2]^T,
P_f = [[0, I_{3(N+1)}], [I_{3(N+1)}, 0]],
x_64 = P_f x_46.
```

No signs are applied to the full local nodal DOFs. For `N=64`, each arm block
has 195 entries and the full system size is `390 x 390`.

The reduced ordering was confirmed as

```text
q = [q_1,int, q_2,int, w_J, theta_t, theta_n]^T,
q_i,int = [w_i,1, psi_i,1, Phi_i,1, ..., w_i,N-1, psi_i,N-1, Phi_i,N-1]^T.
```

Each internal block has 189 entries and the reduced system size is
`381 x 381`. The signed permutation is

```text
P_r = [[0, I_{3(N-1)}, 0],
       [I_{3(N-1)}, 0, 0],
       [0, 0, diag(1,-1,-1)]],
q_64 = P_r q_46.
```

The signs are fixed by the existing endpoint maps, not fitted from assembled
matrices. At `beta=0` they give

```text
[w,psi,Phi]_1 = [w,-theta_n, theta_t],
[w,psi,Phi]_2 = [w, theta_n,-theta_t].
```

After exchanging arms while preserving each arm's local DOFs, `w_J` is
unchanged while `theta_t` and `theta_n` both change sign.

## 20. Matrix-level identities

The hard relative Frobenius and relative maximum-entry threshold is `1e-13`.
All saved absolute and relative residuals below are exactly zero at binary64
saved precision:

| identity | relative Frobenius | relative max entry | threshold | status |
|---|---:|---:|---:|---|
| `P_f^T P_f=I` | `0` | `0` | `1e-13` | PASS |
| `P_f^2=I` | `0` | `0` | `1e-13` | PASS |
| `P_r^T P_r=I` | `0` | `0` | `1e-13` | PASS |
| `P_r^2=I` | `0` | `0` | `1e-13` | PASS |
| `R_64 P_r=P_f R_46` | `0` | `0` | `1e-13` | PASS |
| `K_full,64=P_f K_full,46 P_f^T` | `0` | `0` | `1e-13` | PASS |
| `M_full,64=P_f M_full,46 P_f^T` | `0` | `0` | `1e-13` | PASS |
| reduced `K` forward congruence | `0` | `0` | `1e-13` | PASS |
| reduced `M` forward congruence | `0` | `0` | `1e-13` | PASS |
| reduced `K` reverse congruence | `0` | `0` | `1e-13` | PASS |
| reduced `M` reverse congruence | `0` | `0` | `1e-13` | PASS |

No stiffness, mass, joint-map, or reduction asymmetry was detected within the
declared matrix thresholds. Entrywise first/maximum mismatch fields and DOF
labels are retained in the evidence CSV; they are empty because every matrix
residual is zero.

## 21. Native eigenproblems and backward quality

The two native systems were solved only through the existing
`solve_two_arm_eb_fem(..., num_roots=7)`. The frequencies are:

| mode | `(4,6)` Hz | `(6,4)` Hz | native relative exchange difference |
|---:|---:|---:|---:|
| 1 | 334.385921722132 | 334.385940990400 | `5.762284502891e-8` |
| 2 | 787.866612694194 | 787.866596534094 | `2.051121256078e-8` |
| 3 | 974.239157989062 | 974.239147690443 | `1.057093545693e-8` |
| 4 | 1401.677810883528 | 1401.677811343092 | `3.278670970309e-10` |
| 5 | 1809.159236843100 | 1809.159235557093 | `7.108310922016e-10` |
| 6 | 2309.508076717324 | 2309.508071559309 | `2.233382349116e-9` |
| 7 | 2879.409183334115 | 2879.409183256152 | `2.707628539793e-11` |

The maximum is exactly the saved UT-1 value
`5.762284502890923e-8`; the historical `1e-8` criterion remains failed and
UT-1 remains `PARTIAL_PASS`.

For the required normalized backward residual

```text
||K v_i-lambda_i M v_i||_2 /
max(||K v_i||_2+|lambda_i| ||M v_i||_2,tiny),
```

the native maximum is `7.283856211250112e-7`, above its fixed `1e-8` audit
threshold. The mass-orthonormality residuals are
`1.900779854094985e-15` for `(4,6)` and `1.633656423042813e-15` for
`(6,4)`, both below `1e-10`.

## 22. Transport, Rayleigh, and canonicalized solves

The transported vectors `P_r v_i^46` and `P_r^T v_i^64` were evaluated in
the opposite matrices without manual renormalization. Their maximum backward
residual is `7.283856211249712e-7`, above `1e-8`. The maximum transported
Rayleigh relative difference is `1.479695290188481e-7`, above `1e-10`.
These failures mirror the raw native backward quality and prevent the
`PASS_MATRIX_EQUIVARIANCE` classification even though matrix congruence is
exact at saved precision.

Both canonicalized pairs were solved with the same direct call used by the
production solver:

```python
eigh(K, M, subset_by_index=[0, 6], check_finite=True)
```

The maximum canonicalized spectrum relative difference is `0`, below
`1e-10`. These canonicalized and transported values are diagnostic only; they
do not replace the native FEM frequencies used by UT-1.

## 23. Conditioning and gaps

| ordering | matrix | norm_2 | cond_2 |
|---|---|---:|---:|
| `(4,6)` | `K` | `1.080861893172272e11` | `2.672016397778642e11` |
| `(4,6)` | `M` | `5.923195460186998e-4` | `3.219376675161722e7` |
| `(6,4)` | `K` | `1.080861893172272e11` | `2.672013117214675e11` |
| `(6,4)` | `M` | `5.923195460186998e-4` | `3.219376678387385e7` |

The minimum relative nearest-neighbor generalized-eigenvalue gap is
`0.346005296755862`; the corresponding minimum relative frequency gap is
`0.191300610087940`. No acceptance threshold or rigorous forward-error bound
is assigned to condition numbers or gaps.

## 24. UT-1a classification

The matrix-level gates and canonicalized solve pass, but the prescribed native
and transported backward-residual and transported-Rayleigh gates do not.
Therefore the status is, without threshold adjustment:

```text
Overall unequal-thickness validation: IN_PROGRESS

UT-0 section-routing and beta=0 stepped-reference smoke:
PASS

UT-1 full beta=0 unequal-thickness validation:
PARTIAL_PASS

UT-1a beta=0 EB FEM matrix-level exchange audit:
INCONCLUSIVE_NUMERICAL_AUDIT
```

Interpretation classification: `INCONCLUSIVE`. Assembly or reduction
asymmetry was not detected within the declared matrix thresholds, but the
failed raw numerical-quality gates prevent the stronger conclusion that the
remaining spectral difference is established as a DOF-order-sensitive
eigensolver error. The large condition numbers are consistent diagnostic
context, not a proven forward-error bound.

## 25. UT-1a exclusions and evidence

UT-1a did not change FEM element matrices, joint maps, reduction, boundary
conditions, eigensolver settings, root finder, continuum equations, old
thresholds, or statuses. It did not use balancing, extrapolated frequencies,
Richardson corrections, new mesh sequences, nonzero `beta`, Timoshenko FEM,
3D FEM, MAC, physical mode-shape interpretation, branch tracking, parameter
maps, plots, or PDF.

The four text-only evidence files are under
`results/anisotropic_rods/yartsev_ch2_unequal_thickness_validation/ut1a_fem_exchange_audit/`.
