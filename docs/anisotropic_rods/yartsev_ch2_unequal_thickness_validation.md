# Chapter-2 unequal-thickness validation

Overall unequal-thickness validation: **IN_PROGRESS**

UT-0 section-routing and beta=0 stepped-reference smoke: **PASS**

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
