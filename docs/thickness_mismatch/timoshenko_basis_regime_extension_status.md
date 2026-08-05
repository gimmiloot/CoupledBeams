# Timoshenko basis-regime extension status

## Scope and result

The isotropic variable-length Timoshenko evaluator now represents every
spatial-root regime admitted by the existing physical coefficients at positive
frequency. No Euler--Bernoulli or Timoshenko field equation, joint condition,
unknown ordering, shear coefficient, geometry parameter, `K=10` definition,
`delta_f`, 10% threshold, or root-quality tolerance was changed.

The finite solver-readiness validation is complete:

- `basis_formula_gate = PASS`;
- `old_regime_regression_gate = PASS`;
- `seven_case_reference_gate = PASS`;
- `optimization_recheck_gate = PASS`;
- `optimization_equivalence_gate = PASS` on 24/24 cases;
- `full_grid_solver_readiness_gate = PASS` with 24/24 resolved full references.

The versioned evidence is under
`results/article_epsilon_upper_envelope/solver_readiness_v2/`. This PASS is a
finite numerical readiness gate, not a continuous-domain completeness proof or
authorization to run the 1554-point grid.

## Frozen characteristic problem

For `z=r^2`, the unchanged section characteristic polynomial is

```text
B z^2 + Omega^2 (B m/K + J) z + (J m Omega^4/K - m Omega^2) = 0,
K = kappa G A,  B = E I,  m = rho A,  J = rho I.
```

The state relations used to derive and independently check the basis are

```text
M = B psi',              Q = K (w' - psi),
Q' = -m Omega^2 w,       M' = -Q - J Omega^2 psi.
```

Writing `A0=B m/K`, the discriminant is

```text
D = Omega^4 (A0-J)^2 + 4 B m Omega^2 > 0.
```

The two `z` roots are therefore real. Their sum is negative, so two positive
roots and a two-hyperbolic-pair regime are impossible for the current positive
physical coefficients. The supported regimes are exactly:

- `mixed_hyperbolic_trigonometric`: `z_a=a^2>0`, `z_b=-b^2<0`;
- `cutoff_limit`: `z_a=0`, `z_b=-b^2<0`;
- `two_trigonometric_pairs`: `z_a=-a^2<0`, `z_b=-b^2<0`.

The boundary is `Omega_c=sqrt(K/J)`, where the constant polynomial term
vanishes. It is a basis-regime boundary, not a physical `delta_f` failure.

## Clamped-reduced local columns

Let `alpha=m Omega^2/K`. In the mixed regime define
`h=a+alpha/a`, `q=b-alpha/b`; above cutoff define
`h=a-alpha/a`, `q=b-alpha/b`. The evaluator uses two real columns that already
satisfy `w(0)=psi(0)=0`.

For the two-trigonometric regime, before a nonzero invertible column scaling,

```text
d0 = b^2-a^2,                 d1 = a-(h/q)b,
w0 = (cos(ax)-cos(bx))/d0,    psi0 = (-h sin(ax)+q sin(bx))/d0,
w1 = (sin(ax)-(h/q)sin(bx))/d1,
psi1 = h(cos(ax)-cos(bx))/d1.
```

`w'` and `psi'` are differentiated analytically; the physical endpoint rows
remain `Q=K(w'-psi)` and `M=B psi'`. The coefficient order and the global 6x6
joint matrix rows/signs are unchanged. Each arm classifies its own section, so
one arm may be mixed while the other is above cutoff.

At cutoff the finite canonical limit is

```text
w0 = (1-cos bx)/b^2,          psi0 = (alpha x+q sin bx)/b^2,
w1 = sin(bx)/b,               psi1 = -q(1-cos bx)/b.
```

Half-angle difference formulas avoid cancellation. Explicit finite nonzero
column factors connect these canonical columns to the legacy mixed columns and
keep the local clamp derivative matrix full rank. This scaling changes neither
the determinant zeros nor nullity.

## Independent checks and numerical stabilization

An audit-only first-order state-space oracle uses `scipy.linalg.expm` on
`[w, psi, Q, M]`; it is not a production or grid solver. Below, at, and above
cutoff the maximum scaled closed-form/oracle difference over the seven former
blockers is recorded in `oracle_comparison.csv`. Cutoff continuity and rank are
recorded in `cutoff_continuity.csv`.

The unrelated EB primary/strict failures were traced to cancellation and
column imbalance in the legacy hyperbolic representation. Production root
search now evaluates an exactly column-equivalent EB matrix with positive
transform determinant `exp(-(x1+x2))`. SVD diagnostics use continuous floored
row equilibration; independent scan roots are locally cross-validated with the
unchanged acceptance tolerances. Detection records are deduplicated before
close-seed counting, and a same-sign monotone boundary tail is no longer
reported as an unresolved interior root.

Four disputed EB roots were checked locally at 80 decimal digits. Their maximum
double/high-precision difference is below `3e-12`; high precision is not used
by production searches.

## Validation and limitations

S3 regressions remain `N_true=4`:

```text
S3_12 delta_f,5 = 0.11739469908796035
S3_14 delta_f,5 = 0.10050934855181458
```

Both are within the pre-existing numerical tolerance of the immutable target
values. All seven former blockers now have resolved full K10/root-11 references
and match paired+auto through the required guard. V3 is the one validation case
whose full first 11 Timoshenko roots actually enter
`two_trigonometric_pairs`; paired+auto correctly stops before those unneeded
upper roots after confirming failure at mode 2.

Full K10 remains the default. Paired+auto remains an optional computational
optimization. The coarse grid, Rule A/B/S, mode energy/shapes, FEM/3D FEM, and
anisotropic models were not run or modified in this work.
