# Design note: Timoshenko basis across the first cutoff

## Status and scope

This is the historical design note that preceded the implementation. The
extension has now been implemented and independently validated without changing
the Euler--Bernoulli or Timoshenko equations, joint rows, unknown ordering,
shear coefficient, parameter definitions, tolerances, or production default.
The canonical outcome is the
[implementation status](timoshenko_basis_regime_extension_status.md).

The evidence is recorded under
`results/article_epsilon_upper_envelope/optimization_benchmark/smoke_failure_audit/`.
The preserved observations below explain the original blocker; current
full-grid solver readiness is `PASS`, although a grid run still requires
separate user authorization.

## Current characteristic equation and supported regime

For one rod section, `scripts/lib/variable_length_timoshenko.py::timo_basis`
forms the quadratic polynomial in the squared spatial root `z`:

```text
B z^2 + c2 z + c0 = 0,
c2 = Omega^2 (B m / K + J),
c0 = J m Omega^4 / K - m Omega^2,
K = kappa G A,
B = E I,
m = rho A,
J = rho I.
```

Here `Omega = epsilon Lambda^2`. The current implementation requires one
positive and one negative real root,

```text
z_+ = a^2 > 0,
z_- = -b^2 < 0,
```

and uses one hyperbolic pair and one trigonometric pair. Its reduced clamped-end
basis is written with `cosh(a x)`, `sinh(a x)`, `cos(b x)`, and `sin(b x)`,
with the existing couplings `h` and `q`. This is the regime actually supported
by the current `TimoshenkoBasis(a, b, h, q)` data structure and endpoint-column
formulas.

## Regime boundary and missing basis

The first boundary occurs at

```text
c0 = 0,
Omega_c = sqrt(K / J),
Lambda_c = sqrt(Omega_c / epsilon).
```

Below `Omega_c`, `c0 < 0`, so the two real `z` roots have opposite signs. Above
`Omega_c`, `c0 > 0`; with `B > 0` and `c2 > 0`, the real roots observed in the
audited parameter range are both negative:

```text
z_1 = -b_1^2 < 0,
z_2 = -b_2^2 < 0.
```

The spatial characteristic roots therefore change from

```text
{+a, -a, +i b, -i b}
```

to

```text
{+i b_1, -i b_1, +i b_2, -i b_2}.
```

The missing regime needs two trigonometric pairs. Passing those roots through
the current `sqrt(z_positive)` / `sqrt(-z_negative)` representation is not an
equivalent numerical rearrangement: the type of the spatial basis has changed.
Consequently, changing it merely to make the validation gate pass would be a
scientific implementation extension and is outside the current task.

## Physical meaning

The boundary is the Timoshenko critical frequency already used as a diagnostic
in the repository. Crossing it is not a physical EB/Timoshenko error threshold
and must not be counted as `delta_f > 0.10`. It marks entry into the second
Timoshenko wave/spectrum regime. The local literature notes identify T1 as the
main experimental/theoretical critical-frequency source and T3/T4 as cautions
about above-critical behavior and the second spectrum; see
`docs/literature/timoshenko_shear_sources.md`. The repository currently contains
the mixed-pair formulas but does not contain a verified coupled-rod endpoint
matrix for the two-trigonometric-pair regime.

## Minimal proposed extension

The smallest reviewable implementation would:

1. represent the spatial roots explicitly, including the regime identifier;
2. add a two-trigonometric-pair endpoint basis derived from the unchanged
   Timoshenko field equations and the same clamped-end reduction;
3. add a cutoff-limit representation that remains finite as one `z` root tends
   to zero;
4. dispatch by the signs of the verified real `z` roots without changing the
   characteristic matrix assembly, joint conditions, unknown ordering, or
   tolerances;
5. retain the current mixed basis bit-for-bit away from the boundary.

Before coding, the two-trigonometric endpoint columns and their `w`, `psi`,
`Q`, and `M` couplings must be derived explicitly and checked against the local
literature PDFs. A solver-aware alternative that stops below the first cutoff
could certify a low-frequency prefix only if a separate root-count proof shows
that the K10/right-guard spectrum is complete below that boundary; it would not
make the basis implementation generally above-cutoff capable.

## Required verification before promotion

Any future implementation must include all of the following:

- symbolic substitution of both basis regimes into the unchanged Timoshenko
  differential equations;
- determinant equivalence in the currently supported mixed regime;
- continuity of determinant, roots, and endpoint columns through `c0 = 0`,
  using a nonsingular cutoff-limit representation;
- thin-limit regression;
- unchanged S3_12 and S3_14 roots, `N_true`, and strict status;
- all four targeted smoke cases with at least 14 primary and 14 strict
  candidate slots;
- before/after root, multiplicity, cluster, and strict-gateway comparisons;
- confirmation that the circular shear coefficient and all existing quality
  tolerances remain unchanged.

## Implementation outcome

The regime-complete evaluator, analytic cutoff limit, independent state-space
oracle, EB-equivalent stable columns, and root-inventory fixes now pass the
versioned 24-point validation. See the status note and
`results/article_epsilon_upper_envelope/solver_readiness_v2/`. No production
runner default or automatic full-grid permission follows from that finite PASS.
