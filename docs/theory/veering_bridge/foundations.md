# Foundations

## Scope

This file collects the definitions and conventions actually used by the
current bridge lemmas. It fixes the working vocabulary for the local theorem
line without changing `docs/theory/equations.tex` or
`docs/theory/assumptions.md`.

## Variable and basis convention

The current bridge line uses the following convention.

- `Lambda` is the spectral variable.
- `p` is one selected scalar sweep parameter.
- All other physical parameters are fixed in any local theorem statement.
- The packet and complement matrices

  ```text
  V(p), W(p), Z(p), U(p)
  ```

  depend on `p` only, not on `Lambda`.

This convention matters because the determinant factors coming from
`[V Z]` and `[W U]` then depend only on `p`, so they do not create or remove
zeros as functions of `Lambda`.

`Lambda`-dependent basis choices are not excluded in principle, but they would
belong to a different and stronger framework than the one used here.

## Ambient full matrix

The ambient object is the verified clamped CoupledBeams matrix

```text
M(Lambda,beta,mu,eps)
```

with frozen unknown order

```text
(A_1, B_1, A_2, B_2, P_1, P_2).
```

In local theorem language, write `M(Lambda,p)` after choosing one scalar sweep
parameter `p` and fixing the others.

## Regular local domain

A regular local domain is a local product region

```text
Omega x I
```

where:

1. `Omega` is a bounded local spectral window in the `Lambda`-plane;
2. `I` is a local interval in the chosen scalar parameter `p`;
3. the matrix entries of `M(Lambda,p)` are analytic on a neighborhood of
   `Omega x I`;
4. collapsed or degenerate parameter loci needed to keep the local theorem
   meaningful are excluded.

For the current CoupledBeams matrix, this means in particular:

- exclude collapsed arm lengths such as `1-mu=0` or `1+mu=0` when `p=mu`;
- keep `Omega` bounded;
- keep the region inside the local analytic regime actually used by the
  statement.

Regularity alone does not give isolated clusters, packet data, complements,
or normalized reduced blocks.

## Isolated two-root cluster

For `p in I`, a spectral window `Omega` is an isolated two-root cluster if:

1. `det M(Lambda,p)` has exactly two roots in `Omega`, counted with algebraic
   multiplicity;
2. `det M(Lambda,p)` has no roots on `partial Omega`;
3. those two roots are separated from all roots outside `Omega`.

This is scalar characteristic information only.

## Local right packet

A local right packet is a matrix

```text
V(p) in C^{6x2}
```

with rank two on the local interval, whose column span represents the retained
right-side packet associated with the isolated two-root cluster.

In the simple-root setting, the columns may be chosen from right null vectors
at the two roots. In more general settings, `V(p)` is only a packet-level
retained space unless additional structure is proved.

## Local left packet

A local left packet is a matrix

```text
W(p) in C^{6x2}
```

with rank two on the local interval, associated with the same isolated
two-root cluster and used as the retained left test packet.

## Compatible right/left data

Compatible right/left data consist of a local right packet `V(p)` and a local
left packet `W(p)` such that:

1. both are tied to the same isolated two-root cluster;
2. `rank V(p) = rank W(p) = 2`;
3. they vary with the local regularity required by the theorem;
4. the retained pairing is nondegenerate:

   ```text
   det(W(p)^T V(p)) != 0;
   ```

5. they can be extended to full right/left bases by complements.

This is a theorem-level working object, not a proof that such data exist for
any concrete CoupledBeams packet.

## Regular complement

Given compatible right/left data, choose complements

```text
Z(p), U(p) in C^{6x4}
```

so that `[V Z]` and `[W U]` are invertible local bases. Form

```text
M_tilde(Lambda,p) =
[W U]^T M(Lambda,p) [V Z]
=
[ A  B
  C  D ].
```

The complement is regular on the claimed local region if

```text
det D(Lambda,p) != 0
```

there.

This is stronger than basis completion alone. Invertibility of `[V Z]` and
`[W U]` does not imply invertibility of `D`.

## Projected Schur block

When the complement is regular, the projected Schur block is

```text
S(Lambda,p) = A(Lambda,p) - B(Lambda,p) D(Lambda,p)^(-1) C(Lambda,p).
```

This is the local retained `2x2` characteristic block.

On the claimed local region, the Schur root-capture step identifies the zeros
of `det S` with the retained zeros of `det M`, up to nonzero determinant
factors.

## Exact normalized block

Fix a base point `(Lambda0,p0)` in the local region. If the local reduced
block satisfies the nonsingular spectral-derivative hypothesis

```text
det(partial_Lambda S(Lambda0,p0)) != 0,
```

then, after shrinking if necessary, there exist analytic `2x2` matrices
`R(Lambda,p)` and `K(Lambda,p)` such that

```text
S(Lambda,p) = R(Lambda,p) ((Lambda-Lambda0) I - K(Lambda,p)),
det R(Lambda,p) != 0.
```

The block

```text
F(Lambda,p) = ((Lambda-Lambda0) I - K(Lambda,p))
```

is the exact normalized block.

This is an exact local reformulation up to a nonzero analytic left factor.

## Frozen `p`-only model

From the exact normalized block, define

```text
K0(p) = K(Lambda0,p),
G(Lambda,p) = ((Lambda-Lambda0) I - K0(p)).
```

`G` is the frozen `p`-only first-order reduced model.

It is not exactly equivalent to the exact normalized block in general. It is a
controlled approximation with matrix-level remainder

```text
F(Lambda,p) - G(Lambda,p) = O(|Lambda-Lambda0|).
```

## Fixed-`p` exact and frozen determinants

For one fixed parameter value `p`, define

```text
f_p(Lambda) = det F(Lambda,p),
g_p(Lambda) = det G(Lambda,p).
```

These are one-variable analytic functions of `Lambda`.

Local frozen-model root comparison is stated in terms of `f_p` and `g_p` on a
small closed disc `D` in the `Lambda`-plane.

## Local frozen-model root branch

A local frozen-model root branch is a function

```text
Lambda_fr(p)
```

defined on a local interval `I`, such that

```text
g_p(Lambda_fr(p)) = 0
```

for each `p in I`.

When this object is used in the current bridge line, it is understood as a
root branch of the frozen `p`-only model, not of the exact normalized block.
The simple-root condition for `g_p` is stated separately when needed.

## Local second frozen root branch and frozen-root separation

When the frozen quadratic has a second local root branch on an interval `I`,
denote it by

```text
Lambda_fr,2(p),
```

with

```text
g_p(Lambda_fr,2(p)) = 0,
Lambda_fr,2(p) != Lambda_fr(p).
```

The corresponding frozen-root separation is

```text
sep_fr(p) = |Lambda_fr(p) - Lambda_fr,2(p)|.
```

A uniform frozen-root noncoalescence bound on `I` means that there exists a
constant `sigma_fr > 0` such that

```text
sep_fr(p) >= sigma_fr
```

for all `p in I`.

If the second frozen root stays outside the moving disc `D_r(p)`, then
automatically

```text
sep_fr(p) > r.
```

## Frozen discriminant of the frozen `2x2` model

For the frozen block

```text
G(Lambda,p) = ((Lambda-Lambda0) I - K0(p)),
```

set

```text
T(p) = tr K0(p),
D0(p) = det K0(p).
```

the frozen determinant is the monic quadratic

```text
g_p(Lambda)
=
(Lambda-Lambda0)^2 - T(p) (Lambda-Lambda0) + D0(p).
```

Its frozen discriminant is

```text
Disc_fr(p) = T(p)^2 - 4 D0(p).
```

After shrinking to a sufficiently small interval on which the two frozen roots
can be labeled locally, one has

```text
Disc_fr(p) = (Lambda_fr(p) - Lambda_fr,2(p))^2,
|Disc_fr(p)| = sep_fr(p)^2.
```

## Moving discs

Given a local frozen-model root branch `Lambda_fr(p)` and a radius `r > 0`,
define the moving closed discs

```text
D_r(p) = { Lambda : |Lambda - Lambda_fr(p)| <= r }.
```

The current parameterwise comparison step uses one uniform radius `r` for all
`p` in the chosen local interval.

## Local exact root branch

A local exact root branch is a locally regular function

```text
Lambda_ex(p)
```

defined on a local interval, such that

```text
f_p(Lambda_ex(p)) = 0
```

for each `p` in that interval.

In the current bridge line, `Lambda_ex(p)` is selected from the exact
normalized block and may later be identified with the moving-disc selected
root when uniqueness in each disc is known.

## Local determinant error indicator

Once the exact normalized block and the frozen `p`-only model are both fixed,
define the local determinant discrepancy on a moving disc by

```text
eta_loc(p) = sup_{Lambda in D_r(p)} |f_p(Lambda) - g_p(Lambda)|.
```

Its uniform version on an interval `I` is

```text
eta_* = sup_{p in I} eta_loc(p).
```

This is a scalar determinant-level error indicator. It is not the same object
as the earlier matrix-level remainder for `F-G`, and it does not coincide
with that remainder automatically. The latter may only be used to bound
`eta_loc` through an additional determinant-discrepancy estimate.

## Matrix-level remainder and local matrix controls

Once the exact normalized block and the frozen `p`-only model are both fixed,
set

```text
E(Lambda,p) = F(Lambda,p) - G(Lambda,p).
```

For determinant-discrepancy estimates in the current local `2x2` setting, the
package uses the Frobenius norm. On a moving disc `D_r(p)`, define

```text
M_G(p) = sup_{Lambda in D_r(p)} ||G(Lambda,p)||_F,
M_E(p) = sup_{Lambda in D_r(p)} ||E(Lambda,p)||_F,
M_{E,1}(p) = sup_{Lambda in D_r(p)} ||partial_Lambda E(Lambda,p)||_F.
```

These are local matrix-level control quantities on the moving discs. They are
used to bound `eta_loc(p)` and the derivative discrepancy
`partial_Lambda(f_p-g_p)` through exact `2x2` determinant identities.

## Local constructive branch-shift bound and closure condition

Once the frozen simple-root lower bound `m0` and the local matrix controls
`M_G(p), M_E(p)` are available, define

```text
B_shift(p) = [M_G(p) M_E(p) + (1/2) M_E(p)^2] / m0.
```

Its uniform version on an interval `I` is

```text
B_* = sup_{p in I} B_shift(p).
```

The corresponding moving-disc closure condition is

```text
B_shift(p) < r
```

for one fixed parameter value, or

```text
B_* < r
```

uniformly on an interval.

This is not a new existence theorem. It is a self-consistency condition for
the already chosen moving discs: if it holds, then the current quantitative
branch bound itself guarantees that the exact root stays inside `D_r(p)`.

## Uniform frozen simple-root nondegeneracy

Let `Lambda_fr(p)` be a local frozen-model root branch on an interval `I`, and
let `r > 0` be the moving-disc radius. Uniform frozen simple-root
nondegeneracy on the discs `D_r(p)` means that there exist analytic functions
`a_p(Lambda)` and a constant `m0 > 0` such that

```text
g_p(Lambda) = a_p(Lambda) (Lambda - Lambda_fr(p)),
|a_p(Lambda)| >= m0
```

for all `Lambda in D_r(p)` and all `p in I`.

This is the quantitative form of simple-root isolation used in the current
branch-comparison step. It is equivalent in spirit to a uniform derivative
lower bound for the frozen determinant near the frozen branch, after shrinking
the interval and discs if necessary.

## Uniform boundary-domination language

Let `Lambda_fr(p)` be a local frozen-model root branch on an interval `I`, and
let `r > 0` be fixed. Uniform moving-disc boundary domination means:

1. every disc `D_r(p)` lies inside the claimed local normalization region;
2. `g_p` has no zeros on `partial D_r(p)` for `p in I`;
3. the strict bound

   ```text
   sup_{p in I} sup_{Lambda in partial D_r(p)} |f_p(Lambda) - g_p(Lambda)|
   <
   inf_{p in I} inf_{Lambda in partial D_r(p)} |g_p(Lambda)|
   ```

   holds.

This is stronger than pointwise domination at one parameter value. It is the
uniform hypothesis that lets the same moving-disc comparison argument be
applied for every `p` in the interval.

## Multiplicity convention

Unless explicitly stated otherwise, multiplicity means:

```text
the order of a zero of a one-variable analytic function of Lambda at fixed p.
```

This convention is used for:

- local root capture by `det S`;
- exact normalization;
- fixed-`p` frozen-model root comparison.

No multivariable multiplicity notion on the full `(Lambda,p)` region is used
here.

## Exact versus approximate statements

An exact statement in this package means one of the following:

- an exact algebraic identity;
- an exact local reformulation up to invertible analytic factors;
- an exact equality of local zero sets on a claimed region.

An approximate statement means:

- a controlled remainder formula;
- a local matrix-level approximation;
- or a fixed-`p` spectral comparison that requires extra local smallness and
  boundary hypotheses.

Approximate steps must not be read as exact theorems about the full retained
cluster over a parameter interval.

## Layer discipline

The bridge line currently uses the following layer discipline.

1. scalar root information:
   isolated two-root cluster for `det M`;
2. packet data:
   local right/left packet and compatible right/left data;
3. complement regularity:
   complements `Z,U` and `det D != 0` on the claimed region;
4. Schur root capture:
   `det M` and `det S` have the same local retained zeros;
5. exact `Lambda`-normalization:
   `S = R((Lambda-Lambda0)I-K(Lambda,p))`;
6. frozen first-order model:
   `G = ((Lambda-Lambda0)I-K0(p))` with controlled remainder;
7. fixed-`p` root comparison:
   a simple frozen-model root carries a nearby exact normalized root under
   extra local boundary hypotheses;
8. parameterwise nearby-root comparison:
   a uniform moving-disc argument selects one nearby exact root for each `p`
   in a local interval, and simple exact-root continuation upgrades that
   selection to a local exact branch;
9. quantitative parameterwise comparison:
   under a uniform frozen simple-root lower bound and a determinant error
   indicator, the exact branch is quantitatively close to the frozen branch;
10. constructive determinant-discrepancy control:
    in the local `2x2` setting, the determinant discrepancy can be bounded
    explicitly in terms of the matrix-level remainder of the frozen model;
11. moving-disc closure:
    the constructive branch bound may be checked against the chosen disc
    radius as a self-consistency condition;
12. first-order root-shift comparison:
    under explicit derivative lower bounds and derivative regularity on the
    moving discs, one gets a conditional local first-order shift formula with
    a controlled remainder;
13. derivative control from reduced-model objects:
    the derivative quantities entering the first-order shift formula are
    linked constructively to `G`, `E`, and `partial_Lambda E` on the moving
    discs;
14. frozen-model denominator control from noncoalescence:
    when a local second frozen root branch is available, the remaining
    denominator in the first-order shift formula is identified with the
    frozen-root separation and is therefore controlled by frozen noncoalescence
    or second-root exclusion from the moving disc;
15. frozen discriminant persistence:
    a nonzero frozen discriminant at a base point persists on a smaller local
    interval and, after local labeling of the frozen roots, gives a uniform
    local denominator lower bound there;
16. frozen discriminant persistence from coefficient control:
    explicit bounds on `T(p)` and `D0(p)` can keep the frozen discriminant
    away from zero on an interval and thereby feed the local denominator
    lower bound after shrinking;
17. not yet reached:
    sharper asymptotic or derivative-level branch-shift laws, symmetric or
    self-adjoint normal form, project-defined `delta-kappa`, final veering
    criterion, and branch-case application theorems.
