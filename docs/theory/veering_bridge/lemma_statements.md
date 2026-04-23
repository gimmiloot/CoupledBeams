# Lemma Statements

## Scope

This file collects clean statement-level versions of the current bridge lemmas.
Definitions and conventions are taken from [foundations.md](foundations.md).

## 1. Packet construction in the simple-root setting

Fix a regular local domain `Omega x I` and assume:

1. `Omega` is an isolated two-root cluster for `det M(Lambda,p)` on `I`;
2. both roots in `Omega` are simple for `p in I`;
3. the two roots admit locally regular right and left null-vector
   representatives;
4. the two right directions and the two left directions remain linearly
   independent;
5. the retained pairing is nondegenerate.

Then, after shrinking `I` if necessary, one obtains local rank-2 packets

```text
V(p), W(p) in C^{6x2}
```

associated with the isolated two-root cluster.

This is a simple-root conditional packet lemma only. It does not treat
defective roots, multiple roots, complements, or reduced-block normal forms.

## 2. Complement regularity lemma

Assume a local rank-2 packet `V(p),W(p)` has already been supplied on a
regular local domain. Choose complements `Z(p),U(p)` so that `[V Z]` and
`[W U]` are invertible local bases, and define

```text
D(Lambda,p) = U(p)^T M(Lambda,p) Z(p).
```

If

```text
det D(Lambda,p) != 0
```

on the claimed local region, then the complement is regular there and Schur
elimination is locally meaningful.

This statement separates algebraic basis completion from the extra
regular-complement hypothesis. Basis completion alone does not imply
regularity of `D`.

## 3. Schur root-capture lemma

Assume:

1. compatible right/left data `V,W` and complements `Z,U` are supplied on a
   regular local domain;
2. `[V Z]` and `[W U]` are invertible on the claimed local region;
3. with

   ```text
   M_tilde = [W U]^T M [V Z] = [A B; C D],
   ```

   one has `det D(Lambda,p) != 0` on that region.

Define the Schur block

```text
S(Lambda,p) = A(Lambda,p) - B(Lambda,p) D(Lambda,p)^(-1) C(Lambda,p).
```

Then `det M(Lambda,p)` and `det S(Lambda,p)` have the same local zeros on the
claimed region. For fixed `p`, the same one-variable zero multiplicities in
`Lambda` are preserved.

This is exact local root capture, but only on the region where the basis
factors and `det D` remain nonzero.

## 4. Lambda-normalization lemma

Fix a base point `(Lambda0,p0)` in the claimed local region and assume:

1. the Schur root-capture step has already supplied a local analytic reduced
   block `S(Lambda,p)`;
2. the spectral derivative is nonsingular at the base point:

   ```text
   det(partial_Lambda S(Lambda0,p0)) != 0.
   ```

Then, after shrinking the local region if necessary, there exist analytic
`2x2` matrices `R(Lambda,p)` and `K(Lambda,p)` such that

```text
S(Lambda,p) = R(Lambda,p) ((Lambda-Lambda0) I - K(Lambda,p)),
det R(Lambda,p) != 0.
```

Consequently the exact normalized block

```text
F(Lambda,p) = ((Lambda-Lambda0) I - K(Lambda,p))
```

has the same local zeros as `S`, and therefore the same retained local zeros
as `M`.

## 5. First-order reduced model lemma

Assume the exact normalized block already exists on a local region:

```text
S(Lambda,p) = R(Lambda,p) ((Lambda-Lambda0) I - K(Lambda,p)),
det R(Lambda,p) != 0,
```

and `K(Lambda,p)` is analytic in `Lambda` there.

Define the frozen `p`-only matrix

```text
K0(p) = K(Lambda0,p).
```

Then there exists an analytic matrix `K1(Lambda,p)` such that

```text
K(Lambda,p) = K0(p) + (Lambda-Lambda0) K1(Lambda,p),
```

and hence

```text
((Lambda-Lambda0) I - K(Lambda,p))
=
((Lambda-Lambda0) I - K0(p))
- (Lambda-Lambda0) K1(Lambda,p).
```

On a sufficiently small local region, the remainder is bounded by

```text
O(|Lambda-Lambda0|)
```

at matrix level.

This is a controlled first-order approximation statement, not an exact
spectral equivalence theorem.

## 6. Frozen-model root-comparison lemma

Fix one parameter value `p` and define

```text
f_p(Lambda) = det F(Lambda,p),
g_p(Lambda) = det G(Lambda,p),
```

where

```text
F(Lambda,p) = ((Lambda-Lambda0) I - K(Lambda,p)),
G(Lambda,p) = ((Lambda-Lambda0) I - K0(p)).
```

Let `D` be a closed disc in the `Lambda`-plane inside the claimed local
normalization region. Assume:

1. `g_p` has exactly one zero in `D`, and that zero is simple;
2. `g_p` has no zeros on `partial D`;
3. on `partial D`,

   ```text
   sup_{partial D} |f_p - g_p|
   <
   inf_{partial D} |g_p|.
   ```

Then `f_p` has exactly one zero in `D`, counted with multiplicity. Hence the
exact normalized block has one nearby root in the same disc, unique and
simple in that disc.

This is a narrow fixed-`p` local comparison theorem. It is not a global
continuation statement in `p`.

## 7. Uniform parameterwise moving-disc root-comparison lemma

Let `I` be a local interval in the selected scalar sweep parameter `p`, and
assume:

1. there is a local frozen-model root branch `Lambda_fr(p)` on `I`;
2. one fixed radius `r > 0` is chosen for all `p in I`;
3. every disc

   ```text
   D_r(p) = { Lambda : |Lambda - Lambda_fr(p)| <= r }
   ```

   lies inside the claimed local normalization region;
4. for each `p in I`, the frozen determinant `g_p` has exactly one zero in
   `D_r(p)`, and that zero is simple;
5. for each `p in I`, `g_p` has no zeros on `partial D_r(p)`;
6. the uniform boundary-domination inequality

   ```text
   sup_{p in I} sup_{Lambda in partial D_r(p)} |f_p(Lambda) - g_p(Lambda)|
   <
   inf_{p in I} inf_{Lambda in partial D_r(p)} |g_p(Lambda)|
   ```

   holds.

Then, for each `p in I`, the exact normalized determinant `f_p` has exactly
one zero in `D_r(p)`, counted with multiplicity. Hence for each `p in I`, the
exact normalized block has one root in `D_r(p)`, unique and simple in that
disc.

This is a parameterwise family of local fixed-`p` comparison statements with
a uniform radius and uniform boundary control. It is not yet, by itself, a
regular exact root branch theorem.

## 8. Exact-root branch continuation lemma

Assume the hypotheses and conclusion of the uniform parameterwise moving-disc
root-comparison lemma. Fix one base parameter value `p0 in I`, and let
`Lambda_ex(p0)` be the unique exact normalized root selected in `D_r(p0)`.

Assume further that:

1. the exact normalized determinant is jointly analytic in `(Lambda,p)` on a
   neighborhood of `(Lambda_ex(p0),p0)`;
2. the selected exact root is simple:

   ```text
   partial_Lambda f_{p0}(Lambda_ex(p0)) != 0.
   ```

Then, after shrinking to a smaller local interval `I_ex` around `p0` if
necessary, there exists a local regular exact root branch

```text
Lambda_ex(p)
```

on `I_ex` such that

```text
f_p(Lambda_ex(p)) = 0
```

for `p in I_ex`.

Because each disc `D_r(p)` contains exactly one exact normalized root, this
local branch is identified with the moving-disc selected root for every
`p in I_ex`.

This is a local branch-continuation statement for the exact root. It does not
yet quantify the distance between `Lambda_ex(p)` and `Lambda_fr(p)`.

## 9. Quantitative parameterwise branch-comparison lemma

Assume the hypotheses and conclusion of the exact-root branch continuation
lemma on a local interval `I_ex`. Thus:

- `Lambda_fr(p)` is the local frozen-model root branch;
- `Lambda_ex(p)` is the local exact normalized root branch identified by the
  moving-disc selection;
- each moving disc

  ```text
  D_r(p) = { Lambda : |Lambda - Lambda_fr(p)| <= r }
  ```

  contains exactly one exact normalized root, namely `Lambda_ex(p)`.

Assume further that:

1. there exists a constant `m0 > 0` and analytic functions `a_p(Lambda)` on
   `D_r(p)` such that

   ```text
   g_p(Lambda) = a_p(Lambda) (Lambda - Lambda_fr(p)),
   |a_p(Lambda)| >= m0
   ```

   for all `Lambda in D_r(p)` and `p in I_ex`;
2. the local determinant error indicator

   ```text
   eta_loc(p) = sup_{Lambda in D_r(p)} |f_p(Lambda) - g_p(Lambda)|
   ```

   is finite on `I_ex`.

Then, for every `p in I_ex`,

```text
|Lambda_ex(p) - Lambda_fr(p)| <= eta_loc(p) / m0.
```

Consequently, if

```text
eta_* = sup_{p in I_ex} eta_loc(p)
```

is finite, then

```text
sup_{p in I_ex} |Lambda_ex(p) - Lambda_fr(p)| <= eta_* / m0.
```

This is a local quantitative comparison theorem between the exact and frozen
branches. It is conditional on a uniform frozen simple-root lower bound and a
determinant-level error bound. It does not yet provide a sharper asymptotic
law, a derivative formula for the branch shift, or a veering criterion.

## 10. Determinant-discrepancy bound from matrix remainder lemma

Work on a local interval where the exact normalized block and frozen model are
already defined:

```text
F(Lambda,p) = ((Lambda-Lambda0) I - K(Lambda,p)),
G(Lambda,p) = ((Lambda-Lambda0) I - K0(p)),
E(Lambda,p) = F(Lambda,p) - G(Lambda,p).
```

For one fixed `p`, define

```text
f_p(Lambda) = det F(Lambda,p),
g_p(Lambda) = det G(Lambda,p),
```

and use the Frobenius norm on `2x2` matrices.

Then for every `Lambda`,

```text
f_p(Lambda) - g_p(Lambda)
=
tr(adj(G(Lambda,p))) E(Lambda,p) + det(E(Lambda,p)),
```

and consequently

```text
|f_p(Lambda) - g_p(Lambda)|
<=
||G(Lambda,p)||_F ||E(Lambda,p)||_F + (1/2) ||E(Lambda,p)||_F^2.
```

If `D_r(p)` is a moving disc and

```text
M_G(p) = sup_{Lambda in D_r(p)} ||G(Lambda,p)||_F,
M_E(p) = sup_{Lambda in D_r(p)} ||E(Lambda,p)||_F,
```

then

```text
eta_loc(p) = sup_{Lambda in D_r(p)} |f_p(Lambda) - g_p(Lambda)|
<=
M_G(p) M_E(p) + (1/2) M_E(p)^2.
```

Consequently, if the hypotheses of the quantitative parameterwise
branch-comparison lemma also hold, then

```text
|Lambda_ex(p) - Lambda_fr(p)|
<=
[M_G(p) M_E(p) + (1/2) M_E(p)^2] / m0.
```

This is a constructive local quantitative estimate derived from the matrix
remainder of the frozen model. It is still conditional on the earlier moving-
disc and frozen simple-root hypotheses, and it does not yet yield a
derivative-level or asymptotic shift law.
