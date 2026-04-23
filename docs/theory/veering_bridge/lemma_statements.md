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
