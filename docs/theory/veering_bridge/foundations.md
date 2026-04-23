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
8. not yet reached:
   parameterwise comparison in `p`, symmetric or self-adjoint normal form,
   project-defined `delta-kappa`, final veering criterion, and branch-case
   application theorems.
