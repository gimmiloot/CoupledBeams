# Proof Notes

## Scope

This file is the non-canonical working companion to the clean statements in
[lemma_statements.md](lemma_statements.md).
It records proof attempts, status distinctions, failure modes, and the main
reasons the bridge line can break.

## Status discipline

The package uses five labels:

- `proved_here`
- `accepted_standard_background`
- `explicit_hypothesis`
- `needs_caution`
- `not_reached_yet`

`proved_here` means the package contains the relevant local argument.
`accepted_standard_background` means the step is mathematically standard but
not reproved in full detail inside this package.

## 1. Packet construction in the simple-root setting

### What is actually proved here

- `det M(Lambda_j,p)=0` implies singularity of the square matrix
  `M(Lambda_j,p)`;
- once two right directions and two left directions are supplied at fixed `p`
  with rank two, they define rank-2 packets;
- once continuous representatives exist and the relevant determinants are
  nonzero at a base point, rank and pairing persist after shrinking.

### What is only accepted standard background

- simple determinant root implies one-dimensional right and left kernels;
- local continuation of simple roots;
- local representative selection for one-dimensional right and left kernel
  lines.

### Why this step can break

- roots cease to be simple;
- the two right or left directions lose independence;
- a continuous packet basis is not justified;
- the retained pairing degenerates.

### Why this is still useful

This is the first honest step from scalar root data to packet data. It is
deliberately narrow and does not yet touch complements.

## 2. Complement regularity

### The key separation

Two facts must not be conflated:

1. algebraic basis completion:

   ```text
   V,W -> [V Z], [W U];
   ```

2. regular complement:

   ```text
   det D(Lambda,p) != 0,
   D = U^T M Z.
   ```

The first is ordinary finite-dimensional linear algebra. The second is the
actual theorem-level regularity hypothesis.

### What is proved here

- basis completion exists at fixed `p`;
- the block matrix `M_tilde = [W U]^T M [V Z]` is well defined once the bases
  are supplied;
- basis completion does not imply invertibility of `D`;
- base-point invertibility of `D` gives only local invertibility after
  shrinking.

### Main failure modes

- `D` is already singular at the base point;
- `D` is invertible near the base point but loses invertibility on the claimed
  window;
- one silently treats basis completion as if it already implied regular
  complement.

## 3. Schur root capture

### Three layers

1. block determinant identity:

   ```text
   det M_tilde = det D det S;
   ```

2. basis-factor relation:

   ```text
   det M_tilde = det([W U]^T) det M det([V Z]);
   ```

3. local root capture:
   once all prefactors are nonzero on the claimed region, the local zeros of
   `det M` are exactly the zeros of `det S`.

### What is exact

The determinant factorization is exact once the packet, complements, and
regular complement hypothesis are supplied.

### Caution on multiplicity

Multiplicity is only discussed in the fixed-`p`, one-variable sense in
`Lambda`. The package does not use any multivariable multiplicity notion on
`Omega x I`.

### Main failure modes

- `det D` vanishes somewhere in the claimed region;
- basis determinants vanish;
- the local theorem is overstretched beyond the region where the nonzero
  factors are controlled.

## 4. Exact Lambda-normalization

### What is proved here

Given analytic `S(Lambda,p)` and the explicit hypothesis

```text
det(partial_Lambda S(Lambda0,p0)) != 0,
```

the exact local factorization

```text
S(Lambda,p) = R(Lambda,p) ((Lambda-Lambda0) I - K(Lambda,p))
```

is available after shrinking if needed.

### Why this matters

This is the first exact eigenvalue-style reduced form. It is still general:
`K` may depend on both `Lambda` and `p`.

### What it is not

It is not yet:

- a `p`-only model;
- a symmetric normal form;
- a `delta-kappa` representation.

### Main failure modes

- `partial_Lambda S(Lambda0,p0)` is singular;
- the invertible spectral coefficient is lost on a window larger than the
  region where the proof really works;
- a first-order model is mistaken for this exact normalization.

## 5. First-order freezing

### Exact versus approximate

The exact normalized block is

```text
F(Lambda,p) = ((Lambda-Lambda0) I - K(Lambda,p)).
```

The frozen model is

```text
G(Lambda,p) = ((Lambda-Lambda0) I - K0(p)),
K0(p) = K(Lambda0,p).
```

The package proves only the exact remainder formula

```text
F - G = -(Lambda-Lambda0) K1(Lambda,p),
```

with matrix-level size `O(|Lambda-Lambda0|)`.

### What this step does not prove

It does not prove exact root capture for the frozen model. It only provides a
controlled matrix approximation.

### Main failure modes

- `K` is not regular enough in `Lambda`;
- the local region is too large for the remainder to be informative;
- one silently transfers exact spectral statements from `F` to `G`.

## 6. Fixed-`p` root comparison

### The crucial distinction

Matrix-level approximation does not by itself imply spectral approximation.

To compare roots of `det F` and `det G`, the package fixes one parameter value
`p`, works with one-variable analytic determinants

```text
f_p = det F(.,p),
g_p = det G(.,p),
```

and adds a genuinely spectral local hypothesis on a chosen disc.

### What is proved here

- for `2x2` matrices with `F = G + E`,

  ```text
  det(G+E) - det(G) = tr(adj(G)E) + det(E);
  ```

- therefore the earlier matrix remainder gives a determinant perturbation
  bound on the boundary of a small disc;
- if the frozen determinant has one simple zero in that disc, no boundary
  zero, and the boundary domination inequality holds, then the exact
  normalized determinant has one nearby root in the same disc.

### What is only standard background

The step from boundary domination to equal zero count is imported as a
Rouche-type one-variable analytic root theorem.

### What this still does not give

- any global branch theorem in `p`;
- a quantitative root-shift bound;
- any veering interpretation.

## 7. Parameterwise moving-disc comparison and exact-branch continuation

### Pointwise-in-`p` Rouche versus branch regularity

Two layers must stay separate.

First layer:

```text
for each fixed p,
one Rouche-type comparison on one disc.
```

Second layer:

```text
as p varies,
those selected exact roots fit together into one regular branch.
```

The first layer is a parameterwise family of fixed-`p` arguments. The second
layer needs a separate simple-root continuation theorem.

### Shifted variable around the moving center

Let `Lambda_fr(p)` be a local frozen-model root branch of the frozen model and
set

```text
xi = Lambda - Lambda_fr(p).
```

Then the moving disc condition is simply

```text
|xi| <= r.
```

This makes clear that the real new object is not one fixed disc around one
fixed center, but one family of translated discs moving with the frozen root.

### What is genuinely new here

The real new ingredients are:

- one uniform radius `r > 0`;
- one local interval `I` in the scalar parameter `p`;
- one frozen-model root branch `Lambda_fr(p)`;
- one uniform boundary-domination estimate on all disc boundaries.

Pointwise domination at `p0` alone is not enough. The new theorem step needs a
single radius and a single strict domination margin that survive over a whole
local interval in `p`.

### Why "one root in each disc" is not yet a regular branch

Even if every moving disc contains exactly one exact root, that statement is
still set-theoretic in `p`: it gives a selected root value for each parameter,
but not yet a regular parametrization.

To obtain a local branch `Lambda_ex(p)`, one still needs:

- joint analyticity in `(Lambda,p)`;
- a simple exact root at the base point;
- an implicit-function-style continuation step.

### How uniqueness and continuation combine

The logical order is:

1. uniform moving-disc comparison gives one exact root in each disc;
2. at `p0`, if that exact root is simple, implicit-function-style continuation
   yields a local regular exact branch;
3. uniqueness in each disc prevents ambiguity, so the continued branch must be
   exactly the moving-disc selected root for nearby `p`.

This is the clean bridge from fixed-`p` nearby-root comparison to a local
parameterwise exact-root branch.

### Main failure modes

- one cannot choose a single uniform radius `r`;
- the frozen root is not isolated enough to support one moving disc family;
- domination holds at `p0` but not uniformly over the interval;
- the selected exact root loses simplicity;
- a moving disc captures competing exact or frozen roots.

### What this still does not give

This step still does not provide:

- a quantitative estimate for `|Lambda_ex(p)-Lambda_fr(p)|`;
- a global continuation theorem in `p`;
- a symmetric normal form;
- project-defined `delta-kappa`;
- any veering interpretation.

## 8. Quantitative parameterwise branch comparison

### Exact branch versus frozen branch

At this stage the bridge line has four distinct layers that must stay
separate:

1. exact normalized determinant:

   ```text
   f_p(Lambda) = det F(Lambda,p);
   ```
2. frozen-model determinant:

   ```text
   g_p(Lambda) = det G(Lambda,p);
   ```
3. existence and uniqueness in moving discs:
   for each `p`, one exact root lies in `D_r(p)`;
4. quantitative branch comparison:
   estimate the distance between `Lambda_ex(p)` and `Lambda_fr(p)`.

The new step concerns only layer 4. It must not be confused with the earlier
existence/uniqueness theorem.

### Where simple-root nondegeneracy is used

Two different simple-root inputs appear, and they play different roles.

- Exact simple-root nondegeneracy at the base point is used upstream, in the
  exact-root branch continuation step, to produce the local regular branch
  `Lambda_ex(p)`.
- Frozen simple-root nondegeneracy is used here, in quantitative form, to
  control how far `Lambda_ex(p)` can move away from `Lambda_fr(p)`.

The present estimate uses the second input directly.

### Preferred quantitative line

Work on the interval where the earlier moving-disc selection and exact-branch
continuation already hold. For each `p`, assume the frozen determinant admits
the factorization

```text
g_p(Lambda) = a_p(Lambda) (Lambda - Lambda_fr(p)),
|a_p(Lambda)| >= m0 > 0
```

on `D_r(p)`.

Define the determinant discrepancy by

```text
eta_loc(p) = sup_{Lambda in D_r(p)} |f_p(Lambda)-g_p(Lambda)|.
```

Since `Lambda_ex(p)` lies in `D_r(p)` and satisfies `f_p(Lambda_ex(p))=0`,
one has

```text
|g_p(Lambda_ex(p))|
=
|g_p(Lambda_ex(p)) - f_p(Lambda_ex(p))|
<=
eta_loc(p).
```

Using the frozen factorization at `Lambda_ex(p)` gives

```text
|a_p(Lambda_ex(p))| |Lambda_ex(p)-Lambda_fr(p)|
<=
eta_loc(p),
```

hence

```text
|Lambda_ex(p)-Lambda_fr(p)| <= eta_loc(p) / m0.
```

This is the cleanest honest quantitative statement currently available in the
package.

It is still only a pointwise-in-`p` estimate unless `eta_loc(p)` is itself
controlled uniformly. In the present notation, uniform smallness on an
interval is obtained only after controlling

```text
eta_* = sup_{p in I_ex} eta_loc(p).
```

### What is proved here and what is standard background

- `proved_here`:
  once the factorization and lower bound are assumed, the displacement
  estimate follows immediately by evaluating the determinant discrepancy at
  the exact root.
- `accepted_standard_background`:
  the passage from simple-root isolation or a derivative lower bound to the
  factorization

  ```text
  g_p(Lambda) = a_p(Lambda)(Lambda-Lambda_fr(p))
  ```

  with `a_p` nonvanishing on a small neighborhood.

If one prefers to formulate the hypothesis through
`|partial_Lambda g_p(Lambda_fr(p))| >= m0`, then converting that derivative
bound into the factorized lower bound above is standard local analytic-root
theory, not something reproved in full detail here.

### Why the derivative lower bound matters

The estimate divides by the frozen nondegeneracy constant. If the lower bound
for `|a_p|` deteriorates, then the comparison constant blows up and the
quantitative control becomes useless even if nearby-root selection still
holds.

This is exactly where the quantitative step is stricter than the earlier
existence theorem.

### Main failure modes

- the exact branch exists, but the frozen simple-root lower bound is too weak
  or collapses somewhere on the interval;
- the moving-disc radius is large enough that the factorized frozen
  coefficient is no longer uniformly bounded away from zero;
- the determinant discrepancy `eta_loc(p)` is finite but too large to give a
  useful comparison;
- one silently treats this quantitative estimate as if it were already a
  derivative formula or a gap law.

### What this still does not give

This step still does not provide:

- an asymptotic formula for `Lambda_ex(p)-Lambda_fr(p)`;
- a derivative identity for `d/dp (Lambda_ex-Lambda_fr)`;
- a symmetric or self-adjoint reduced normal form;
- project-defined `delta-kappa`;
- any veering criterion.

## 9. Determinant-discrepancy control from the matrix remainder

### Exact identity versus scalar bound

At the matrix level one has

```text
F(Lambda,p) = G(Lambda,p) + E(Lambda,p).
```

For `2x2` matrices, the determinant discrepancy is governed by the exact
identity

```text
det(G+E)-det(G)=tr(adj(G)E)+det(E).
```

This identity is exact. The scalar bound derived from it is a second step.

### Why the Frobenius norm is convenient here

In the current local `2x2` setting, the Frobenius norm gives a clean bound:

```text
|tr(adj(G)E)| <= ||adj(G)||_F ||E||_F = ||G||_F ||E||_F,
```

because for `2x2` matrices one has `||adj(G)||_F = ||G||_F`. Together with

```text
|det(E)| <= (1/2) ||E||_F^2,
```

this yields

```text
|f_p(Lambda)-g_p(Lambda)|
<=
||G(Lambda,p)||_F ||E(Lambda,p)||_F + (1/2)||E(Lambda,p)||_F^2.
```

### Moving-disc consequence

On a moving disc `D_r(p)`, define

```text
M_G(p)=sup_{Lambda in D_r(p)} ||G(Lambda,p)||_F,
M_E(p)=sup_{Lambda in D_r(p)} ||E(Lambda,p)||_F.
```

Then the determinant error indicator satisfies

```text
eta_loc(p)
<=
M_G(p) M_E(p) + (1/2) M_E(p)^2.
```

This is the missing constructive bridge from matrix-level freezing error to
the scalar discrepancy used in the earlier quantitative branch-comparison
lemma.

### Combined branch-distance estimate

Once the earlier `QPBC` hypotheses also hold, substitution into

```text
|Lambda_ex(p)-Lambda_fr(p)| <= eta_loc(p)/m0
```

gives the constructive bound

```text
|Lambda_ex(p)-Lambda_fr(p)|
<=
[M_G(p) M_E(p) + (1/2) M_E(p)^2] / m0.
```

This is still conditional. It does not yet produce a derivative law, because
it only controls the size of the branch shift through local matrix norms.

### Main failure modes

- `M_E(p)` is finite but not small enough for the bound to be useful;
- `M_G(p)` becomes large on the moving discs, weakening the estimate;
- the matrix bound is available only on discs that are too small to support
  the earlier moving-disc comparison cleanly;
- one confuses the exact determinant identity with the derived norm bound;
- one treats this constructive estimate as if it were already an asymptotic
  expansion.

### What this still does not give

This step still does not provide:

- a rate sharper than quadratic-linear control in the local matrix norms;
- a derivative formula for the branch shift;
- a symmetric or self-adjoint reduced normal form;
- project-defined `delta-kappa`;
- a veering criterion.

## 10. Ideal-model algebra stays downstream

The earlier `ideal symmetric 2x2` algebra remains useful only as an ideal
model. It must not be transferred to the CoupledBeams bridge line until a
later symmetric or self-adjoint reduced normal form is actually justified.

Therefore the current package does not claim:

- project-defined `delta(p)` and `kappa(p)`;
- ideal gap formulas for the actual reduced block;
- a direct bridge from the present reduced block to a final veering
  criterion.

## 11. Constructive remarks

The most realistic constructive packet path remains:

1. isolate a two-root cluster;
2. build pointwise right/left root data in the simple-root setting;
3. promote those data to rank-2 packets;
4. test complement regularity;
5. only then use Schur reduction and normalization.

This is still a constructive program, not a theorem that any named packet in
CoupledBeams already satisfies the hypotheses.

## 12. Current next obstruction

The bridge line has reached exact local reduced-block statements, local exact
normalization, controlled freezing, local parameterwise nearby-root
selection, exact-branch continuation, a conditional local quantitative
comparison estimate between the exact and frozen branches, and a constructive
determinant-discrepancy bound from the matrix remainder. The next real
obstruction is:

```text
turn this constructive local branch bound into a sharper branch-shift law or
a more structural reduced description.
```

That step must be handled before any honest move to:

- symmetric or self-adjoint reduced normal form;
- project-defined `delta-kappa`;
- final veering criteria;
- branch-case application theorems.
