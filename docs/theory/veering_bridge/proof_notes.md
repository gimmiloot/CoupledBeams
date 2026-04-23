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

## 7. Ideal-model algebra stays downstream

The earlier `ideal symmetric 2x2` algebra remains useful only as an ideal
model. It must not be transferred to the CoupledBeams bridge line until a
later symmetric or self-adjoint reduced normal form is actually justified.

Therefore the current package does not claim:

- project-defined `delta(p)` and `kappa(p)`;
- ideal gap formulas for the actual reduced block;
- a direct bridge from the present reduced block to a final veering
  criterion.

## 8. Constructive remarks

The most realistic constructive packet path remains:

1. isolate a two-root cluster;
2. build pointwise right/left root data in the simple-root setting;
3. promote those data to rank-2 packets;
4. test complement regularity;
5. only then use Schur reduction and normalization.

This is still a constructive program, not a theorem that any named packet in
CoupledBeams already satisfies the hypotheses.

## 9. Current next obstruction

The bridge line has reached exact local reduced-block statements plus a narrow
fixed-`p` comparison with the frozen model. The next real obstruction is:

```text
parameterwise comparison in p.
```

That step must be handled before any honest move to:

- symmetric or self-adjoint reduced normal form;
- project-defined `delta-kappa`;
- final veering criteria;
- branch-case application theorems.
