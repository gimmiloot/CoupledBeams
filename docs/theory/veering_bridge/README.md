# Veering Bridge Package

## Purpose

This package is the canonical local theorem line for the CoupledBeams
veering bridge. It asks when the verified full matrix problem

```text
M(Lambda,p) c = 0
```

admits a local reduced `2x2` characteristic description near an isolated
two-root cluster, and what has actually been proved so far about that line.

## Source of truth

When this package depends on project mathematics, use the project priority:

1. `docs/theory/equations.tex`
2. `docs/theory/assumptions.md`
3. `docs/literature/source_index.md`
4. `src/my_project/analytic/formulas.py`

This package does not modify the verified CoupledBeams matrix entries, sign
conventions, or unknown ordering.

## How to read the package

Read in this order:

1. [foundations.md](foundations.md)
2. [theorem_map.md](theorem_map.md)
3. [lemma_statements.md](lemma_statements.md)
4. [proof_notes.md](proof_notes.md)

`foundations.md` fixes the objects and conventions. `theorem_map.md` is the
main status file. `lemma_statements.md` gives clean statements only.
`proof_notes.md` collects proof attempts, cautions, and failure modes.
`migration_log.md` is a maintenance/history note, not part of the main reading
path.

## What has been reached

The current line has honestly reached:

- exact local reduced-block root capture via Schur elimination, under explicit
  packet/complement hypotheses;
- exact local `Lambda`-normalization of the reduced block, under an explicit
  nonsingular spectral-derivative hypothesis;
- a controlled frozen first-order `p`-only model at matrix level;
- a narrow fixed-`p` nearby-root comparison between the exact normalized block
  and the frozen model, under extra local boundary hypotheses;
- a local parameterwise nearby-root selection over a small interval in `p`
  under uniform moving-disc hypotheses, together with local exact-root branch
  continuation under a simple-root condition;
- a conditional local quantitative comparison estimate between the exact
  branch and the frozen branch, under explicit uniform frozen simple-root and
  determinant-discrepancy hypotheses;
- a constructive determinant-discrepancy bound from the matrix-level
  remainder of the frozen model in the local `2x2` setting;
- a moving-disc closure check for that constructive bound, and a conditional
  local first-order root-shift formula with an explicit remainder bound;
- constructive derivative-control identities and bounds for the quantities
  entering that first-order shift formula in the local `2x2` setting;
- a conditional denominator-control step for the frozen model, reducing the
  remaining derivative lower bound to a frozen-root separation or second-root
  exclusion hypothesis when a local second frozen root branch is available;
- a local coefficient-level frozen discriminant persistence step, together
  with a simple trace/determinant variation criterion that preserves frozen
  noncoalescence on a small interval when the base-point discriminant has
  enough margin;
- a ready-to-use conditional local first-order shift corollary that packages
  moving-disc closure, determinant-discrepancy, derivative-control, and frozen
  discriminant persistence into one local tool;
- a conditional parameterwise derivative-comparison step for the exact and
  frozen branch equations, with the needed denominator and `p`-derivative
  discrepancy controls kept explicit;
- a constructive `p`-derivative determinant-discrepancy bound in the local
  `2x2` setting, replacing the abstract `H_pDelta` input by reduced-model
  matrix controls involving `G`, `E`, `partial_p G`, and `partial_p E`.
- frozen coefficient-derivative identities and Frobenius bounds in the local
  `2x2` setting, identifying `M_{G,p}(p)` with `||K0'(p)||_F` and bounding
  the frozen coefficient derivatives `T'(p)` and `D0'(p)` through `K0(p)` and
  `K0'(p)`;
- a frozen-coefficient ready-to-use derivative-comparison corollary that
  substitutes those bounds into the `RUPDC` derivative estimate while keeping
  `K0'(p)`, `M_{E,p}(p)`, the existing matrix controls, and frozen
  discriminant margins explicit.
- an exact-remainder `p`-derivative control step, reducing `M_{E,p}(p)` to a
  mixed-derivative bound for the exact normalized matrix `K(Lambda,p)` on the
  segments from `Lambda0` to the moving discs;
- an exact-remainder ready-to-use derivative-comparison corollary that
  substitutes this mixed-derivative control into the previous derivative
  estimate.
- a local CoupledBeams quantitative-regime package `Q_der(I_Q)` and a
  conditional derivative-comparison corollary showing that, under that
  explicit regime, `RUPDC-ER` can be cited without restating every remaining
  bound.
- a branch-ready `Q_der` verification protocol for checking one concrete
  candidate interval/packet/window before any branch-case theorem is claimed.
- a grounded C1 theorem-line candidate record, fixing at candidate level the
  base parameter, spectral-window target, and retained-pair target for
  `bending_desc_01` / `bending_desc_04`.

## What has not been reached

The package does not yet prove:

- a project-specific or asymptotic branch-shift law for
  `Lambda_ex(p)-Lambda_fr(p)` with the needed `K0'(p)`, mixed derivative of
  `K(Lambda,p)`, matrix-control, and frozen-discriminant bounds verified;
- verification that the local regime `Q_der(I_Q)` holds for any actual
  CoupledBeams branch or interval;
- packet construction, certified spectral-window endpoints, complement
  regularity, or `Q_der` verification for the current C1 candidate object
  `bending_desc_01` / `bending_desc_04` on the documented `beta=15 deg`,
  physical-radius-`5 mm`, `0.20 <= mu <= 0.35` window;
- further abstract derivative-level progress beyond `RUPDC-ER` without
  project-specific regime assumptions, derivative control from a concrete
  normalization construction, or a different reduction strategy;
- a symmetric or self-adjoint reduced normal form;
- project-defined `delta(p)` and `kappa(p)`;
- a final veering criterion for CoupledBeams;
- any branch-case application theorem.

## Canonical structure

The canonical theorem-line package is:

- `README.md`
- `foundations.md`
- `theorem_map.md`
- `lemma_statements.md`
- `proof_notes.md`

`migration_log.md` records package maintenance history and is not part of the
main theorem line.
