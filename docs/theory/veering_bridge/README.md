# Veering Bridge Package

## Purpose

This package is the canonical local theorem line for the CoupledBeams
veering bridge. It asks when the verified full matrix problem

```text
M(Lambda,p) c = 0
```

admits a local reduced `2x2` characteristic description near an isolated
two-root cluster, and what has actually been proved so far about that line.

The package is intentionally compact. Older bridge notes, registries, and
status files have been absorbed into the canonical files listed below.

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
5. [migration_log.md](migration_log.md)

`foundations.md` fixes the objects and conventions. `theorem_map.md` is the
main status file. `lemma_statements.md` gives clean statements only.
`proof_notes.md` collects proof attempts, cautions, and failure modes.
`migration_log.md` records how the previous file forest was absorbed.

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
  enough margin.

## What has not been reached

The package does not yet prove:

- a derivative-level or asymptotic branch-shift law for
  `Lambda_ex(p)-Lambda_fr(p)`;
- a symmetric or self-adjoint reduced normal form;
- project-defined `delta(p)` and `kappa(p)`;
- a final veering criterion for CoupledBeams;
- any branch-case application theorem.

## Canonical structure

The canonical package is now:

- `README.md`
- `foundations.md`
- `theorem_map.md`
- `lemma_statements.md`
- `proof_notes.md`
- `migration_log.md`

No local archive folder is kept inside the canonical package.

## Editorial policy

This package is canonical and intentionally compact.

### 1. Allowed files

All new bridge mathematics must be integrated into these files only:

- `README.md`
- `foundations.md`
- `theorem_map.md`
- `lemma_statements.md`
- `proof_notes.md`

### 2. No file proliferation

Do not create new narrow notes such as:

- `*_status.md`
- `*_registry.md`
- `*_lemma.md`
- `*_next_step.md`
- local archives inside `veering_bridge/`

Historical preservation should rely on git history, not on a second markdown forest.

### 3. Role of each canonical file

- `foundations.md`: definitions, conventions, notation, layer discipline.
- `theorem_map.md`: what has been reached, status labels, dependencies, limits.
- `lemma_statements.md`: clean statement-level theorem/lemma/corollary text only.
- `proof_notes.md`: proof attempts, cautions, failure modes, constructive remarks.
- `README.md`: short overview of the package and current global status.

### 4. Exact vs approximate discipline

Every new result must be classified explicitly as one of:

- exact algebra / exact local reformulation;
- conditional theorem under explicit hypotheses;
- controlled approximation;
- accepted standard background.

Do not mix exact root-capture statements with approximate frozen-model statements.

### 5. Scope discipline

Do not introduce symmetric-normal-form, `delta/kappa`, veering-criterion, or
branch-application language before the package has actually reached that layer.

### 6. Update rule for every new step

Each genuine new step must update, in this order:

1. `foundations.md` — only if new objects or notation are needed;
2. `lemma_statements.md` — clean statement;
3. `theorem_map.md` — status, dependencies, exact limits;
4. `proof_notes.md` — proof line and failure modes;
5. `README.md` — only if the package-level overview has genuinely changed.

If a step does not change one of these layers, do not edit that file.
