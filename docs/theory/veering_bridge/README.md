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
  and the frozen model, under extra local boundary hypotheses.

## What has not been reached

The package does not yet prove:

- parameterwise comparison or continuation in `p`;
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
