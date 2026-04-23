# Theory AGENTS

## Scope
This file applies only to `docs/theory/`.
It governs work on `equations.tex`, `assumptions.md`, derived theory notes, and notation-consistency tables.

## Theory Source Of Truth
For theory work, treat the following as the source of truth, in this order:
1. Locally verified formulas in `docs/theory/equations.tex`
2. Explicit assumptions in `docs/theory/assumptions.md`
3. Source-specific notes fixed in `docs/literature/source_index.md`

If theory and code disagree, do not silently "fix" one from the other. Record the mismatch explicitly first.

## Notation
- Orient base notation to `docs/literature/pdf/Статья-Дорофеев-2025.pdf`.
- Local extensions are allowed only when they are explicitly explained.
- Record any new notation in `docs/theory/assumptions.md` or in a dedicated notation note.
- Do not perform mass renaming of symbols unless the user explicitly requests it.

## Formula Edit Rules
Treat verified determinant-like expressions, matrices, unknown ordering, signs, and coefficients as frozen.

Do not:
- change determinant entries
- change signs or coefficients
- reorder unknowns
- rewrite formulas into an "equivalent" but less transparent form
- hide sign structure through algebraic reshuffling

If a place looks suspicious:
1. Do not fix it immediately.
2. Point to the exact file and line.
3. Compare it against code and literature.
4. Explain why it is suspicious.
5. Only then propose the smallest possible correction.

## Symbolic Checking
- Prefer explicit symbolic checks for nontrivial derivation steps.
- Always check matrix dimensions, unknown ordering, and limiting cases such as `mu = 0`, `beta = 0`, the symmetric limit, and the fixed-fixed limit.
- If Lean is available and the user asks for strict verification, use Lean for proposition-level transitions, reduction claims, and delicate symbolic steps.
- If Lean is unavailable, say so explicitly and use symbolic or reproducible numerical checks instead.
- Do not claim a step is correct without a derivation, a symbolic verification, or a reproducible numerical check.

## Physical Meaning
Check every new hypothesis, simplification, or proposal for physical meaning.

Always ask:
- what physical regime is being described?
- is dimensional consistency preserved?
- are the required boundary conditions still satisfied?
- what happens in the symmetric limit?
- does this contradict verified code?
- is this a physical claim or only a technical convenience?

Record every new assumption in `docs/theory/assumptions.md`.

## Theory And Code
When working on theory, always look for the corresponding implementation in:
- `src/my_project/analytic/`
- `src/my_project/fem/` when relevant

For important formulas, check:
- parameter consistency
- unknown ordering
- determinant structure
- at least one control limiting case

## Literature Warning
`docs/literature/pdf/2003JSVb.pdf` contains a known sign issue in a determinant-like matrix record.
Do not copy its sign pattern blindly into local theory.
When comparing against that source, always account for the fact that some signs in local formulas may already be absorbed by arguments such as `-\Lambda(...)` or by moving terms to the left-hand side.

## Veering bridge package maintenance
This subsection applies to `docs/theory/veering_bridge/`.

Canonical files for `docs/theory/veering_bridge/` are:
- `README.md`
- `foundations.md`
- `theorem_map.md`
- `lemma_statements.md`
- `proof_notes.md`

Do not create file proliferation in this package, including:
- `*_status.md`
- `*_registry.md`
- `*_lemma.md`
- `*_next_step.md`
- local archives inside `docs/theory/veering_bridge/`

Historical preservation must use git history, not a new markdown archive.
The existing `migration_log.md` is a maintenance/history note, not part of the canonical theorem line.
`migration_log.md` update only when the package structure changes materially; do not edit it for ordinary theorem-line progress.

Canonical file roles:
- `foundations.md`: definitions, conventions, notation, layer discipline.
- `theorem_map.md`: what has been reached, status labels, dependencies, limits.
- `lemma_statements.md`: clean statement-level theorem, lemma, and corollary text only.
- `proof_notes.md`: proof attempts, cautions, failure modes, constructive remarks.
- `README.md`: short reader-facing overview of the package and current global status.

Every new result must be classified explicitly as one of:
- exact algebra / exact local reformulation;
- conditional theorem under explicit hypotheses;
- controlled approximation;
- accepted standard background.

Do not mix exact root-capture statements with approximate frozen-model statements.

Do not introduce premature symmetric-normal-form, `delta-kappa`, veering-criterion, or branch-application language before the package has actually reached that layer.

Each genuine new theorem step must update, in this order:
1. `foundations.md` only if new objects are needed.
2. `lemma_statements.md` with the clean statement.
3. `theorem_map.md` with status, dependencies, and exact limits.
4. `proof_notes.md` with the proof line and failure modes.
5. `README.md` only if the package-level overview genuinely changes.

If a ready-to-use corollary follows from an already proved step and can be used directly later, state it explicitly in `lemma_statements.md` and reflect it in `theorem_map.md`; do not leave it only implicit in `proof_notes.md`.

## Edit Style
- Prefer small, reviewable diffs.
- Before a large theory edit, briefly explain the plan.
- After each edit, always report:
  - what changed
  - what was checked
  - what remains doubtful
  - what still needs manual review
