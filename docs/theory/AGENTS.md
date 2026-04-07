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

## Edit Style
- Prefer small, reviewable diffs.
- Before a large theory edit, briefly explain the plan.
- After each edit, always report:
  - what changed
  - what was checked
  - what remains doubtful
  - what still needs manual review
