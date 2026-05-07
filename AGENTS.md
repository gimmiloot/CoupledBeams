# AGENTS.md

## Project overview

This repository contains the `CoupledBeams` research project.

Main components:
- theory and derivations: `docs/theory/`
- literature and source notes: `docs/literature/`
- project journal: `docs/project_log/journal.md`
- analytic Python code: `src/my_project/analytic/`
- FEM Python code: `src/my_project/fem/`
- runnable scripts: `scripts/`
- tests: `tests/`

## Source of truth

When there is a conflict, use this priority:
1. verified formulas in `docs/theory/equations.tex`
2. explicitly recorded assumptions in `docs/theory/assumptions.md`
3. source notes in `docs/literature/source_index.md`
4. code in `src/my_project/`

Do not silently "fix" mathematics in code if the theory files say otherwise.
Do not silently "fix" theory if the code differs.
Always report the mismatch first.

## Notation policy

Base notation should follow:
- `docs/literature/pdf/Статья-Дорофеев-2025.pdf`

Known warning:
- `docs/literature/pdf/2003JSVb.pdf` contains a known sign issue in the determinant-like matrix formula discussed in the source notes.
- Do not copy sign patterns from that source blindly.
- Always compare against the local verified theory and code.

Do not introduce new notation unless necessary.
If notation must change, update:
- `docs/theory/equations.tex`
- `docs/theory/assumptions.md`
- relevant code comments/docstrings
- `README.md` if the change affects project usage

## Branch identity rule

For analytic branches, `branch_id` is defined by continuation from the base point
`beta = 0`, `mu = 0` for each `epsilon` independently.
`current_sorted_index` is only the current position of that branch in the sorted
spectrum and must not be used as branch identity. Analytic frequency plots and
analytic shape plots must use `scripts/lib/analytic_branch_tracking.py` for
branch selection. Tracking CSV/JSON files are optional debug artifacts only and
must not be treated as the source of truth.
Branch tracking results used in figures must be regression-tested when they
resolve previous inconsistencies. Low-MAC assignments are not canonical unless
they are explicitly accepted in a diagnostic run.
Adaptive refinement and failure CSVs are diagnostic machinery only; they do not
change the definition of branch identity.

## Hard rules for mathematical edits

Treat verified formulas, determinants, matrix entries, sign conventions, and unknown ordering as frozen unless the user explicitly asks for a mathematical revision.

Do not:
- change determinant entries
- change signs or coefficients
- reorder unknowns
- rewrite formulas into "equivalent-looking" forms
- simplify symbolic expressions in ways that hide sign structure

If a formula looks suspicious:
1. do not patch it immediately
2. identify the exact location
3. compare with code and literature
4. state why it is suspicious
5. propose the smallest possible correction only after evidence is given

## Symbolic and formal checking

For nontrivial algebraic transitions:
- prefer explicit symbolic verification over intuition
- check matrix dimensions and unknown ordering
- check limiting cases (`mu = 0`, `beta = 0`, symmetric limits, fixed-fixed limits, etc.)
- when feasible, verify reductions with symbolic tools or reproducible scripts

If a Lean workflow is available or requested, use Lean to verify delicate symbolic transitions, proposition statements, or reductions.
If Lean is not available, state that clearly and fall back to symbolic or numerical verification.

Never claim a symbolic reduction is correct without either:
- direct derivation,
- symbolic verification,
- or numerical consistency checks in representative limits.

## Physical meaning checks

Every new hypothesis, proposition, or simplification must be checked for physical meaning.

When introducing a new assumption, always ask:
- what physical regime does it represent?
- is it dimensionally consistent?
- does it preserve the intended boundary conditions?
- what happens in the symmetric limit?
- does it contradict the verified code or reference theory?
- is it only a computational convenience, or a physical statement?

Do not accept a mathematically convenient assumption if its physical interpretation is unclear or implausible.
Record all new assumptions in `docs/theory/assumptions.md`.

## Theory/code consistency workflow

When touching theory or analytic code:
1. identify the corresponding block in `docs/theory/equations.tex`
2. identify the corresponding block in `src/my_project/analytic/`
3. compare formulas, parameter definitions, and unknown ordering
4. check at least one limiting case
5. report consistency or mismatch explicitly

If code is refactored:
- preserve formulas exactly
- preserve numerical meaning
- preserve output unless the user asked otherwise
- confirm that the determinant and root structure are unchanged

## FEM workflow

Treat newly added FEM programs as baseline implementations unless the user explicitly asks for model changes.

For FEM code:
- first make sure the baseline version runs
- do not change the discretization, constitutive model, or boundary conditions without explicit instruction
- separate technical cleanup from mathematical/model changes
- document required inputs, outputs, and run steps

## Documentation rules

Update documentation when relevant:
- `docs/theory/assumptions.md` for new assumptions
- `docs/literature/source_index.md` for source-specific warnings or corrections
- `CHANGELOG.md` after every meaningful project change, including theory consistency audits, assumption or notation policy changes, code refactors, added baseline scripts, changed run instructions, and added tests
- `README.md` whenever it becomes relevant to user-visible changes, especially project structure, run commands, dependencies, baseline analytic or FEM workflows, source-of-truth policy, or important literature warnings

Do not update `README.md` mechanically if nothing user-visible changed.
Keep changelog entries factual, concise, and tied to actual edits.

At the end of any substantial Codex task, before finishing:
- check whether `CHANGELOG.md` needs an entry
- check whether `README.md` has become stale
- report explicitly whether `README.md` and `CHANGELOG.md` were updated or intentionally left unchanged

## Change discipline

Prefer small, reviewable diffs.
Before major edits, state what will be changed.
After edits, report:
- what changed
- what was verified
- what remains uncertain
- what needs manual review

## Testing and verification

For code changes, prefer lightweight verification first:
- import test
- smoke run
- determinant spot checks
- root consistency checks
- one or two known benchmark limits

If a test cannot be run, say so explicitly.
Never claim "works" without stating what was actually checked.

## External article writing

For external-facing texts, do not use internal project-note style.

For article work in `paper_dorofeev_style/`, follow:
- `skills/external_writing_general.md`
- `skills/paper_dorofeev_style_article_writing.md`

Keep `AGENTS.md` short; detailed writing rules belong in `skills/`.
