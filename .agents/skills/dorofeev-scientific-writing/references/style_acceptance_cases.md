# Style acceptance cases

Use these cases to forward-test the skill after material changes. Judge preservation of meaning and scope, not exact wording.

## 1. Direct abstract opening

Prompt: draft the first two sentences of an abstract from a supplied object, model, and verified result.

Pass conditions: the object appears immediately; there is no general motivation, advertising, or invented application.

## 2. Literature-to-purpose transition

Prompt: connect a grouped literature paragraph to a supplied unresolved question and purpose.

Pass conditions: the gap follows from the supplied source roles; the purpose is concrete and not exaggerated; no source is claimed to omit something without evidence.

## 3. Parameter introduction

Prompt: introduce a supplied dimensionless parameter already named in the manuscript.

Pass conditions: the approved name is repeated; definition precedes use; no synonym or new physical interpretation appears.

## 4. Characteristic equation transition

Prompt: connect supplied boundary and coupling conditions to an unchanged determinant equation.

Pass conditions: steps are sequential; the role of non-trivial solvability is clear; no program matrix or altered sign, coefficient, or unknown order appears.

## 5. Overloaded sentence

Input: a 50--70-word sentence containing the task, parameter list, graph observation, explanation, and conclusion.

Pass conditions: two or three simpler sentences; all original facts remain; no new cause or qualification is introduced.

## 6. Invented term

Input: `ориентационно-индуцированная спектральная реорганизация`.

Pass conditions: replace it with an approved current term, or a plain expression such as `изменение спектра при изменении угла ориентации армирующих волокон`; do not create a second neologism.

## 7. Number-by-number graph description

Input: a paragraph with several constructions `изменилась с ... до ...` for one figure.

Pass conditions: state the main trend and relevant difference between curves; retain no more than one decisive value unless more are required; add a mechanism only when evidence is supplied.

## 8. Table retelling

Input: a row-by-row paraphrase of a table containing an applicability threshold and one exception.

Pass conditions: state the general dependence, threshold outcome, and exception or maximum; do not repeat every cell.

## 9. Universal claim from a finite grid

Input: `Эта зависимость выполняется для всех параметров`, supported only by a stated finite sample.

Pass conditions: restrict the conclusion to the sampled values, region, or grid.

## 10. Local replacement

Prompt: replace one specified sentence in a correct paragraph.

Pass conditions: output only the replacement; do not rewrite the paragraph, add commentary, or modify adjacent claims.

## 11. LaTeX preservation

Input: a paragraph containing inline and display mathematics, citation commands, non-breaking spaces, and escaped symbols.

Pass conditions: formulas and commands are byte-for-byte unchanged unless the user explicitly authorises formatting changes; prose alone is revised.

## 12. References and labels

Input: a LaTeX block containing `\cite`, `\ref`, `\eqref`, and `\label`.

Pass conditions: every key and label is preserved; citations are not added, removed, or reassigned without a content request.

## 13. Exact term repetition

Input: a paragraph that repeats `собственные частоты` three times accurately.

Pass conditions: repetition may be reduced by restructuring, but the term is not replaced with `спектральные показатели` or another decorative synonym.

## 14. Sorted frequencies are not tracked modes

Input: ordered frequency curves without continuation, mode shapes, MAC, or equivalent evidence.

Pass conditions: do not call them tracked branches, mode exchange, mode transformation, or localization; retain an observation about ordered frequencies only.

## 15. Observation versus mechanism

Input: a visible non-monotonic curve and an unsupported proposed mechanical cause.

Pass conditions: separate the observed non-monotonicity from the unconfirmed interpretation, or omit the cause.

## 16. Audit-only request

Prompt: audit a conclusion but do not rewrite it.

Pass conditions: each item gives location, problem type, reason, and minimal correction direction; no replacement section is produced.

## 17. Approved source term exception

Input: a phrase listed as disfavoured by this skill but explicitly defined in the current primary source and approved manuscript.

Pass conditions: keep the term and note no style defect merely because the linter flags it; a local allowlist may suppress the exact warning.

## 18. Current manuscript overrides the corpus

Input: the current manuscript uses `угол ориентации армирующих волокон`, while an older corpus article uses a different expression.

Pass conditions: preserve the current term throughout and do not blend the alternatives.
