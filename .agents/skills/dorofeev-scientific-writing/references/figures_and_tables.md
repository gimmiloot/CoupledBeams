# Figures and tables

## Core principle

A figure or table already contains numerical data. The prose should explain the main dependence, the difference between cases, and the supported mechanical meaning instead of repeating values visible to the reader.

## Do not

- list graph coordinates in sequence;
- retell a table row by row or column by column;
- replace analysis with repeated constructions `изменилось с ... до ...`;
- repeat the caption;
- give several numbers without stating the dependence they support;
- call a visually apparent change a mechanical mechanism;
- treat lines as the same mode when mode tracking was not performed.

## Identify what matters

Depending on the actual evidence, discuss:

- monotonic or non-monotonic behaviour;
- convergence or separation of curves;
- a maximum, minimum, or change in ordering;
- a difference between models;
- strengthening or weakening of an effect;
- compliance with or violation of a criterion;
- stability across parameter cases;
- an exception to the main dependence.

## Match the claim to the evidence

Do not treat diagnostic quantities as interchangeable. The minimum relevant data depend on the claim:

- branch continuation -> mode shapes along a declared path and a project-approved matching rule; use modal subspaces when individual matching is ambiguous;
- mode composition -> component fields, energy fractions, partial-mode projections, coupling coefficients, or another explicitly defined composition measure;
- localization -> spatial distributions of displacement or energy between or along structural elements;
- mutual transformation -> tracked forms or modal subspaces together with evidence that their character changes;
- curve convergence or separation -> sorted frequencies are sufficient only for this spectral observation.

Partial problems may support interpretation of composition and interaction. Equations, limiting cases, and additional calculations may support a causal interpretation. None of these automatically replaces branch tracking for identity or spatial evidence for localization.

## Use a number only when it does work

Retain a value when it defines a threshold, applicability boundary, global extremum, verification accuracy, method comparison, count of passing cases, or another principal quantitative result used by the conclusion. Usually one decisive number is stronger than a catalogue.

Do not calculate a new difference, rate, or other derived quantity during writing unless the user requests it or supplies it as a verified result. When the trend and comparison answer the question, omit endpoint values that do not establish a necessary threshold, extremum, or validation result.

## Paragraph sequence

1. Identify the figure or table and conditions.
2. Explain lines, points, or cases if not already clear.
3. State the main dependence.
4. Compare the relevant cases.
5. Add a supported mechanism or explicitly keep it as an observation.
6. Limit the conclusion to the displayed domain.

## Example transformation

Avoid:

`Первая частота изменилась с 3,1 до 3,8, вторая -- с 4,2 до 4,9, а третья -- с 5,0 до 5,3.`

Prefer, if the data support it:

`При увеличении параметра все три частоты возрастают, однако низшие частоты изменяются заметнее. Следовательно, рассматриваемая часть спектра не смещается равномерно.`

The second sentence is admissible only if non-uniform change is the intended supported conclusion.
