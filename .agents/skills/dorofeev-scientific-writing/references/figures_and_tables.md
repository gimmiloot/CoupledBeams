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

Name localization, shear influence, rotary inertia, bending--torsion coupling, or another mechanism only when it is supported by mode shapes, energy, MAC, partial problems, a limiting case, an additional calculation, or a direct consequence of equations.

## Use a number only when it does work

Retain a value when it defines a threshold, applicability boundary, global extremum, verification accuracy, method comparison, count of passing cases, or another principal quantitative result used by the conclusion. Usually one decisive number is stronger than a catalogue.

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
