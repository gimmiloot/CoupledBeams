# Structure templates

Use these as functional sequences, not mandatory heading sets. Preserve the current manuscript's structure when revising it.

## Abstract

Preferred sequence:

1. object of study;
2. varied parameters or studied effect;
3. model;
4. solution or comparison method;
5. principal supported results.

Omit a literature review, generic significance claims, formulas, long parameter lists, and unsupported applications.

## Introduction

Preferred sequence:

1. immediate scientific field;
2. related studies grouped by function;
3. what those studies establish;
4. the question that remains or the aspect isolated here;
5. purpose of the paper;
6. parameters, models, or effects considered.

Do not build a catalogue of disconnected sentences beginning `В работе [N]`. Do not broaden the review beyond the purpose of the paper.

## Problem statement

Use the applicable elements in this order:

1. object and geometry;
2. supports and coupling;
3. coordinate systems;
4. unknown fields;
5. geometric parameters;
6. material properties;
7. mechanical model;
8. governing equations;
9. boundary conditions;
10. coupling conditions;
11. characteristic equation or other final mathematical problem.

Do not expose program matrices, storage layouts, solver settings, or code-oriented names unless they are independently part of the scientific method.

## Transition through a derivation

Use a sequence of small logical steps:

1. define the new quantity;
2. state the condition or substitution;
3. write the resulting equation or system;
4. state its role;
5. give only the consequence needed next.

Useful transitions include `Обозначим...`, `Тогда...`, `Из условия ... следует...`, `После подстановки получим...`, and `Корни этого уравнения определяют...`.

After a formula, explain new symbols and the formula's function. Do not repeat the full formula in words.

## Results paragraph

Choose only the elements needed for one question:

1. what is shown;
2. conditions of the comparison;
3. meaning of lines, points, or cases;
4. main dependence;
5. difference between cases;
6. mechanical interpretation, if supported;
7. domain of the conclusion.

One paragraph should not simultaneously introduce the model, enumerate every plotted value, propose a mechanism, and conclude the entire section.

## Conclusion

Preferred sequence:

1. object and problem studied;
2. model or comparison used;
3. main results;
4. parameter domain and limitations;
5. at most one justified next question.

Do not summarise the paper section by section, advertise significance, or introduce a result absent from the body.

## Local replacement

When the user requests one sentence or one LaTeX block, preserve the surrounding paragraph's function and return only the replacement. Do not add a new subsection or repair neighbouring prose without permission.
