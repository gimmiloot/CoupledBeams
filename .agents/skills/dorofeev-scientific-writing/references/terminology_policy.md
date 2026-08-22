# Terminology policy

## Priority

Use terminology in this order:

1. the user's current instruction;
2. the latest approved manuscript;
3. current project theory, assumptions, source notes, and notation policy;
4. the primary source for the model;
5. this skill's corpus-derived phrases;
6. a neutral descriptive expression.

Older corpus wording never overrides an approved current term.

## Exact repetition is acceptable

Repeat a precise term whenever it denotes the same object. Do not substitute decorative alternatives such as:

- `собственные частоты` → `спектральные показатели`;
- `изменение спектра` → `спектральная реорганизация`;
- `сближение кривых` → `кластеризация ветвей`;
- `угол ориентации армирующих волокон` → `ориентационный параметр`.

## Before using a new term

Check whether the needed wording exists in the current manuscript, project documentation, corpus, or primary source. If not, prefer a plain description. If a technical term is genuinely needed, keep it outside the finished passage until the user approves it or introduce it explicitly with its basis.

## Keep models distinct

Do not mix terminology belonging to Euler--Bernoulli, Rayleigh, Timoshenko, anisotropic, three-dimensional, FEM, or reduced models. Do not describe one model as true, complete, or exact unless the current source of truth explicitly does so.

Distinguish:

- physical frequency, angular frequency, and a dimensionless frequency parameter;
- difference, relative difference, and error when no reference truth is defined;
- sorted spectral index, current sorted position, and continuation-defined `branch_id`;
- observed curve proximity, veering, crossing, localization, and tracked mode exchange.

## Branch and mode language

In CoupledBeams, branch identity follows continuation from the base point as specified by project rules. A sorted frequency line is not a tracked branch. Do not claim mode transformation, exchange, or localization from ordered frequencies alone.

## Source terms

Do not replace a term taken from a primary source with a more elegant synonym. When source terminology conflicts with the current project's approved terminology, report the mismatch rather than blending the terms silently.
