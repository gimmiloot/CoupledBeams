---
name: dorofeev-scientific-writing
description: Draft, locally revise, or audit Russian-language mechanics paper text in the approved CoupledBeams style. Use for scientific prose, LaTeX sections, abstracts, derivations, results, figure or table discussions, and conclusions when terminology, formulas, notation, values, citations, and claim strength must remain unchanged. Do not use for changing scientific models, deriving new results, translating into another language, or journal-layout work alone.
---

# Dorofeev Scientific Writing

Write Russian mechanics papers as clear sequences of supported statements. Reproduce the stable working method of the four-paper corpus, not isolated wording or its editorial defects. The result should need only local author revision.

## Preserve authority and meaning

Use sources in this order:

1. the user's instruction for the current task;
2. the latest approved version of the current manuscript;
3. the current project's terminology, notation, theory, and assumptions;
4. the primary scientific sources behind the model;
5. the phrase families in this skill;
6. a neutral descriptive formulation.

Treat formulas, symbols, numerical values, parameter ranges, citations, labels, scientific assumptions, and claim strength as frozen during stylistic work. Do not invent results, mechanisms, criteria, terms, or evidence. If support is insufficient, keep the wording neutral or state that the claim is not established.

## Choose the mode

### Draft

Write a new requested fragment.

1. Read the current manuscript and relevant project documents.
2. Identify the function of each requested paragraph, its evidence, and approved terms.
3. Make a short internal paragraph plan; do not show it unless asked.
4. Write only the requested fragment.

### Revise

Revise the smallest requested span. Preserve scientific meaning, formulas, notation, values, citations, LaTeX commands, and scope. Split overloaded sentences, remove empty evaluation, and replace invented wording with an already approved term. Do not rewrite adjacent correct text for variety.

### Audit

Do not rewrite unless asked. Report the location, problem type, reason, and smallest correction direction. Separate style defects from unsupported scientific claims.

## Read only the references needed

- For any Draft, Revise, or Audit task, read [terminology policy](references/terminology_policy.md), [simplicity rules](references/simplicity_rules.md), and [claim discipline](references/claims_and_caution.md).
- For section planning, read [structure templates](references/structure_templates.md).
- Before choosing recurring transitions, read [approved phrase families](references/approved_phrase_families.md). Reuse a precise term instead of creating synonyms.
- For prose about figures or tables, read [figures and tables](references/figures_and_tables.md).
- For Audit or final cleanup, read [forbidden patterns](references/forbidden_patterns.md).
- Read [corpus map](references/corpus_map.md), [evidence matrix](references/corpus_evidence.md), and [corpus examples](references/corpus_examples.md) only when maintaining, explaining, or revalidating this skill. The original articles are not runtime dependencies.
- Read [acceptance cases](references/style_acceptance_cases.md) only when testing a change to the skill.

Use `scripts/style_lint.py` as an optional warning-only check for a file or stdin. Treat its findings as prompts for review, not as scientific judgments.

## Writing rules

- Give each paragraph one clear function and each sentence one main thought.
- Prefer ordinary sentences of about 12--25 words. Review sentences longer than 35 words and split them when the logical connection permits.
- Prefer a direct verb to a chain of verbal nouns. Do not make the prose telegraphic.
- Repeat exact scientific terms when needed. Never synonymize them merely for variety.
- Introduce a parameter or symbol immediately before its first use. Move between formulas step by step and explain only new quantities, the formula's role, and the consequence used next.
- Describe a figure or table through the main dependence, difference between cases, supported mechanical meaning, and actual domain of the conclusion. Use numbers as evidence for a material conclusion, not as a substitute for analysis.
- Distinguish an observation in data, a conclusion supported by analysis, and a hypothesis. Limit finite-grid conclusions to the sampled domain.
- Do not call a model exact, complete, true, or an error reference without an explicit basis.
- Do not infer modal tracking, mode transformation, or localization from sorted frequencies alone. Require the evidence stated in the current project, such as continuation, mode shapes, MAC, energy, partial problems, limiting cases, or direct equations.

## Output contract

Follow the requested scale exactly. Return only a replacement block when one is requested, and only an audit when review without rewriting is requested. Preserve LaTeX, references, `\label`, equations, and established notation. Default to Russian and to LaTeX when editing a `.tex` manuscript. Do not prepend a process explanation or provide multiple variants unless asked.

## Final self-check

Before returning scientific prose, correct any failure in this list:

1. Does every paragraph have one clear function?
2. Can any overloaded sentence be divided without losing a necessary relation?
3. Did any term appear that is absent from the current text, project documents, corpus, or primary source?
4. Was an exact term replaced only for verbal variety?
5. Can an abstract noun chain be replaced with a concrete action?
6. Is one conclusion repeated in different words?
7. Is every conclusion limited to the domain actually studied?
8. Was a cause added that the evidence does not establish?
9. Does the prose merely enumerate values visible in a figure or table?
10. Is every retained number needed to support a material conclusion?
11. Is a sorted spectral line being presented as a tracked mode without tracking evidence?
12. Are formulas, references, labels, symbols, and numerical values unchanged?
13. Were any unrequested neighbouring passages changed?
14. Is the first sentence understandable without the hidden drafting plan?
15. Does the result read as a mechanics argument rather than a report about the writing process?

## Invocation examples

- `Use $dorofeev-scientific-writing in Draft mode to write the Russian discussion of Figure 4 from the supplied verified results.`
- `Use $dorofeev-scientific-writing in Revise mode to simplify only this LaTeX paragraph without changing formulas or values.`
- `Use $dorofeev-scientific-writing in Audit mode to identify unsupported claims and terminology drift in this conclusion.`
