# Forbidden and disfavoured patterns

Avoid these expressions when they merely inflate a simple statement:

- `в рамках настоящего комплексного исследования`;
- `следует особо подчеркнуть`;
- `представляется целесообразным отметить`;
- `всесторонне исследовано`;
- `уникальный подход`;
- `инновационный метод`;
- `значительный научный вклад`;
- `открывает широкие перспективы`;
- `обеспечивает глубокое понимание`;
- `результаты однозначно доказывают`;
- `комплексная спектральная картина`;
- `многофакторная параметрическая обусловленность`;
- `спектральная реорганизация`;
- `модально-спектральная перестройка`;
- `ориентационно-индуцированная связанность`;
- `геометрически обусловленная модальная трансформация`;
- `параметрически индуцированный эффект`;
- `мера продольной геометрической асимметрии`.

Also remove internal workflow labels from external prose, including `TODO`, `FIXME`, `TBD`, `NEEDS_CHECK`, `proved_here`, `needs_caution`, `not_reached_yet`, `branch-ready`, `ceiling diagnosis`, and publisher notes such as `здесь будет` or `вставить`.

## Apply the rule by meaning

Do not ban a phrase mechanically when it is an exact term in the current manuscript or primary source. The defect is unjustified substitution of inflated wording for an accepted precise expression. Record a justified exception in a local allowlist when using `style_lint.py`.

An allowlist is a UTF-8 text file. Each non-comment line contains a warning code, a tab, and the exact normalised warning subject, for example:

```text
TERM001	поворот материальных осей
```

Pass it with `--allowlist FILE`. Use `*` instead of the code only when the same exact subject is intentionally accepted for every warning type. The allowlist cannot disable an entire rule without naming a subject.

## Common repair directions

- Replace an evaluation with the concrete result.
- Replace a coined noun phrase with a direct dependence: `изменение X при изменении Y`.
- Replace `однозначно доказывает` with the evidence-appropriate `показывает` or `наблюдается`.
- Delete a meta-comment instead of translating it into article prose.
