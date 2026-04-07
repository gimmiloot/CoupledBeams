# AGENTS.md

## Project purpose
Этот репозиторий посвящен второй научной задаче по механике/устойчивости.

## Source of truth
- Теория: `docs/theory/`
- Литература: `docs/literature/`
- Журнал проекта: `docs/project_log/journal.md`

## Code layout
- Аналитические программы: `src/my_project/analytic/`
- FEM-код на Python: `src/my_project/fem/`
- Запускные скрипты: `scripts/`

## Working rules
- Не менять обозначения без согласования с `docs/theory/`.
- При новых допущениях обновлять `docs/theory/assumptions.md`.
- Не складывать временные результаты в `src/`.
- Новые вычислительные результаты сохранять в `results/`.

## Output style
- Предпочитать минимальные, аккуратные изменения.
- Перед крупной правкой кратко описывать, что будет изменено.