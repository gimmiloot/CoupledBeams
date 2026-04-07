# CHANGELOG

## 2026-04-07

- Проведена проверка согласованности формул между `docs/theory/equations.tex` и аналитическими программами в `src/my_project/analytic/`.
- Зафиксировано замечание по `docs/literature/pdf/2003JSVb.pdf`: determinant-like запись матрицы `T` в формуле `(39)` требует ручной сверки знаков с условиями непрерывности и не должна переноситься в код без проверки.
- Выполнен структурный рефакторинг `src/my_project/analytic/FreqFromAngle.py` и `src/my_project/analytic/FreqFromMu.py` без изменения формул, determinant, порядка неизвестных и математических соглашений.
- Общий математический слой вынесен в `src/my_project/analytic/formulas.py` и `src/my_project/analytic/solvers.py`.
- Добавлен smoke test `tests/test_analytic_smoke.py` для проверки контрольных значений после рефакторинга.
