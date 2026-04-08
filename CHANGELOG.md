# CHANGELOG

## 2026-04-07

- Проведена проверка согласованности формул между `docs/theory/equations.tex` и аналитическими программами в `src/my_project/analytic/`.
- Зафиксировано замечание по `docs/literature/pdf/2003JSVb.pdf`: determinant-like запись матрицы `T` в формуле `(39)` требует ручной сверки знаков с условиями непрерывности и не должна переноситься в код без проверки.
- Выполнен структурный рефакторинг `src/my_project/analytic/FreqFromAngle.py` и `src/my_project/analytic/FreqFromMu.py` без изменения формул, determinant, порядка неизвестных и математических соглашений.
- Общий математический слой вынесен в `src/my_project/analytic/formulas.py` и `src/my_project/analytic/solvers.py`.
- Добавлен smoke test `tests/test_analytic_smoke.py` для проверки контрольных значений после рефакторинга.
- Добавлена baseline-программа `src/my_project/analytic/FreqMuNet.py` для `mu`-графика в переменной `Lambda` с дополнительными reference-кривыми одиночного стержня CS.
- `src/my_project/analytic/FreqFromMu.py` и `src/my_project/analytic/FreqMuNet.py` переведены на общий `mu`-sweep слой в `src/my_project/analytic/solvers.py` без изменения determinant, формул и численного смысла; для `FreqMuNet.py` сохранён baseline greedy branch tracking, для `FreqFromMu.py` — текущий shared tracking mode.
- `README.md` обновлён под новый analytic workflow и запуск `FreqMuNet.py`.
- Добавлен сценарий `scripts/plot_freq_mu_vs_fem.py` для наложения аналитических ветвей `FreqFromMu` на tracked FEM-частоты из `results/fem_spectrum.csv` и сохранения comparison PNG в `results/`.
- Обновлено оформление `scripts/plot_freq_mu_vs_fem.py`: FEM-маркеры теперь рисуются цветом соответствующей аналитической моды, при этом легенда comparison-графика остаётся компактной.
