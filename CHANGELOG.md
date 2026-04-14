# CHANGELOG

## 2026-04-13

- Added `scripts/plot_mu_sweep_beta_fixed_four_radii_compare.py` to generate shared `2x2` fixed-`beta` mu-sweep figures in `Lambda` for `beta = 7.5 deg` and `beta = 15 deg` across radii `r = 0.005, 0.01, 0.015, 0.02`, reusing the established positive-`beta` low-branch comparison workflow together with the same presentation style and muted CS reference lines as the existing `beta = 0` four-radii figure.
- Added `scripts/plot_mu_sweep_beta0_four_radii_compare.py` for a shared `2x2` analytic-vs-FEM `beta = 0` mu-sweep figure in `Lambda` at radii `r = 0.005, 0.01, 0.015, 0.02`, reusing the existing type-aware bending matching and muted CS single-rod dashed reference-line logic and saving one common PNG together with per-radius and combined CSV tables in `results/`.
- Добавлен отдельный сценарий `scripts/compare_beta0_analytic_vs_fem.py` для сравнения determinant-based аналитических частот с baseline FEM при `beta = 0` без изменения формул, determinant или solver logic.
- Новый сценарий использует существующий `results/fem_spectrum.csv`, если он согласуется с текущим baseline FEM при `mu = 0`; в противном случае пересчитывает `beta = 0` sweep из `src/my_project/fem/python_fem.py` в памяти.
- В сценарий добавлены воспроизводимые выгрузки в `results/`: таблица сравнения при `mu = 0`, таблица классификации bending/axial мод по FEM-формам при `mu = 0`, таблицы radius-sensitivity при `mu = 0` и comparison plot для первых `beta = 0` ветвей.
- Сценарий `scripts/compare_beta0_analytic_vs_fem.py` расширен type-aware сопоставлением мод: отдельный `mu`-график для bending-ветвей, отдельные графики первой axial-ветви по `mu` и по `r`, а также CSV-таблицы с `analytic_branch_id`, `fem_mode_id`, `mode_type`, частотами и относительной ошибкой.
- Добавлен сценарий `scripts/compare_beta_positive_type_aware.py` для multibeta type-aware comparison при `beta > 0`: ветви переносятся продолжением от проверенного случая `beta = 0`, FEM-потомки сопоставляются по MAC+frequency continuity, а в `results/` сохраняются multibeta CSV-таблицы и графики bending/axial ветвей и смешения по `fem_axial_fraction`.
- Добавлен сценарий `scripts/analyze_branchwise_fem_spectrum.py` для FEM-only branch-wise анализа полного отслеживаемого спектра на сетке `(beta, mu)`: сохраняются branch-id-таблица с `current_sorted_index`, `axial_fraction`, соседними gap-ами и численными `df/dmu`, `df/dbeta`, а также summary/extrema CSV и графики branch-colored spectrum, axial-fraction heatmaps и карты минимальных gap-ов.
- Добавлен сценарий `scripts/plot_beta_sweep_mu0_compare.py` для презентационного сравнения analytic vs FEM при `mu = 0` на beta-sweep `0..90°` и радиусах `r = 0.005` и `r = 0.01`: в `results/` сохраняются отдельные PNG-графики в переменной `Λ` и CSV-таблицы с `beta`, `branch_id`, `analytic_hz`, `fem_hz`, `relative_error` и matched `fem_mode_id`.
- Сценарий `scripts/plot_beta_sweep_mu0_compare.py` обновлён для сглаженного presentation sweep без не физических high-branch spikes: в проблемных окнах по `beta` добавлены локально уплотнённые точки, analytic acquisition расширен до 16 candidate roots, локальный `scan_step` по `Λ` уменьшен, а для тесных пар используется local continuation поверх существующего determinant/root-finding слоя без изменения формул, shared solver logic или baseline FEM.
- Problem-window smoothing из `scripts/plot_beta_sweep_mu0_compare.py` расширен на радиусы `r = 0.015` и `r = 0.02`, чтобы presentation beta-sweep оставался гладким и без branch-swap artifacts и для дополнительных high-branch close-pair зон.
- Добавлен сценарий `scripts/plot_beta_sweep_mu0_four_radii_compare.py` для общей `2x2` figure analytic vs FEM при `mu = 0` и радиусах `r = 0.005, 0.01, 0.015, 0.02`; в `results/` сохраняются общий PNG, combined CSV по всем 4 радиусам и per-radius CSV для новых случаев `r = 0.015` и `r = 0.02`.

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
