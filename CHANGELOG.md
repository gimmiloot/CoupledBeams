# CHANGELOG

## 2026-05-02

- Updated `src/my_project/analytic/FreqMuNet.py` as the standalone fixed-`beta`
  `Lambda(mu)` plotting entrypoint: added CLI controls for `--beta`,
  `--epsilon`, `--num-modes`, dashed CS reference-line counts, `mu` range, and
  output path, without changing determinant assembly, shared root-finding, or
  FEM/data files.

## 2026-04-22

- Integrated four new veering/localization literature sources into the project
  literature system: added notes, bibliography entries, and source-index blocks
  for Manconi--Mace 2017, Ehrhardt et al. 2018, Lacarbonara et al. 2005, and
  Fontanela et al. 2021; updated the veering assessment/terminology only to
  clarify mechanism-vs-geometry priority and the weak-coupling versus
  strong-coupling distinction, with no solver-code, script, test, data, figure,
  or manuscript edits.

## 2026-04-21

- Added `nair_1973_quasi_degeneracies` to the project literature system and
  recorded its use as a mechanism-level source for quasi-degeneracy / rapid
  modal-reorganization vocabulary; updated veering notes with a conservative
  `bending_desc_04` tracked-pair search from existing `beta = 15 deg`,
  `r = 5 mm` data only, without new solves, code changes, or manuscript edits.

- Синхронизирована проектная теория с текущим текстом статьи в `paper_dorofeev_style/`: в `docs/theory/main_note.md` и `docs/theory/assumptions.md` зафиксированы смысл параметра `\mu`, нормировка `\Lambda`, эталонная роль `\mu=0`, пунктирные CS reference-линии и необходимость интерпретировать `\Lambda(\mu)` вместе с формами колебаний без изменения формул, determinant, solver logic, FEM baseline или текста статьи.

## 2026-04-20

- Переписано введение статьи в `paper_dorofeev_style/sections/introduction.tex`: убраны неподтверждённые обобщения, встроены реальные ссылки по месту, а работы `perkins1986`, `pierre1988` и `liu2002` теперь используются как источники по общим спектральным механизмам `veering`/`mode localization`, а не как прямые аналоги текущей геометрии.
- Обновлена локальная библиография статьи в `paper_dorofeev_style/bib/references.bib`: добавлены записи `liu2002derivatives` и `bauer2025coupledrods`, причём ссылка на предыдущую нашу статью теперь привязана к DOI `10.24412/0136-4545-2025-3-73-81`.

## 2026-04-19

- Rebuilt article Figure 2 in `paper_dorofeev_style/generate_article_spectral_figures.py` around branch identity instead of the smoothed sorted-spectrum CSV: the figure now seeds working branch numbers from the `beta = 0` bending descendants, tracks FEM descendants across `0 <= beta <= 90 deg` by MAC+frequency continuity, matches analytic roots to those tracked descendants at each `beta`, and therefore preserves branch colors/numbers through real crossings without changing the determinant, shared solver logic, or baseline FEM.
- Updated `paper_dorofeev_style/sections/baseline_spectral_picture_mu0.tex` so the Figure 2 caption explicitly states that the right-edge numbers denote the working branch numbering and that color/number stay attached to a branch even when its current place in the spectrum changes.
- Reworked the article mode-shape presentation layer for Figures 4 and 5: `paper_dorofeev_style/generate_targeted_analysis_figures.py` now saves silent per-panel PNGs with only the undeformed geometry and the two deformed arm curves, while `paper_dorofeev_style/manuscript.tex` and `paper_dorofeev_style/sections/mode_shapes_and_frequency_shift.tex` assemble those panels through `subfigure` blocks and move the panel labels `а)`–`е)` / `а)`–`в)` into LaTeX without changing the tracked states, branch logic, or data.

- Fixed the article mode-shape sign presentation in `paper_dorofeev_style/generate_targeted_analysis_figures.py`: within each figure, each branch family now uses a reference panel plus scalar-product sign alignment, with the reference orientation anchored by the short-arm deflection direction, so panel-to-panel sign flips no longer create artificial visual "rearrangements" of the same tracked mode.

## 2026-04-18

- Added `scripts/analyze_target_descendants_beta15_r5.py` for targeted analysis-only tracking of `bending_desc_01`, `bending_desc_02`, and `bending_desc_04` at `beta = 15 deg`, `r = 0.005`, and `0 <= mu <= 0.9`: the script reuses the established FEM continuation and MAC-based descendant tracking, counts local half-wave order on the left and right arms from reconstructed beam-shape profiles, compares each tracked state against single-rod CS reference families, and saves compact CSV/PNG outputs in `results/` without changing the determinant, shared solver logic, or baseline FEM.
- Updated `paper_dorofeev_style/sections/analytic_model.tex` and the `Lambda(mu)` figure caption to document the compact coupled-system characteristic determinant and to identify the dashed reference curves as single-rod clamped-supported families with $l_1=l(1-\mu)$ and $l_2=l(1+\mu)$, without changing the verified mathematics or determinant code.

## 2026-04-15

- Added `scripts/plot_flat_mu_bending_desc_01_mu0_0p1_0p2_ru.py` to generate the clean Russian `2x3` mode-shape figure for the `beta = 15 deg` descendant branch `bending_desc_01` at `r = 0.005, 0.02` and `mu = 0, 0.1, 0.2`, reusing the established FEM MAC-tracking plotting workflow without changing the determinant, shared solver logic, or baseline FEM.
- Added `scripts/plot_flat_mu_bending_desc_02_mu0_0p1_0p2_ru.py` to generate the clean Russian `2x3` mode-shape figure for the `beta = 15 deg` descendant branch `bending_desc_02` at `r = 0.005, 0.02` and `mu = 0, 0.1, 0.2`, reusing the established FEM MAC-tracking plotting workflow without changing the determinant, shared solver logic, or baseline FEM.

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
