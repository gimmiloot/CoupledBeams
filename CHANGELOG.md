# CHANGELOG

## 2026-05-06

- Added adaptive refinement and failure diagnostics to `scripts/lib/analytic_branch_tracking.py`: low-MAC transitions now bisect beta/mu steps before failing, missed sign-change roots can be recovered locally by minimizing the determinant matrix smallest singular value for tracking candidates, failure reports include per-candidate MAC/frequency/cost CSVs under `results/debug/`, and the desc05 plotter exposes refinement controls, with no determinant, FEM baseline, formula, solver, physical-model, or eigenfrequency changes.
- Improved low-MAC diagnostics for analytic branch tracking and `scripts/analysis/plot_desc05_full_shapes_beta15_eps_sweep.py`: beta/mu continuation steps are explicit user parameters, fail-fast messages now include branch id, beta, mu, MAC, current sorted index, Lambda, and a refinement recommendation, exploratory `--allow-low-mac` runs are labeled in stdout, and tracking debug CSVs include an explicit `Lambda` column under `results/debug/`, with no determinant, FEM baseline, formula, solver, physical-model, or eigenfrequency changes.
- Fixed `scripts/analysis/plot_desc05_full_shapes_beta15_eps_sweep.py` so the historical known-contradiction diagnostic is opt-in behind `--check-known-contradiction` and skips gracefully when the current tracking result does not contain `beta = 15`, `epsilon = 0.0025`, `mu = 0.8`, with no determinant, FEM baseline, formula, solver, physical-model, or eigenfrequency changes.
- Added `tests/test_analytic_branch_tracking_regression.py` for the `bending_desc_05`, `beta = 15 deg`, `epsilon = 0.0025` tracking inconsistency, made low-MAC analytic assignments fail fast unless `allow_low_mac` / `--allow-low-mac` is explicitly requested, added branch-Lambda/current-index consistency assertions, and documented that low-MAC results are not canonical, with no determinant, FEM baseline, formula, solver, physical-model, or eigenfrequency changes.
- Added `scripts/run/run_lambda_mu_fixed_beta_analytic.py` as a convenient fixed-`beta` analytic `Lambda(mu)` plotting entrypoint for arbitrary `--num-modes`, preserving the existing `FreqMuNet.py` visual style, using shared canonical analytic branch tracking, and writing a PNG plus summary CSV such as `results/lambda_mu_beta15_eps0p0025_modes7_dash7.*`, with no determinant, FEM baseline, formula, solver, physical-model, or eigenfrequency changes.
- Introduced `scripts/lib/analytic_branch_tracking.py` as the shared source-of-truth helper for analytic branch identity from `beta = 0`, `mu = 0` for each `epsilon`, separated `branch_id` from `current_sorted_index`, made tracking CSVs optional debug artifacts only, synchronized the desc05 analytic shape sweep and `FreqMuNet.py` Lambda(mu) branch selection through the same helper, and documented the rule in `AGENTS.md` and `scripts/README.md`, with no determinant, FEM baseline, formula, solver, FEM tracking, or eigenfrequency changes.
- Updated `scripts/analysis/plot_desc05_full_shapes_beta15_eps_sweep.py` so the fifth-descendant epsilon-sweep figures use analytic-only branch tracking from `beta = 0`, `mu = 0` for each `epsilon` independently, select roots by sampled-shape MAC instead of FEM Lambda, label `current_sorted_index` instead of branch identity as a root number, and save full tracking-path CSV diagnostics only behind `--save-tracking-debug`, with no determinant, FEM baseline, formula, solver, FEM tracking, or eigenfrequency changes.
- Updated `scripts/analysis/validate_article_shape_cases_beta15.py` to use explicit corrected direct FEM element-stiffness energy fields, add transverse-only and axial-only shape metrics, split frequency/transverse/full-component/energy flags, and write `results/article_shape_validation_beta15_report.md`, with no determinant, FEM baseline, formula, solver, root, tracking, or eigenfrequency changes.
- Fixed FEM arm-wise energy diagnostics by making `arm_energy_diagnostics(...)` use direct local element-stiffness energy blocks for both arms, preserving the legacy nodal-gradient implementation as `legacy_arm_energy_diagnostics(...)` for reference and adding a regression smoke test for `bending_desc_05`, `beta = 30 deg`, `mu = 0`, `epsilon = 0.0025`, with no FEM baseline, determinant, formula, solver, root, tracking, or eigenfrequency changes.
- Added `scripts/analysis/validate_article_shape_cases_beta15.py` for batch analytic-vs-FEM component and energy validation of the article `beta = 15 deg`, `epsilon = 0.0025` small-angle cases (`bending_desc_01`, `bending_desc_02`, `mu = 0, 0.1, 0.2` by default), writing `results/article_shape_validation_beta15_summary.csv` plus per-case overlay/diagnostic files with no determinant, FEM, formula, solver, root, tracking, or eigenfrequency changes.
- Added `scripts/run/run_analytic_coupled_rods_mode_shape_ru.py` and `scripts/lib/analytic_coupled_rods_shapes.py` for analytic mode-shape reconstruction directly from the determinant nullspace, including Russian-labeled full/transverse/components plots, samples CSV, diagnostics CSV, endpoint residuals, and analytic arm-energy shares, with no determinant, FEM, formula, solver, root, tracking, or eigenfrequency changes.
- Added a direct FEM element-stiffness energy check to `scripts/analysis/compare_analytic_fem_tracked_descendant_shape.py`, writing `*_fem_direct_energy_check.csv` behind `--check-fem-direct-energy` and comparing `arm_energy_diagnostics(...)` against energies recomputed from `fem.elem_K(Le)` local axial/bending blocks, with no FEM baseline or `python_fem.py` changes.
- Added analytic-vs-FEM arm-wise axial/bending energy comparison diagnostics to `scripts/analysis/compare_analytic_fem_tracked_descendant_shape.py`, writing `*_energy_comparison.csv` behind `--compare-energies` and appending key energy fields to the diagnostics CSV, with no determinant, FEM, formula, solver, root, tracking, or eigenfrequency changes.
- Added endpoint/matrix consistency diagnostics to `scripts/analysis/compare_analytic_fem_tracked_descendant_shape.py`, writing `*_endpoint_consistency.csv` and an optional text summary for analytic null-vector rows, analytic basis columns, external clamps, and FEM joint kinematics, with no determinant, FEM, formula, solver, root, tracking, or eigenfrequency changes.
- Added an orientation/sign-convention scan to `scripts/analysis/compare_analytic_fem_tracked_descendant_shape.py`, writing `*_orientation_scan.csv`, printing the top matching variants, and keeping the main overlay on the current convention unless `--use-best-orientation` is explicitly supplied, with no determinant, FEM, formula, solver, root, or tracking changes.
- Added `scripts/analysis/compare_analytic_fem_tracked_descendant_shape.py` for analytic-vs-FEM local component comparison of tracked descendants; it reconstructs the analytic mode shape from the determinant nullspace and changes no determinant, FEM, formulas, root-finding, or branch tracking.

## 2026-05-05

- Added `--diagnostics-level all` to the single tracked bending descendant shape runner, with local projection reconstruction checks, strain-like local derivative diagnostics, and arm-wise axial/bending diagnostic energies for the selected mode shape, without changing determinant, FEM, formulas, root-finding, branch tracking, or eigenfrequencies.
- Added an editable `USER PARAMETERS` block to the single tracked bending descendant shape runner so no-argument runs use file-local working defaults while CLI arguments remain available as overrides, without changing determinant, FEM, formulas, root-finding, branch tracking, or eigenfrequencies.
- Expanded the single tracked bending descendant shape runner CLI with explicit branch defaults, `--l-total`, `--dpi`, `--figsize`, `--show`, `--normalize auto`, deterministic scale-aware output names, and optional one-row diagnostics CSV export, without changing determinant, FEM, formulas, root-finding, branch tracking, eigenfrequencies, or parameter meanings.
- Added `full`, `transverse`, and `components` plot modes plus `--mode-scale` and normalization controls to the single tracked bending descendant shape runner; the new diagnostics report local axial/transverse amplitudes without changing determinant, FEM, formulas, root-finding, branch tracking, eigenfrequencies, or baseline mode-shape normalization.
- Generalized Russian title-label inference for tracked bending descendants so any `bending_desc_NN` branch id renders as `потомок N-й изгибной ветви`, with a unit-style fallback check and no determinant, FEM, formula, root-finding, tracking, normalization, layout, or color changes.
- Cleaned up the single-case tracked bending descendant shape PNG layout by removing its top legend, shortening the title to the requested beta/mu/epsilon and Russian branch label metadata, and reducing the single-case figure height, without changing tracking, normalization, scale, colors, formulas, determinant, root-finding, or FEM logic.
- Split tracked bending descendant shape plotting into a reusable `scripts/lib/tracked_bending_descendant_shapes.py` core used by both the single-shape runner and the multi-panel builder; multi-panel figures now call shared in-process helpers instead of depending on single-case CLI behavior, with no determinant, FEM, formula, root-finding, or branch-tracking changes.
- Added `scripts/run/run_tracked_bending_descendant_shape_ru.py` as a single-shape user-facing runner for tracked bending descendant mode shapes; it exposes descendant branch number/id, `mu`, `epsilon`, and `beta`, without changing determinant, FEM, formula, root-finding, or branch-tracking logic.
- Introduced `scripts/plot_tracked_bending_descendant_shapes_ru.py` as the parameterized tracked bending-descendant mode-shape plotting entrypoint; converted the old `scripts/plot_flat_mu_bending_desc*.py` scripts into compatibility wrappers that preserve their branch ids and output PNG paths, without changing determinant, FEM, formula, or branch-tracking logic.
- Added a `scripts/README.md` script map, introduced clear user-facing wrappers in `scripts/run/`, and documented analysis/helper/legacy groupings while preserving old root-level command paths; mathematical model, determinant/root logic, and FEM baseline were unchanged.

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
