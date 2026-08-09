> Этот набор рисунков относится к реализации теории моноклинных прямоугольных стержней из главы 2 монографии Ярцева. Он не относится к параллельной статье о сравнении теорий Эйлера–Бернулли и Тимошенко для круглых изотропных стержней. Совпадение терминов «Эйлер–Бернулли», «Тимошенко», `beta`, `mu` или `Lambda` не означает совпадения физических моделей, геометрии, материалов или нормировок.

# Рисунки для отчёта научному руководителю

Статус отдельного presentation/report workflow:

```text
Supervisor figure workflow: PASS
Fast beta-sweep optimization: PASS
Extended supervisor figures 5–8: PASS
Small-theta supervisor figures: PASS
```

Workflow создаёт двенадцать воспроизводимых рисунков по реализованной теории главы
2 и не изменяет научные модели, результаты UT-0–UT-3, завершённую 1D phase или
статус 3D-A0. Entry point намеренно отделён от большого unequal-thickness CLI:

```text
python scripts/analysis/anisotropic_rods/plot_yartsev_ch2_supervisor_figures.py
python scripts/analysis/anisotropic_rods/plot_yartsev_ch2_supervisor_figures.py --reuse-data
python scripts/analysis/anisotropic_rods/plot_yartsev_ch2_supervisor_figures.py --validate-fast-solver --jobs 1
python scripts/analysis/anisotropic_rods/plot_yartsev_ch2_supervisor_figures.py --benchmark-fast-solver --jobs 1
python scripts/analysis/anisotropic_rods/plot_yartsev_ch2_supervisor_figures.py --figure new --solver-mode fast --resume --jobs 1
python scripts/analysis/anisotropic_rods/plot_yartsev_ch2_supervisor_figures.py --figure theta-small --solver-mode fast --resume --jobs 1
python scripts/analysis/anisotropic_rods/plot_yartsev_ch2_supervisor_figures.py --figure theta-small --reuse-data --jobs 1
python scripts/analysis/anisotropic_rods/plot_yartsev_ch2_supervisor_figures.py --figure all --reuse-data
```

Результаты сохраняются только в
`results/anisotropic_rods/yartsev_ch2_supervisor_figures/`. Этот каталог не
входит ни в одно workspace статьи.

## Состав рисунков

### Figure 1 — воспроизведение Figure 2.2

`figure_01_yartsev_fig_2_2_reproduction` повторно использует verified evidence
canonical workflow
`scripts/analysis/anisotropic_rods/reproduce_yartsev_fig_2_2.py`:

- `complex_roots.csv`, SHA-256
  `35d69671127e3f74e7305a3cdace1636d23e568b3183340cbacc32f826d8db41`;
- `figure_2_2_digitized_calculated_curves.csv`, SHA-256
  `72aa2e96f58792c325ef07d3ee955ffbe819bc3d644457b1696db437708e43c6`;
- `figure_2_2_digitized_comparison.csv`, SHA-256
  `9d6b1b24c2e18b01b4ef3cff262802142da3d007278d8389b5644c96ad0071ac`.

Источник имеет статус `PASS_WITHIN_GRAPH_RESOLUTION`. Используются
`BookMaterial()`, `book_complex`, canonical geometry и `state_corrected`.
Новая подгонка материала и повторная оцифровка не выполняются. Сохранены
четыре исходные панели: частоты и модальные коэффициенты потерь в книжных
координатах. Линии показывают calculated curves, маркеры — существующие
округлённые значения этих кривых, считанные с исходного рисунка. Отдельные
экспериментальные маркеры монографии не были оцифрованы и не показаны.

### Figure 2 — два вида внешней заделки

`figure_02_clamp_comparison_lambda_vs_beta_book_material` использует только
новую Chapter-2 модель: `state_corrected` Timoshenko bending, generalized
rectangular torsion и неизменённый `J_book(beta)`.

Параметры:

```text
material = BookMaterial() = T-53(VM)-78/PN-609-21M
material_mode = elastic
a1 = a2 = 9.76e-3 m
b1 = b2 = 2.524e-2 m
L1 = L2 = 0.295 m
theta1 = theta2 = 0 deg
mu = 0
beta = 0...90 deg, step 0.5 deg
```

Сплошные линии соответствуют `book_slope_clamp`, то есть
`w=0, w'=0, Phi=0`. Пунктирные линии соответствуют
`timoshenko_section_clamp`, то есть `w=0, psi=0, Phi=0`. EB theory на этом
рисунке отсутствует.

Section-clamp matrix является приватным helper plotting script:

```text
D_section = J_book(beta) blockdiag(T1 B1_section, T2 B2_section)
B_i_section = cantilever_clamp_matrix(
    point_i, "timoshenko_section_clamp", scaled=False
)
```

Он повторно использует production `physical_state_transfer_matrix`,
`equilibrate_matrix`, `cantilever_clamp_matrix` и `joint_matrix_book`; ни один
из этих helpers не изменён.

### Figures 3–4 — Chapter-2 Timoshenko и rectangular EB comparator

Оба рисунка используют одну и ту же прямоугольную постановку:

```text
material = hms_dx_209_material() = HMS/DX-209
material_mode = elastic
a1 = a2 = 0.005 m
b1 = b2 = 0.020 m
b/a = 4
L1 = L2 = 0.400 m
theta1 = theta2 = 0 deg
mu = 0
beta = 0...90 deg, step 0.5 deg
```

На `figure_03_timoshenko_vs_eb_book_slope_lambda_vs_beta` внешние концы имеют
`book_slope_clamp`. На
`figure_04_timoshenko_vs_eb_section_clamp_lambda_vs_beta` новая модель имеет
`timoshenko_section_clamp`. Сплошные линии на обоих рисунках — модель
моноклинного прямоугольного стержня Тимошенко с обобщённым кручением.
Пунктирные линии — прямоугольная ортотропная модель Эйлера–Бернулли с
кручением Сен-Венана. Для EB выполняется `w'=psi`, поэтому её массив ровно один
и тот же на Figures 3–4; сохранённая regression-проверка `array_equal` прошла.

Это сравнение не относится к аналитической модели круглых изотропных стержней.
Круглое сечение, isotropic steel defaults и физические формулы параллельной
статьи не используются.

### Figure 5 — разные длины плеч

`figure_05_timoshenko_vs_eb_unequal_lengths_book_slope` использует
`hms_dx_209_material()`, elastic mode, `L1=0.300 m`, `L2=0.500 m`,
`l=0.400 m`, `mu=0.25`, `a1=a2=0.005 m`, `b1=b2=0.020 m`,
`theta1=theta2=0` и `book_slope_clamp`. Сплошные линии — Chapter-2
`state_corrected` Timoshenko model с generalized rectangular torsion;
пунктирные — существующая rectangular orthotropic EB model с кручением
Сен-Венана.

Готовая подпись: «Сравнение моделей для прямоугольных стержней HMS/DX-209 с
различными длинами плеч: \(L_1=0.3\) м, \(L_2=0.5\) м,
\(a_1=a_2=5\) мм. Сплошные линии — модель моноклинного прямоугольного
стержня Тимошенко с обобщённым кручением; пунктирные линии — прямоугольная
ортотропная модель Эйлера–Бернулли с кручением Сен-Венана. Используется
`book_slope_clamp`. Показаны первые шесть sorted spectral positions при каждом
фиксированном beta; это не отслеженные модальные ветви».

### Figure 6 — разные длины и толщины

`figure_06_timoshenko_vs_eb_unequal_lengths_and_thickness_book_slope`
использует тот же материал, `L1=0.300 m`, `L2=0.500 m`, `l=0.400 m`,
`mu=0.25`, `a1=0.004 m`, `a2=0.006 m`, `b1=b2=0.020 m`,
`theta1=theta2=0` и `book_slope_clamp`. Более длинное второе плечо является
более толстым. Это direct geometry: mass-preserving parameterization,
перестановка толщин и старый thickness parameter не применяются.

Готовая подпись: «Сравнение моделей для прямоугольных стержней HMS/DX-209 с
различными длинами и толщинами плеч: \(L_1=0.3\) м, \(L_2=0.5\) м,
\(a_1=4\) мм, \(a_2=6\) мм. Геометрия задана непосредственно без
массо-сохраняющей параметризации. Сплошные линии — модель Тимошенко с
обобщённым кручением; пунктирные линии — модель Эйлера–Бернулли с кручением
Сен-Венана. Показаны первые шесть sorted spectral positions при каждом
фиксированном beta; это не отслеженные модальные ветви».

### Figure 7 — слабая моноклинность и ортотропное приближение

`figure_07_monoclinic_theta5_vs_orthotropic_eb_approximation` использует
HMS/DX-209, `L1=L2=0.400 m`, `a1=a2=0.005 m`, `b1=b2=0.020 m` и
`book_slope_clamp`. Сплошной spectrum рассчитан полной Chapter-2
`state_corrected` моделью при `theta1=theta2=5 deg`; он включает поворот
эффективных свойств, generalized torsion и `Sbar16` coupling. Пунктирный
orthotropic EB spectrum при `theta=0` не пересчитывается: это точная копия
validated EB array Figure 3, подтверждённая `array_equal`.

Готовая подпись: «Сопоставление полной моноклинной модели главы 2 при
\(\theta_1=\theta_2=5^\circ\) с ортотропным приближением Эйлера–Бернулли
при \(\theta=0^\circ\). Сплошные линии учитывают поворот эффективных свойств
и материальную изгибно-крутильную связь; пунктирные линии соответствуют
приближению, в котором эта связь полностью отсутствует. Рисунок является
диагностикой применимости ортотропного приближения, а не чистым сравнением
только сдвиговой и вращательно-инерционной поправок. Показаны первые шесть
sorted spectral positions при каждом фиксированном beta; это не отслеженные
модальные ветви».

Различие Figure 7 одновременно содержит эффекты shear deformation, rotary
inertia, generalized versus Saint-Venant torsion, rotation of effective
properties и отсутствия материальной bending–torsion coupling в comparator.
Максимум нельзя интерпретировать как чистую EB-vs-Timoshenko truncation error
или как общий запрет EB для всех моноклинных стержней.

### Figure 8 — поворот материальных осей внутри Chapter-2 theory

`figure_08_chapter2_theta15_vs_theta0` использует ту же геометрию и
`book_slope_clamp`. Обе группы кривых являются Chapter-2 `state_corrected`
Timoshenko spectra с generalized torsion. Сплошной spectrum рассчитан при
`theta1=theta2=15 deg`; пунктирный `theta1=theta2=0` spectrum не
пересчитывается и точно совпадает с validated Timoshenko book-slope array
Figure 3.

Готовая подпись: «Влияние поворота материальных осей в полной модели главы 2.
Сплошные линии — \(\theta_1=\theta_2=15^\circ\); пунктирные линии —
\(\theta_1=\theta_2=0^\circ\). Обе группы кривых рассчитаны по одной модели
Тимошенко с обобщённым кручением и `book_slope_clamp`. Показаны первые шесть
sorted spectral positions при каждом фиксированном beta; это не отслеженные
модальные ветви». Различие измеряет объединённое влияние rotated effective
properties и `Sbar16` coupling внутри одной Chapter-2 model; Figure 8 не
является EB-vs-Timoshenko comparison.

### Figures 9–12 — small-theta continuation Figure 7

Figures 9–12 являются одной параметризованной серией, а не четырьмя
независимыми scientific code paths:

| Figure | material-axis orientation |
|---:|---:|
| 9 | `theta1=theta2=1 deg` |
| 10 | `theta1=theta2=2 deg` |
| 11 | `theta1=theta2=3 deg` |
| 12 | `theta1=theta2=4 deg` |

Здесь `theta_i` — угол ориентации главных материальных осей относительно
локальной продольной оси соответствующего стержня. Это не геометрический угол
между стержнями: угол жёсткого сопряжения по-прежнему обозначается `beta`.
Во всех четырёх случаях используются HMS/DX-209, elastic material line,
`a1=a2=0.005 m`, `b1=b2=0.020 m`, `L1=L2=0.400 m`, `mu=0`,
`book_slope_clamp`, `beta=0...90 deg` с шагом `0.5 deg`, семь принятых
положительных корней и первые шесть independently sorted positions.

Solid spectra являются новыми решениями неизменённого Chapter-2 determinant:
`state_corrected`, rotated material properties, generalized rectangular
torsion и `J_book(beta)`. Интерполяция, scaling Figure 3/7, polynomial fitting
и подмена через `E_x(theta)` не применяются. Dashed spectrum не решается
повторно: для каждого рисунка поля Figure-3 rectangular orthotropic EB при
`theta=0` и Saint-Venant torsion копируются из `figure_03_data.csv`; полный
массив `(beta, sorted_position, frequency, Lambda)` проходит `array_equal`.

Все четыре рисунка имеют тот же `ylim`, что сохранённый Figure 7, размер
`7.2 x 4.8 in`, 12 линий (`6 solid + 6 dashed`), общий deterministic color
cycle и не имеют legend. Relative diagnostic называется
`relative_difference_to_theta0_EB_baseline`, поскольку одинаковый индекс
обозначает только одинаковую sorted spectral position.

Для `beta=0` файл `theta_small_beta0_mode_character.csv` содержит по семь
строк для `theta=1,2,3,4 deg`. Fractions bending/shear/torsion вычисляются из
существующих Chapter-2 strain measures и generalized-torsion energy density
по физическим состояниям обоих плеч. Label `bending-like` или `torsion-like`
назначается только при fraction не ниже `0.6`; остальные случаи помечаются
`mixed`. Это диагностическая классификация, не MAC и не mode tracking.

Equal color/index across the two models denotes the same sorted spectral
position, not necessarily the same physical modal descendant.

В рассчитанной beta-zero diagnostic table порядок `torsion-like/bending-like`
для positions 5–6 меняется уже между `theta=0` и `1 deg`; порядок для
positions 3–4 меняется между `theta=3` и `4 deg`. Это локализация по
дискретной сетке углов и energy-character labels, а не доказательство
modal-descendant exchange: MAC и shape tracking в этой серии не выполняются.

## Нормировка Figures 2–4

Для всех кривых применяется только rectangular Chapter-2 normalization:

\[
\Lambda =
\left(\frac{\rho A\omega^2l^4}{E_xI_y}\right)^{1/4}
=
l\left(\frac{\rho A}{E_xI_y}\right)^{1/4}\sqrt{\omega},
\qquad
\omega=2\pi f,
\]

\[
l=\frac{L_1+L_2}{2},
\qquad
A=ab,
\qquad
I_y=\frac{a^3b}{12}.
\]

Здесь `E_x` — elastic modulus при `theta=0`. Одна reference-section
нормировка применяется к обоим clamps, обеим теориям и всем шести показанным
модам. Эквивалентность двух записей формулы покрыта unit test. Старый
thickness-mismatch параметр `eta` отсутствует. Буквенное обозначение
коэффициента потерь на Figure 1 имеет другой, книжный смысл и не участвует в
Figures 2–4.

## Fixed reference normalization Figures 5–12

Все новые curves используют один fixed rectangular reference, в том числе
обе actual-section curves Figure 6 и повернутые материалы Figures 7–12:

\[
a_0=0.005\ \mathrm{m},\qquad b_0=0.020\ \mathrm{m},\qquad
A_0=a_0b_0,\qquad I_{y0}=\frac{a_0^3b_0}{12},
\]

\[
\Lambda_{\mathrm{ref}}=
\left(\frac{\rho A_0\omega^2l^4}{E_{x0}I_{y0}}\right)^{1/4},
\qquad l=\frac{L_1+L_2}{2}=0.400\ \mathrm{m},
\qquad E_{x0}=E_x(\theta=0).
\]

Нормировка не берёт свойства одного фактического плеча Figure 6 и не меняет
`E_x0` при `theta=5` или `15 deg`.

## Fast beta-sweep design и legacy fallback

Diagnostic-only coordinator находится в
`scripts/lib/yartsev_ch2_fast_beta_sweep.py`; физических уравнений, boundary
matrices или нового determinant в нём нет. Plotting entry point передаёт ему
callbacks существующих Chapter-2 builders и root-quality primitives.

- `beta=0` и обязательные anchors `0,15,...,90 deg` выполняют global
  inventory.
- На втором шаге predictor равен предыдущему sorted spectrum, затем применяется
  линейный predictor по двум предыдущим beta. Это только подсказка для окон, а
  не mode tracking.
- Relative neighbor gap `<1e-3` объединяет connected cluster; cluster roots
  восстанавливаются одним локальным inventory и снова сортируются.
- Wrong count, rejected/duplicate roots, predictor crossing, cluster-topology
  change, excessive predictor residual, bracket failure или exception вызывают
  автоматический global fallback. Predictor продолжается от принятого global
  spectrum.
- Legacy `_sweep` и production root finder не переписаны и доступны через
  `--solver-mode legacy`; validated fast path выбирается через
  `--solver-mode fast` и является default.

Transfer matrices `physical_state_transfer_matrix` и
`eb_state_transfer_matrix` кэшируются bounded LRU. Ключ включает model type,
точное `omega` без округления и полный immutable point/material/geometry
fingerprint. Это позволяет identical arms и clamp variants переиспользовать
ровно одинаковый transfer, но не создаёт approximate determinant cache.

Каждая завершённая spectral family атомарно сохраняется в
`fast_family_checkpoints/`. `--resume` принимает checkpoint только при полном
числе строк и совпадении fingerprint, включающего material, geometry, theta,
clamp/model, beta grid, root count, solver version/settings и normalization.

## Oracle-validation и benchmark

До Figures 5–8 быстрый solver был проверен на пяти distinct families
сохранённых Figures 2–4: два clamps Figure 2, два Timoshenko clamps HMS/DX-209
и один shared EB spectrum Figures 3–4. Проверены все
`5 x 181 x 7 = 6335` roots:

```text
maximum relative frequency error = 6.503683351288e-09
maximum relative Lambda error    = 3.251841724285e-09
quality acceptance               = identical for all roots
fast sequential runtime          = 139.714119 s
legacy recorded runtime          = 860.774803 s
speedup                          = 6.160972
mandatory global anchors         = 35
global fallbacks                 = 19
fallback rate                    = 0.0209945
transfer expm, fast              = 310445
transfer expm, legacy estimate   = 5589198
transfer cache hit rate          = 0.775474
```

Hard oracle tolerance `1e-8` и performance target (`<=180 s` либо speedup
`>=5`) пройдены, поэтому сохранён статус `FAST_SOLVER_PASS`. Полная таблица
находится в `fast_solver_validation.csv`, counters и per-family fallback reasons
— в `fast_solver_benchmark.json`. Global fallback является ожидаемой частью
алгоритма, а не замаскированным пропуском.

## Спектр и root-quality gate

При каждом фиксированном `beta` независимо формируются первые семь
положительных корней. Все семь проходят quality gate, первые шесть строятся, а
седьмой служит completeness guard. В CSV и отчёте действует контракт:

```text
The plotted curves are sorted spectral positions 1–6 at every beta,
not tracked modal descendants.
```

Предыдущий отсортированный спектр применяется только как численная подсказка
для узких локальных окон при проверке близких корней. Окончательный спектр
каждого угла вновь объединяется, очищается от дублей и сортируется по частоте.
Это не branch continuation, MAC или shape tracking.

Quality имеет `PASS`, если проходит scaled branch или normalized physical-raw
branch. Неизменённые пороги:

```text
normalized determinant residual <= 1e-8
relative singular residual <= 1e-8
root status does not start with rejected
```

Пропуски не интерполируются. Дополнительная completeness-проверка применяет
существующий canonical root finder с шагами сканирования 2 Hz и 0.5 Hz и
локальными окнами около близких отсортированных положений; уравнения, signs,
thresholds и production root solver не меняются. Максимальные принятые
residuals полного расчёта:

| Figure/model | determinant residual | singular residual |
|---|---:|---:|
| 2, book slope | `2.64e-14` | `2.70e-11` |
| 2, section clamp | `2.52e-14` | `6.95e-13` |
| 3, Timoshenko | `1.92e-14` | `1.34e-10` |
| 3–4, shared EB | `1.97e-14` | `1.42e-10` |
| 4, Timoshenko | `1.99e-14` | `9.13e-13` |
| 5, Timoshenko | `5.05e-15` | `4.65e-13` |
| 5, EB | `4.78e-15` | `4.77e-13` |
| 6, Timoshenko | `1.34e-14` | `5.22e-13` |
| 6, EB | `1.28e-14` | `3.68e-12` |
| 7, Timoshenko theta=5 | `3.17e-15` | `5.82e-13` |
| 7, reused EB theta=0 | `1.97e-14` | `1.42e-10` |
| 8, Timoshenko theta=15 | `1.82e-15` | `7.84e-13` |
| 8, reused Timoshenko theta=0 | `1.92e-14` | `1.34e-10` |
| 9, Timoshenko theta=1 | `1.66e-14` | `5.41e-13` |
| 10, Timoshenko theta=2 | `8.70e-15` | `4.71e-13` |
| 11, Timoshenko theta=3 | `8.23e-15` | `7.94e-13` |
| 12, Timoshenko theta=4 | `6.25e-15` | `7.33e-13` |
| 9–12, exactly reused EB theta=0 | `1.97e-14` | `1.42e-10` |

Каждая строка охватывает `181 x 7 = 1267` принятых корней соответствующей
модели.

## Независимая straight-reference проверка

Для `beta=0`, равных плеч, равных сечений и `theta=0` coupled section-clamp
спектр сравнивается с независимо собранным однородным стержнем длины `2L`:

```text
C_section = [[1,0,0,0,0,0],
             [0,1,0,0,0,0],
             [0,0,1,0,0,0]]
D_straight,section = C_section T_2L B_section
```

В independent reference `J_book` не используется. Для первых семи корней
получено:

- BookMaterial, Figure 2: maximum relative residual
  `4.597081368430e-12`;
- HMS/DX-209, Figures 3–4: maximum relative residual
  `6.399356074405e-11`.

Оба значения ниже hard criterion `1e-8`.

## Оформление

Figures 2–12 имеют оси `beta, degrees` и `Lambda`, major ticks
`0,15,...,90`. Цвет однозначно кодирует sorted position 1–6 и одинаков на
трёх рисунках. Для counterpart curves используется тот же цвет: solid
linewidth `1.6`, dashed linewidth `1.3`. Figures 3–4 имеют одинаковые размер
`7.2 x 4.8 in`, `xlim`, ticks и общий `ylim = [1.6771721625,
7.3506453369]`. Figure 2 имеет собственный `ylim = [1.7005200703,
9.0889311315]`. Figure 1 имеет размер `11 x 8.2 in`.

Ни на одном axes и ни на одном figure нет legend object. Заголовки внутри
рисунков отсутствуют. PNG сохраняются при 300 dpi на белом фоне, PDF — в
vector format.

Figures 5–6 имеют общий `ylim = [1.6459362339, 8.2269789489]`, Figures 7–8 —
общий `ylim = [1.2085724051, 7.8080869689]`; читаемость не потребовала
раздельных limits. Figures 9–12 используют строго тот же Figure-7 `ylim =
[1.2085724051, 7.8080869689]`. Все новые figures имеют размер `7.2 x 4.8 in`.

## Численные diagnostics

Различия ниже являются diagnostics, а не acceptance thresholds и не задают
ранжирование физической правильности моделей.

| Figure | maximum relative difference | beta, deg | mode |
|---|---:|---:|---:|
| 2, section clamp vs book slope | `6.735357686257e-03` | 0 | 6 |
| 3, Timoshenko vs EB | `1.063085751210e-02` | 0 | 6 |
| 4, Timoshenko vs EB | `1.416300859176e-02` | 0 | 6 |
| 5, unequal-length Timoshenko vs EB | `1.051903119037e-02` | 0 | 6 |
| 6, unequal-length/thickness Timoshenko vs EB | `1.051083268607e-02` | 0 | 5 |
| 7, theta=5 full model vs theta=0 EB approximation | `9.442988425292e-02` | 0 | 2 |
| 8, Chapter-2 theta=15 vs theta=0 | `3.306307663917e-01` | 0 | 1 |
| 9, theta=1 vs theta=0 EB baseline | `1.100143905580e-02` | 0 | 4 |
| 10, theta=2 vs theta=0 EB baseline | `2.481075259751e-02` | 0 | 4 |
| 11, theta=3 vs theta=0 EB baseline | `4.556772791475e-02` | 0 | 4 |
| 12, theta=4 vs theta=0 EB baseline | `6.673192621912e-02` | 80 | 4 |

Per-mode maxima для Figures 2, 3 и 4 полностью сохранены в `report.md` и
`plot_manifest.json`. Scientific/root runtime полного расчёта составил
`860.7748033 s`. Matrix evaluation counts также сохранены в manifest, отдельно
для scaled и physical-raw quality branches и straight references.

Минимальные relative neighbor gaps новых figures: Figure 5 —
`6.953517569581e-03` при `beta=0`, pair 5–6, Timoshenko; Figure 6 —
`6.090494268655e-03` при `beta=0`, pair 5–6, EB; Figure 7 —
`2.435600145749e-04` при `beta=79.5`, pair 4–5, reused EB; Figure 8 —
`1.747369722404e-04` при `beta=75.5`, pair 4–5, reused theta=0 Timoshenko.

Scientific runtimes Figures 5–8 составили соответственно `55.2674584`,
`52.4492923`, `20.6939058` и `21.7066812 s` (`150.1173377 s` суммарно).
Global fallback count новых families равен нулю. Fast transfer `expm` counts
по figures: `234556`, `228918`, `58542`, `63064`; соответствующие cache hit
rates: `0.511142`, `0.517218`, `0.755790`, `0.746487`. Обязательные семь
global anchors выполнены для каждой из шести новых calculated families.
Per-mode maxima и все отдельные performance counters сохранены в manifest и
`report.md`.

Для Figures 9–12 minimum relative neighboring gaps solid spectra равны:
`1.165407892956e-02` (theta=1, beta=0, pair 5–6),
`4.870138623548e-02` (theta=2, beta=64, pair 4–5),
`5.513255616911e-02` (theta=3, beta=0, pair 3–4) и
`4.211246777191e-02` (theta=4, beta=0, pair 3–4). Scientific runtimes:
`25.0078239`, `24.8618400`, `24.7267013`, `24.8123398 s`; global fallback
count равен нулю для каждого семейства. Fast transfer `expm` counts:
`58685`, `57867`, `57837`, `58198`.

При `beta=0` diagnostic energy labels показывают, что positions 5–6 меняют
порядок `torsion-like/bending-like` между theta `0` и `1 deg`, а positions
3–4 — между theta `3` и `4 deg`. Полные семь частот и fractions для каждого
малого угла находятся в `theta_small_beta0_mode_character.csv` и `report.md`.
Это не mode tracking и не утверждение о modal descendants.

## Выходы и режимы воспроизведения

Создаются двенадцать PDF, двенадцать PNG, двенадцать CSV,
`theta_small_beta0_mode_character.csv`, `plot_manifest.json` и `report.md`.
CSV Figure 1 содержит 637 строк. Каждый CSV Figures 2–12 содержит
1267 строк, включая root 7 как guard. Обычный import не выполняет расчётов и
не создаёт figure objects.

Новые outputs:

```text
figure_05_timoshenko_vs_eb_unequal_lengths_book_slope.{pdf,png}
figure_05_data.csv
figure_06_timoshenko_vs_eb_unequal_lengths_and_thickness_book_slope.{pdf,png}
figure_06_data.csv
figure_07_monoclinic_theta5_vs_orthotropic_eb_approximation.{pdf,png}
figure_07_data.csv
figure_08_chapter2_theta15_vs_theta0.{pdf,png}
figure_08_data.csv
figure_09_monoclinic_theta1_vs_orthotropic_eb_approximation.{pdf,png}
figure_09_data.csv
figure_10_monoclinic_theta2_vs_orthotropic_eb_approximation.{pdf,png}
figure_10_data.csv
figure_11_monoclinic_theta3_vs_orthotropic_eb_approximation.{pdf,png}
figure_11_data.csv
figure_12_monoclinic_theta4_vs_orthotropic_eb_approximation.{pdf,png}
figure_12_data.csv
theta_small_beta0_mode_character.csv
fast_solver_validation.csv
fast_solver_benchmark.json
```

`--reuse-data` запрещает scientific recomputation и строит только из
сохранённых CSV. `--force-recompute` явно пересчитывает данные; эти опции
mutually exclusive. `--figure new` означает Figures 5–8. Повторный
`--figure theta-small --reuse-data --jobs 1` строит Figures 9–12 только из
сохранённых CSV и не вызывает scientific solver. Повторный
`--figure all --reuse-data` занял `3.000039 s` и сохранил точные SHA-256 всех
новых CSV и PNG; byte identity PDF не является требованием из-за допустимых
backend metadata.

## Ограничения интерпретации

- Figure 1 — проверенное воспроизведение отдельного книжного single-rod
  результата, а не новый coupled-rod результат.
- Figures 2–8 показывают sorted spectral positions, а не физически
  отслеженные modal descendants.
- Различия clamps и theories не являются критерием выбора «правильной» модели.
- Workflow не включает MAC, shapes, FEM, 3D FEM, mesh refinement или широкое
  parameter study.
- Научные builders происходят только из
  `yartsev_ch2_monoclinic_rod.py`, `yartsev_ch2_coupled_rods.py` и
  `yartsev_ch2_rectangular_eb.py`.
- Модели, UT statuses, 1D validation phase и 3D-A0 status не изменяются.
- Figure 7 является намеренной диагностикой применения orthotropic theta=0 EB
  approximation к weakly monoclinic rod; её maximum difference не является
  чистой EB-vs-Timoshenko truncation error.
- Figure 8 сравнивает два угла внутри одной Chapter-2 theory и измеряет
  совокупный эффект rotated properties и `Sbar16` coupling.
