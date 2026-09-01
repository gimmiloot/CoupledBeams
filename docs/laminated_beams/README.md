# Симметрично слоистая балка Reddy

Раздел содержит первую автономную модель одного прямого симметрично
слоистого стержня. Реализованы изгибная FSDT-подсистема, независимая
продольная подсистема и их блочно-диагональное объединение.

All frequency-versus-parameter maps in this direction inherit the project-wide
[`frequency-map-v1` policy](../numerics/frequency_map_computation_policy.md)
unless a local note records an explicit override.

Сохранённый исторический результат RLB-0:

```text
RLB-0-SOURCE: PASS_WITH_DOCUMENTED_SOURCE_RECONSTRUCTION
RLB-0A-BENDING: PARTIAL_PASS
RLB-0B-AXIAL: PASS
RLB-0C-COMBINED-UNION: PASS
OVERALL: PARTIAL_PASS
```

Отдельный предварительный gate следующего этапа:

```text
RLB-1G-COORDINATES: PASS
```

После координатного gate выполнен узкий прямой coupled-пилот:

```text
RLB-1J-JOINT-THEORY: PASS
RLB-1A-BETA0-HOMOGENEOUS: PASS
RLB-1B-BETA0-STEPPED: PASS
RLB-1-BETA0-ROOT-INVENTORY: PASS
OVERALL: PASS
```

Последующий независимый двухплечевой Ritz-gate остановлен в прямом пределе:

```text
RLB-1C0-BETA0-RITZ-BRIDGE: FAIL
RLB-1C-NONZERO-BETA-SPECTRUM: FAIL
RLB-1C-RITZ-NATURAL-EQUILIBRIUM: PARTIAL_PASS
RLB-1D-MODE-SHAPES: FAIL
RLB-1S-SYMMETRIES: FAIL
OVERALL: FAIL
```

При разрешённом guard `N=16` верхние позиции полиномиального Ritz inventory
не прошли пороги сходимости и transfer/Ritz для полного first-13 набора.
Поэтому воспроизводимый validation workflow не запускал спектральные расчёты
при `beta!=0`. Предварительный read-only feasibility audit успел выполнить
несохранённые Ritz-only пробы при `beta=30/90`; transfer search и generated
nonzero outputs отсутствуют. Этот статус не изменяет прошедший transfer
baseline RLB-1.

Отдельный четырёхслойный изотропный предел закрыт по двум независимым
coupled determinant inventories:

```text
RLB-ISO-4PLY-CONSTITUTIVE: PASS
RLB-ISO-SECTION-REDUCTION: PASS
RLB-ISO-LEGACY-RECTANGULAR-ADAPTER: PASS
RLB-ISO-LOCAL-ARM-EQUIVALENCE: PASS
RLB-ISO-COUPLED-SPECTRUM: PASS
RLB-ISO-MODE-SHAPES: PASS
SCIENTIFIC_OVERALL: PASS_WITH_AUXILIARY_NUMERICAL_QUALIFICATIONS
```

Все 104 сопоставления первых 12 корней и root 13 guard проходят порог
`1e-8`. Direct `beta=0` 3x3 diagnostic исключён из scientific status из-за
обусловленности его mixed-regime basis. Единственное вспомогательное
превышение legacy arm exchange проходит после локального уточнения двух roots
только внутри их собственных frozen brackets. Исторический Ritz FAIL остаётся
отдельным результатом и в этот overall не входит.

RLB-1G фиксирует физический трёхмерный базис двух плеч. RLB-1J выводит
матрицу идеального жёсткого узла из глобальных физических условий и
независимо сверяет её с замкнутой формой. RLB-1A и RLB-1B проверяют только
прямой предел `beta=0`: однородный fixed--fixed стержень,
материальный интерфейс и зеркальную перестановку материалов. Этот PASS не
является проверкой спектра углового узла.

`PARTIAL_PASS` для RLB-0A связан с 22 расхождениями с Table 4.3.3 за
пределами точности трёх напечатанных десятичных знаков. Аналитические
изгибные gates, Eq. (4.3.51), продольная модель и объединение спектров
прошли. Параметры, граничные условия и допуски не подбирались по отдельным
случаям.

## RLB-2E: карта контраста жёсткости

Конечная карта RLB-2E исследует только контраст жёсткости и положение
материалов по толщине четырёхслойных стержней при фиксированных
`mu=tau=0` и `beta=30 deg`:

```text
RLB-2E: PASS
```

Для трёх конфигураций рассчитаны 41/41 значений `chi`, первые восемь
независимо упорядоченных частот и root 9 guard. Конститутивная проверка,
контроль соседних значений и перестановка плеч прошли. Формы, MAC,
идентичность ветвей, veering и распределение энергии не определялись.
Полный контракт и ограничения приведены в
[`reddy_stiffness_layout_contrast_sweep.md`](reddy_stiffness_layout_contrast_sweep.md).

## RLB-2F: одноплечевой слоистый контраст

Этап RLB-2F строит одну конечную карту по знаковому параметру `xi` при
фиксированных `mu=tau=0` и `beta=30 deg`. Первое плечо имеет стопку
`M_outer/M_inner/M_inner/M_outer`, второе остаётся однородным
`M0/M0/M0/M0`. Отрицательные `xi` соответствуют более жёстким внутренним
слоям, положительные — более жёстким наружным слоям.

```text
RLB-2F: PASS
```

Карта содержит 81/81 значение `xi`, 729 строк `BASE`, независимо
упорядоченные позиции 1--8 и 81 root-9 guard. Все 13 отмеченных neighbour
audit точек воспроизведены локальными проверками; unresolved points нет.
Проверка перестановки плеч прошла. Полный контракт, численная диагностика и
ограничения приведены в
[`reddy_one_arm_layered_contrast_sweep.md`](reddy_one_arm_layered_contrast_sweep.md).

## RLB-2G: двойственность компоновок по массе

Этап RLB-2G строит два конечных density-аналога карт RLB-2E и RLB-2F при
фиксированных `mu=tau=0` и `beta=30 deg`. Упругие свойства, статические
жёсткости и погонная масса сохраняются. Компоновки различаются только
вращательной инерцией `J` вследствие распределения плотности по толщине.

```text
RLB-2G: PASS
```

Эксперимент A содержит 123/123 группы и 1107 строк `BASE`; эксперимент B —
81/81 группу и 729 строк `BASE`. Общая точка нулевого контраста рассчитана
один раз для четырёх логических групп. Все 204 root-9 guards прошли, neighbour
flags и unresolved points отсутствуют, а arm-swap diagnostic получил статус
`PASS`. Карты содержат только независимо упорядоченные позиции 1--8 и не
задают идентичность модальных ветвей. Полный контракт и ограничения приведены
в [`reddy_mass_layout_duality.md`](reddy_mass_layout_duality.md).

## RLB-2H: видимость продольной жёсткости

Этап RLB-2H строит две конечные карты при `beta=0,30 deg` по сетке
`alpha_A=A/A0=0.70:0.02:1.30`. Перераспределение in-plane жёсткости между
четырьмя 0°-слоями меняет только reduced axial stiffness `A`; значения
`D`, `S`, `m` и `J` сохраняются.

```text
RLB-2H: PASS
```

Обе карты содержат 31/31 группу, всего 558 строк `BASE` и 62 root-9 guards.
Все 14 neighbour flags воспроизведены локальными пересчётами; unresolved
points нет. Beta=0 reference использует завершённую coupled characteristic
group при `alpha_A=1` и exact axial frequencies. Direct fixed-fixed boundary
matrices проверялись на сингулярность на заданных union frequencies; эти
частоты не локализовались независимо. При `beta=30 deg` все восемь sorted
positions не убывают на принятой сетке, однако их относительные изменения
различаются. Полный контракт и ограничения приведены в
[`reddy_axial_stiffness_visibility.md`](reddy_axial_stiffness_visibility.md).

## RLB-2I: глобально эквивалентные шестислойные ламинаты

Этап RLB-2I задаёт 51 симметричный шестислойный ламинат с тремя
зеркальными парами 0°-слоёв. Точное направление множителей (2,-3,1)
сохраняет полные матрицы A и D; неизменными остаются также S, m и J.

```text
RLB-2I-EXACT-EQUIVALENCE: PASS
RLB-2I-LAYERWISE-RECONSTRUCTION: PASS
RLB-2I-BOUNDARY-MATRIX-EQUIVALENCE: PASS
RLB-2I-SPECTRAL-SPOT-CHECK: PASS
RLB-2I-FIGURES: PASS
OVERALL: PASS
```

Послойное восстановление показывает изменение долей энергии и переход
максимального изгибного напряжения между промежуточной и наружной парами при
zeta=-1/9, несмотря на постоянную глобальную изгибную жёсткость. Частотная
карта по zeta намеренно не строилась: выполнены только шесть spot checks при
zeta=-0.25,0,+0.25 и beta=0,30 deg, roots 1--9. Полный контракт и ограничения
приведены в
[`reddy_six_ply_equivalent_laminates.md`](reddy_six_ply_equivalent_laminates.md).

## RLB-2J: попарный перенос жёсткости в шестислойном ламинате

Этап RLB-2J сравнивает переносы CENTER--MIDDLE, MIDDLE--OUTER и
CENTER--OUTER в симметричном шестислойном 0°-ламинате. Все три направления
сохраняют reduced (A,S,m,J), но изменяют (D/D_0) с точными рычагами
1:2:3.

```text
RLB-2J: PASS
```

Все три `frequency-map-v1` карты содержат 81/81 точки, независимо
отсортированные positions 1--8 и root 9 как completeness guard. Все 25
neighbour flags воспроизведены локально, unresolved points отсутствуют.
Спектры в точно совпадающих (D/D_0) collapse до относительной разности
не более (1.94\cdot10^{-12}). Полный контракт и ограничения приведены в
[`reddy_six_ply_pairwise_stiffness_transfer.md`](reddy_six_ply_pairwise_stiffness_transfer.md).

## Соглашения

У Reddy координата по толщине положительна вниз, а слои нумеруются сверху
вниз. Python API хранит стопку снизу вверх в проектной координате
`z_project = -z_Reddy`. Различие явно записано в source contract. Оно не
меняет текущие палиндромные симметричные benchmark-стопки, но существенно
для будущей модели с \(B\ne0\).

Коэффициент сдвига `K` является обязательным входом. Benchmark Table 4.3.3
задаёт `K=5/6` с provenance
`INFERRED_FROM_INTERNAL_NUMERICAL_CONSISTENCY`; это значение не считается
универсальным свойством произвольного ламината.

## Файлы

- `reddy_ch4_source_contract.md` — формулы, страницы, нормировка и политика
  реконструкции источника;
- `reddy_symmetric_single_beam_validation.md` — численные результаты и
  ограничения;
- `reddy_inplane_coordinate_contract.md` — физический координатный контракт
  RLB-1G и сопоставление со старыми display/FEM соглашениями;
- `reddy_table_4_3_3_discrepancy_audit.md` — независимая
  Rayleigh--Ritz-проверка RLB-0A-S1 и классификация ограничений source table;
- `reddy_symmetric_rigid_joint.md` — инвариантный вывод, скалярные
  условия, полная матрица узла и virtual-work gate;
- `reddy_symmetric_coupled_beta0_validation.md` — homogeneous,
  stepped, reflection и root-inventory результаты RLB-1;
- `reddy_symmetric_coupled_nonzero_beta_validation.md` — независимая
  двухплечевая Ritz-постановка, результат обязательного beta=0 bridge и
  зафиксированная остановка до ненулевого угла;
- `reddy_four_ply_isotropic_limit_validation.md` — аналитический
  четырёхслойный изотропный предел, независимый rectangular Timoshenko
  comparator, frozen spectra, формы и auxiliary qualifications;
- `reddy_stiffness_layout_contrast_sweep.md` — конечная карта RLB-2E по
  контрасту жёсткости для трёх расположений четырёх 0°-слоёв;
- `reddy_one_arm_layered_contrast_sweep.md` — контракт и область RLB-2F для
  знакового контраста одного слоистого плеча и однородного reference-плеча;
- `reddy_mass_layout_duality.md` — два density-layout эксперимента RLB-2G,
  аналитические формулы для `J/J0` и конечные sorted-spectrum результаты;
- `reddy_axial_stiffness_visibility.md` — A-only контракт RLB-2H, beta=0
  reference и конечные sorted-spectrum результаты при `beta=0,30 deg`;
- `reddy_six_ply_equivalent_laminates.md` — точное семейство RLB-2I,
  послойные напряжения, доли энергии и ограниченная спектральная регрессия;
- `reddy_six_ply_pairwise_stiffness_transfer.md` — три finite-карты RLB-2J,
  точные изгибные рычаги 1:2:3 и matched-(D) collapse;
- `scripts/lib/reddy_symmetric_laminated_beam.py` — узкий вычислительный API;
- `scripts/lib/reddy_inplane_geometry.py` — изолированный physical-coordinate
  helper;
- `scripts/lib/reddy_symmetric_coupled_beams.py` — узкий coupled helper,
  переиспользующий single-beam и coordinate APIs;
- `scripts/lib/reddy_symmetric_coupled_beams_ritz.py` — независимая
  constrained Rayleigh--Ritz-модель, не импортирующая transfer joint helper;
- `scripts/lib/isotropic_rectangular_timoshenko_coupled_beams.py` — независимый
  closed-form comparator для прямоугольного сечения с circular-backcompat;
- `scripts/analysis/laminated_beams/validate_reddy_symmetric_single_beam.py`
  — воспроизводимый CLI;
- `scripts/analysis/laminated_beams/pilot_reddy_symmetric_coupled_beams_beta0.py`
  — seed-free диагностический pilot только для `beta=0`;
- `scripts/analysis/laminated_beams/validate_reddy_symmetric_coupled_beams_nonzero_beta.py`
  — RLB-1C CLI с обязательной остановкой после failed beta=0 bridge;
- `scripts/analysis/laminated_beams/sweep_reddy_stiffness_layout_contrast.py`
  — устойчивый RLB-2E entry point с `missing-only`, `resume` и `plot-only`;
- `scripts/analysis/laminated_beams/sweep_reddy_one_arm_layered_contrast.py`
  — устойчивый RLB-2F entry point с `missing-only` и `plot-only`;
- `scripts/analysis/laminated_beams/sweep_reddy_mass_layout_duality.py`
  — устойчивый RLB-2G entry point с `missing-only`, `plot-only` и
  `manifest-only`;
- `scripts/analysis/laminated_beams/sweep_reddy_axial_stiffness_visibility.py`
  — устойчивый RLB-2H entry point с `missing-only`, `plot-only` и
  `manifest-only`;
- `scripts/analysis/laminated_beams/analyze_reddy_six_ply_equivalent_laminates.py`
  — RLB-2I entry point с полным, `plot-only` и `manifest-only` режимами;
- `scripts/analysis/laminated_beams/sweep_reddy_six_ply_pairwise_stiffness_transfer.py`
  — RLB-2J entry point с `missing-only`, `plot-only` и `manifest-only`;
- `tests/test_reddy_symmetric_laminated_beam.py` — целевые регрессии;
- `tests/test_reddy_inplane_geometry.py` — coordinate-gate regressions;
- `tests/test_reddy_symmetric_coupled_beams_beta0.py` — joint,
  virtual-work, transfer, root-inventory и direct-reference regressions;
- `tests/test_reddy_symmetric_coupled_beams_ritz.py` — basis, quadrature,
  constrained matrices, beta=0 bridge и natural-equilibrium regressions;
- `tests/test_reddy_four_ply_isotropic_limit.py` — constitutive, section,
  legacy-adapter, local-space, frozen-evidence и closing regressions;
- `tests/test_reddy_stiffness_layout_contrast.py` — material, constitutive,
  policy, root-inventory, neighbour-audit и output regressions RLB-2E;
- `tests/test_reddy_one_arm_layered_contrast.py` — targeted contract,
  constitutive, inventory и output regressions RLB-2F;
- `tests/test_reddy_mass_layout_duality.py` — density, constitutive,
  shared-anchor, inventory, plot-only и output regressions RLB-2G;
- `tests/test_reddy_axial_stiffness_visibility.py` — A-only constitutive,
  beta=0 reference, inventory, plot-only и output regressions RLB-2H;
- `tests/test_reddy_six_ply_equivalent_laminates.py` — exact-equivalence,
  layerwise reconstruction, hotspot, spectral и zero-calculation regressions
  RLB-2I;
- `tests/test_reddy_six_ply_pairwise_stiffness_transfer.py` — geometry,
  constitutive, inventory, matched-(D), neighbour-audit и output regressions
  RLB-2J;
- `tests/data/reddy_ch4_table_4_3_3.json` — машинная транскрипция источника.

Generated CSV, JSON, report и две диагностические фигуры находятся в
`results/laminated_beams/reddy_symmetric_single_beam/`. Каталог игнорируется
Git. Отдельные RLB-1 generated data находятся в
`results/laminated_beams/reddy_symmetric_coupled_beta0_pilot/`; figures
в этом каталоге не создаются. Frozen RLB-1C-ISO evidence и closing report
находятся в игнорируемом каталоге
`results/laminated_beams/reddy_four_ply_isotropic_limit_validation/`.
CSV, manifest, report и трёхпанельный рисунок RLB-2E находятся в
`results/laminated_beams/reddy_stiffness_layout_contrast_sweep/`.
CSV, manifest, report и однопанельный рисунок RLB-2F находятся в
`results/laminated_beams/reddy_one_arm_layered_contrast_sweep/`.
Две частотные карты, диагностический рисунок `J/J0`, CSV, manifest и report
RLB-2G находятся в
`results/laminated_beams/reddy_mass_layout_duality/`.
Двухпанельная частотная карта, beta=0 reference figure, CSV, manifest и report
RLB-2H находятся в
`results/laminated_beams/reddy_axial_stiffness_visibility/`.
Четыре CSV, три диагностических рисунка, manifest и report RLB-2I находятся
в `results/laminated_beams/reddy_six_ply_equivalent_laminates/`.
Таблицы, два диагностических рисунка, manifest и report RLB-2J находятся в
`results/laminated_beams/reddy_six_ply_pairwise_stiffness_transfer/`.

## Запуск

Из корня репозитория:

```powershell
python scripts/analysis/laminated_beams/validate_reddy_symmetric_single_beam.py --source-check-only
python -m pytest -q -p no:cacheprovider tests/test_reddy_symmetric_laminated_beam.py
python scripts/analysis/laminated_beams/validate_reddy_symmetric_single_beam.py --full-validation
python scripts/analysis/laminated_beams/validate_reddy_symmetric_single_beam.py --plot-only
python scripts/analysis/laminated_beams/pilot_reddy_symmetric_coupled_beams_beta0.py --manifest-only
python scripts/analysis/laminated_beams/pilot_reddy_symmetric_coupled_beams_beta0.py
python -m pytest -q -p no:cacheprovider tests/test_reddy_symmetric_coupled_beams_beta0.py
python scripts/analysis/laminated_beams/validate_reddy_symmetric_coupled_beams_nonzero_beta.py --manifest-only
python scripts/analysis/laminated_beams/validate_reddy_symmetric_coupled_beams_nonzero_beta.py
python -m pytest -q -p no:cacheprovider tests/test_reddy_symmetric_coupled_beams_ritz.py
python -m pytest -q -p no:cacheprovider tests/test_reddy_four_ply_isotropic_limit.py
python scripts/analysis/laminated_beams/reddy_four_ply_isotropic_postprocess.py --close-existing-results
python scripts/analysis/laminated_beams/sweep_reddy_stiffness_layout_contrast.py --missing-only
python scripts/analysis/laminated_beams/sweep_reddy_stiffness_layout_contrast.py --plot-only
python -m pytest -q -p no:cacheprovider tests/test_reddy_stiffness_layout_contrast.py
python scripts/analysis/laminated_beams/sweep_reddy_one_arm_layered_contrast.py --missing-only
python scripts/analysis/laminated_beams/sweep_reddy_one_arm_layered_contrast.py --plot-only
python -m pytest -q -p no:cacheprovider tests/test_reddy_one_arm_layered_contrast.py
python scripts/analysis/laminated_beams/sweep_reddy_mass_layout_duality.py --missing-only
python scripts/analysis/laminated_beams/sweep_reddy_mass_layout_duality.py --plot-only
python scripts/analysis/laminated_beams/sweep_reddy_mass_layout_duality.py --manifest-only
python -m pytest -q -p no:cacheprovider tests/test_reddy_mass_layout_duality.py
python scripts/analysis/laminated_beams/sweep_reddy_axial_stiffness_visibility.py --missing-only
python scripts/analysis/laminated_beams/sweep_reddy_axial_stiffness_visibility.py --plot-only
python scripts/analysis/laminated_beams/sweep_reddy_axial_stiffness_visibility.py --manifest-only
python -m pytest -q -p no:cacheprovider tests/test_reddy_axial_stiffness_visibility.py
python scripts/analysis/laminated_beams/analyze_reddy_six_ply_equivalent_laminates.py
python scripts/analysis/laminated_beams/analyze_reddy_six_ply_equivalent_laminates.py --plot-only
python scripts/analysis/laminated_beams/analyze_reddy_six_ply_equivalent_laminates.py --manifest-only
python -m pytest -q -p no:cacheprovider tests/test_reddy_six_ply_equivalent_laminates.py
python scripts/analysis/laminated_beams/sweep_reddy_six_ply_pairwise_stiffness_transfer.py --missing-only
python scripts/analysis/laminated_beams/sweep_reddy_six_ply_pairwise_stiffness_transfer.py --plot-only
python scripts/analysis/laminated_beams/sweep_reddy_six_ply_pairwise_stiffness_transfer.py --manifest-only
python -m pytest -q -p no:cacheprovider tests/test_reddy_six_ply_pairwise_stiffness_transfer.py
```

Команда `reddy_four_ply_isotropic_postprocess.py --close-existing-results`
только проверяет уже замороженные inventories, добавляет closing provenance и
выполняет два разрешённых bracket-local refinements. Она не запускает global
root search, direct fixed--fixed solver или Ritz solve.

`--source-check-only` не импортирует научный solver. `--plot-only` читает
только сохранённые CSV/JSON и не выполняет matrix exponential или поиск
корней. Режим `--smoke` намеренно возвращает для неисполненного RLB-0C статус
`PARTIAL_PASS`, а внутренний статус — `NOT_RUN`.

## Границы этапа

Продольная подсистема является project derivation from symmetric CLT и
стандартной динамики стержня, а не результатом, напечатанным в §4.3.4.
Combined state не содержит искусственной продольно-изгибной связи. На этом
этапе не реализованы \(B\ne0\), torsion, damping, FEM и complex roots.
RLB-1 реализует идеальный точечный узел и прошедший transfer spectrum при
`beta=0`. В RLB-1C реализован независимый двухплечий Ritz solver, однако его
полный first-13 beta=0 bridge не прошёл при максимальном разрешённом порядке.
Официальные spectra/shapes при `beta!=0`, parameter sweep и branch tracking
поэтому не вычислялись; отдельная предварительная Ritz-only проба отмечена в
validation note как процедурное отклонение.
