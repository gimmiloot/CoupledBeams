# Симметрично слоистая балка Reddy

Раздел содержит первую автономную модель одного прямого симметрично
слоистого стержня. Реализованы изгибная FSDT-подсистема, независимая
продольная подсистема и их блочно-диагональное объединение.

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
- `tests/test_reddy_symmetric_laminated_beam.py` — целевые регрессии;
- `tests/test_reddy_inplane_geometry.py` — coordinate-gate regressions;
- `tests/test_reddy_symmetric_coupled_beams_beta0.py` — joint,
  virtual-work, transfer, root-inventory и direct-reference regressions;
- `tests/test_reddy_symmetric_coupled_beams_ritz.py` — basis, quadrature,
  constrained matrices, beta=0 bridge и natural-equilibrium regressions;
- `tests/test_reddy_four_ply_isotropic_limit.py` — constitutive, section,
  legacy-adapter, local-space, frozen-evidence и closing regressions;
- `tests/data/reddy_ch4_table_4_3_3.json` — машинная транскрипция источника.

Generated CSV, JSON, report и две диагностические фигуры находятся в
`results/laminated_beams/reddy_symmetric_single_beam/`. Каталог игнорируется
Git. Отдельные RLB-1 generated data находятся в
`results/laminated_beams/reddy_symmetric_coupled_beta0_pilot/`; figures
в этом каталоге не создаются. Frozen RLB-1C-ISO evidence и closing report
находятся в игнорируемом каталоге
`results/laminated_beams/reddy_four_ply_isotropic_limit_validation/`.

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
```

Последняя команда только проверяет уже замороженные inventories, добавляет
closing provenance и выполняет два разрешённых bracket-local refinements. Она
не запускает global root search, direct fixed--fixed solver или Ritz solve.

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
